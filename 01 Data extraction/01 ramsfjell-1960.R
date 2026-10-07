# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Ramsfjell-1960
# Original author: Julianne Jager
# Fernando Rayón Viña
# Original version: 2022-05-26 Julianne Jager

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(readxl)
library(here)
library(lubridate)
library(janitor)
library(dplyr)
library(DBI)
library(RSQLite)

# LOAD PATHS ────────────────────────────────────────────────────────────────

hereHelper <- function(dataset) {
  paths <- list.dirs(here::here(), recursive = TRUE, full.names = TRUE)
  paths <- paths[basename(paths) == dataset]
  if (length(paths) != 1L) {
    stop(sprintf("Expected one input directory named %s.", dataset))
  }
  paths
}

directory <- hereHelper("Ramsfjell-1960")

fn <- here(directory, "original-data-Ramsfjell-1960.xlsx")

sheet_names <- excel_sheets(fn)

# IMPORT DATA ───────────────────────────────────────────────────────────────

d1 <- read_excel(fn, sheet_names[4])

d2 <- read_excel(fn, sheet_names[5])

d3 <- read_excel(fn, sheet_names[6])

station_metadata <- d1 |>
  mutate(date = as_date(date)) |>
  mutate(
    date = ymd(date),
    day = day(date),
    month = month(date),
    year = year(date)
  ) |>
  pivot_longer(
    starts_with("temperature"),
    names_to = "depth_m",
    values_to = "temperature_C"
  ) |>
  pivot_longer(
    starts_with("salinity"),
    names_to = "depth_m2",
    values_to = "salinity"
  ) |>
  mutate(depth_m = str_remove(depth_m, "temperature_")) |>
  mutate(depth_m2 = str_remove(depth_m2, "salinity_")) |>
  subset(depth_m2 == depth_m) |>
  select(
    station_id,
    latitude,
    longitude,
    year,
    month,
    day,
    depth_m,
    temperature_C,
    salinity
  )

station_metadata <- clean_names(station_metadata) |>
  mutate(station_id = as.character(station_id), depth_m = as.numeric(depth_m))

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data1 <- d2 |>
  mutate(across(everything(), as.character)) |>
  pivot_longer(
    starts_with("station"),
    names_to = "station_id",
    values_to = "abundance_l"
  ) |>
  filter(!is.na(abundance_l)) |>
  select(taxon, depth_m, station_id, abundance_l) |>
  mutate(presence = ifelse(abundance_l == 0, FALSE, TRUE)) |>
  mutate(
    abundance_l = as.numeric(abundance_l),
    depth_m = as.numeric(depth_m)
  ) |>
  mutate(
    station_id = str_remove(station_id, "station_"),
    station_id = str_remove(station_id, "\\.\\.\\.[1-9][0-9][0-9]?$")
  ) |>
  left_join(station_metadata) |>
  select(
    year,
    month,
    day,
    depth_m,
    station_id,
    latitude,
    longitude,
    taxon,
    abundance_l,
    presence
  )

count_data1 <- clean_names(count_data1)

count_data2 <- d3 |>
  mutate(across(everything(), as.character)) |>
  pivot_longer(
    starts_with("depth_"),
    names_to = "depth_m",
    values_to = "abundance_l"
  ) |>
  filter(!is.na(abundance_l)) |>
  mutate(presence = ifelse(abundance_l == 0, FALSE, TRUE)) |>
  mutate(
    abundance_l = as.numeric(abundance_l),
    depth_m = as.numeric(str_remove(depth_m, "depth_")),
    abundance_l = as.numeric(abundance_l)
  ) |>
  left_join(station_metadata) |>
  select(
    year,
    month,
    day,
    depth_m,
    station_id,
    latitude,
    longitude,
    taxon,
    abundance_l,
    presence
  )

count_data <- mutate(bind_rows(count_data1, count_data2)) |>
  distinct()

# SET PROVENANCE ────────────────────────────────────────────────────────────

count_data <- count_data |> mutate(source = "ramsfjell-1960")

# PREPARE COUNT DATA ───────────────────────────────────────────────────────────
count_data <- count_data |>
  filter(
    !is.na(year),
    !is.na(latitude),
    !is.na(longitude),
    !is.na(taxon),
    !is.na(abundance_l)
  ) |>
  mutate(
    datetime = sprintf("%04d-%02d-%02d", year, month, day),
    n_uncertain = 0L,
    n_ltr = 0L,
    n_flag_for_removal = 0L
  )


# PREPARE DATABASE RECORDS ─────────────────────────────────────────────────────
count_data <- count_data |> mutate(across(where(is.factor), as.character))

column_aliases <- c(
  taxon = "t_taxon",
  abundance_l = "abundance",
  abundance_L = "abundance",
  depth_m = "depth",
  datetime_utc = "datetime",
  count_method = "s_count_method",
  taxon_original = "t_taxon",
  taxon_clean = "t_taxon_clean",
  notes = "n_notes"
)

for (old_name in intersect(names(column_aliases), names(count_data))) {
  new_name <- column_aliases[[old_name]]
  if (new_name %in% names(count_data)) {
    old_value <- as.character(count_data[[old_name]])
    new_value <- as.character(count_data[[new_name]])
    conflict <- !is.na(old_value) & !is.na(new_value) & old_value != new_value
    if (any(conflict)) {
      stop(sprintf("Conflicting columns: %s and %s.", old_name, new_name))
    }
    missing_value <- is.na(count_data[[new_name]])
    count_data[[new_name]][missing_value] <- count_data[[old_name]][
      missing_value
    ]
    count_data[[old_name]] <- NULL
  } else {
    names(count_data)[names(count_data) == old_name] <- new_name
  }
}

protocol_fields <- c(
  "s_count_method",
  "s_enumeration_method",
  "s_sampling_method",
  "s_sampling_method_notes",
  "s_concentration_method",
  "s_concentration_method_notes",
  "s_volume_sampled",
  "s_microscopy_method",
  "s_microscopy_notes",
  "s_preservation",
  "s_volume_concentrated",
  "s_volume_counted",
  "s_magnification",
  "s_counting_chamber"
)
count_data <- count_data |>
  mutate(across(any_of(protocol_fields), as.character)) |>
  mutate(across(where(is.factor), as.character))

if ("datetime" %in% names(count_data)) {
  count_data$datetime <- as.character(count_data$datetime)
}
if (anyDuplicated(names(count_data))) {
  stop("Duplicate output column names.")
}
required_fields <- c("source", "t_taxon", "abundance")
if (!all(required_fields %in% names(count_data))) {
  stop("Output must contain source, t_taxon, and abundance.")
}

# EXPORT TO DATABASE ──────────────────────────────────────────────────────────
db_path <- Sys.getenv("PHYTOMAP_DB")

con <- dbConnect(RSQLite::SQLite(), db_path, flags = RSQLite::SQLITE_RW)
tryCatch(
  {
    if (!dbExistsTable(con, "abundance_raw")) {
      stop("The working database must already contain abundance_raw.")
    }
    dbWithTransaction(con, {
      new_fields <- setdiff(
        names(count_data),
        dbListFields(con, "abundance_raw")
      )
      for (field in new_fields) {
        field_type <- dbDataType(con, count_data[[field]])
        dbExecute(
          con,
          paste(
            "ALTER TABLE abundance_raw ADD COLUMN",
            dbQuoteIdentifier(con, field),
            field_type
          )
        )
      }
      dbWriteTable(
        con,
        name = "abundance_raw",
        value = count_data,
        append = TRUE,
        row.names = FALSE
      )
    })
  },
  finally = dbDisconnect(con)
)
rm(con)
