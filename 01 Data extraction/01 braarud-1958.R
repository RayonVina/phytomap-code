# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Braarud-1958
# Original author: Julianne Jager
# Fernando Rayón Viña
# Original version: 2022-06-21 Julianne Jager

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(readxl)
library(here)
library(lubridate)
library(janitor)
library(stringr)
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

directory <- hereHelper("Braarud-1958")

fn <- here(directory, "original-data-Braarud-1958.xlsx")

sheet_names <- excel_sheets(fn)

# IMPORT DATA ───────────────────────────────────────────────────────────────

d1 <- read_excel(fn, sheet_names[4])

d2 <- read_excel(fn, sheet_names[5])

d3 <- read_excel(fn, sheet_names[6])

station_metadata <- d1 |>
  mutate(
    date = ymd(date),
    year = year(date),
    month = month(date),
    day = day(date)
  ) |>
  rename(salinity = `Salinity_‰`, temperature_C = Temp_C) |>
  clean_names() |>
  select(
    year,
    month,
    day,
    latitude,
    longitude,
    depth_m,
    station_id,
    temperature_c,
    salinity
  )

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data1 <- d2 |>
  pivot_longer(
    starts_with("date_"),
    names_to = "date",
    values_to = "abundance_l"
  ) |>
  mutate(
    date = str_remove(date, "date_"),
    date = str_remove(date, "\\.\\.\\.[1-9][0-9]?$")
  ) |> # repeated dates -- how do we distinguish them?
  mutate(
    date = mdy(date),
    year = year(date),
    month = month(date),
    day = day(date)
  ) |>
  clean_names() |>
  left_join(station_metadata) |>
  select(
    year,
    month,
    day,
    latitude,
    longitude,
    depth_m,
    station_id,
    taxon,
    abundance_l
  ) |>
  filter(!is.na(abundance_l)) |>
  clean_names()

count_data2 <- d3 |>
  pivot_longer(
    ends_with("_m"),
    names_to = "depth_m",
    values_to = "abundance_l"
  ) |>
  mutate(
    date = ymd(date),
    year = year(date),
    month = month(date),
    day = day(date),
    depth_m = str_remove(depth_m, "_m"),
    depth_m = as.numeric(depth_m)
  ) |>

  mutate(
    presence = ifelse(abundance_l == "+", TRUE, FALSE),
    abundance_l = as.numeric(abundance_l)
  ) |>
  left_join(station_metadata) |>
  select(
    year,
    month,
    day,
    latitude,
    longitude,
    depth_m,
    station_id,
    taxon,
    abundance_l
  ) |>
  filter(!is.na(abundance_l)) |>
  clean_names()

count_data <- bind_rows(count_data1, count_data2) |>
  distinct()

# STATION METADATA ──────────────────────────────────────────────────────────

station_coords <- station_metadata |>
  distinct(station_id, depth_m, latitude, longitude)

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data1 <- d2 |>
  pivot_longer(
    starts_with("date_"),
    names_to = "date",
    values_to = "abundance_l"
  ) |>
  mutate(
    date = str_remove(date, "date_"),
    # drop suffixes like .1, .2, etc.
    date = str_remove(date, "\\.[0-9]+$"),
    date = mdy(date),
    year = year(date),
    month = month(date),
    day = day(date)
  ) |>
  clean_names() |>
  # join only on station + depth
  left_join(station_coords, by = c("station_id", "depth_m")) |>
  select(
    year,
    month,
    day,
    latitude,
    longitude,
    depth_m,
    station_id,
    taxon,
    abundance_l
  ) |>
  filter(!is.na(abundance_l))

count_data2 <- d3 |>
  pivot_longer(
    ends_with("_m"),
    names_to = "depth_m",
    values_to = "abundance_l"
  ) |>
  mutate(
    date = ymd(date),
    year = year(date),
    month = month(date),
    day = day(date),
    depth_m = str_remove(depth_m, "_m"),
    depth_m = as.numeric(depth_m),
    presence = ifelse(abundance_l == "+", TRUE, FALSE),
    abundance_l = as.numeric(abundance_l)
  ) |>
  clean_names() |>
  left_join(station_coords, by = c("station_id", "depth_m")) |>
  select(
    year,
    month,
    day,
    latitude,
    longitude,
    depth_m,
    station_id,
    taxon,
    abundance_l
  ) |>
  filter(!is.na(abundance_l))

count_data <- bind_rows(count_data1, count_data2) |>
  distinct()

count_data <- count_data |>
  mutate(
    latitude = case_when(
      station_id == "UTSIRA" ~ 59.32,
      station_id == "SOGNESJØEN" ~ 61.00,
      station_id == "SKROVA" ~ 68.00,
      station_id == "EGGUM" ~ 68.38,
      TRUE ~ latitude
    ),
    longitude = case_when(
      station_id == "UTSIRA" ~ 4.98,
      station_id == "SOGNESJØEN" ~ 4.83,
      station_id == "SKROVA" ~ 14.65,
      station_id == "EGGUM" ~ 13.63,
      TRUE ~ longitude
    )
  )

# SET PROVENANCE ────────────────────────────────────────────────────────────

count_data <- count_data |> mutate(source = "braarud-1958")

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
