# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Beers-et-al-1977
# Original author: Julianne Jager
# Fernando Rayón Viña
# Original version: 2022-07-06 Julianne Jager

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(readxl)
library(here)
library(lubridate)
library(janitor)
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

directory <- hereHelper("Beers-et-al-1977")

fn <- here(directory, "original-data-Beers-et-al-1977.xlsx")

# IMPORT DATA ───────────────────────────────────────────────────────────────

d1 <- read_excel(fn, "station-environment-data") |> clean_names()

d2 <- read_excel(fn, "count-data") |> clean_names()

station_metadata <- d1

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data <- d2 |>
  filter(!is.na(count)) |>
  left_join(
    station_metadata,
    by = c("expedition", "station_id", "day", "month", "year")
  ) |>
  select(
    expedition,
    year,
    month,
    day,
    depth_m,
    station_id,
    latitude,
    longitude,
    taxon,
    abundance_l
  )

reshape_one <- function(count_data) {
  count_data |>
    pivot_wider(
      values_from = abundance_l,
      names_from = taxon,
      values_fill = 0,
      names_prefix = "taxon_"
    ) |>
    pivot_longer(
      starts_with("taxon_"),
      names_to = "taxon",
      names_prefix = "taxon_",
      values_to = "abundance_l"
    )
}

bind_rows(
  count_data |> filter(expedition == "Cato, leg 1") |> reshape_one(),
  count_data |> filter(expedition == "Climax VII") |> reshape_one(),
  count_data |> filter(expedition == "Dramamine II") |> reshape_one(),
  count_data |> filter(expedition == "South tow, leg 13") |> reshape_one(),
  count_data |> filter(expedition == "Tasaday, leg 11") |> reshape_one()
) -> count_data

# SET PROVENANCE ────────────────────────────────────────────────────────────

count_data <- count_data |> mutate(source = "beers-et-al-1977")

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
