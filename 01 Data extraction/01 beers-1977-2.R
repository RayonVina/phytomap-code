# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Beers-1977-2
# Original author: Catherine Savard
# Fernando Rayón Viña
# Original version: 2022-07-06 Catherine Savard

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

directory <- hereHelper("Beers-1977-2")

fn <- here(directory, "original-data-Beers-1977-2.xlsx")

sheet_names <- excel_sheets(fn)

# IMPORT DATA ───────────────────────────────────────────────────────────────

d1 <- read_excel(fn, sheet_names[4])

d2 <- read_excel(fn, sheet_names[5])

capitalize <- function(x) {
  i <- utf8ToInt(x)[1]
  if (i > 96) {
    substr(x, 1, 1) <- intToUtf8(i - 32)
  }
  x
}

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data <- d2 |>
  filter(!is.na(taxon)) |>
  mutate(station_ID = as.numeric(station_ID)) |>
  rename(abundance_L = abundance_cells_L) |>
  mutate(year = year(date), month = month(date), day = day(date)) |>
  rowwise() |>
  mutate(taxon = capitalize(taxon)) |>
  ungroup()

count_data <- count_data |>
  left_join(
    d1 |> select(-year, -month, -day) |> remove_empty(),
    by = "station_ID"
  ) |>
  select(
    year,
    month,
    day,
    depth_m,
    station_ID,
    latitude,
    longitude,
    taxon,
    abundance_L
  ) |>
  filter(!is.na(abundance_L))

count_data <- clean_names(count_data)

# SET PROVENANCE ────────────────────────────────────────────────────────────

count_data <- count_data |> mutate(source = "beers-1977-2")

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
