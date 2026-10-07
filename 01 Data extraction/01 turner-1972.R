# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Turner-1972
# Fernando Rayón Viña
# last update: 20260416T15:50 ADT
#
# Turner, J. T. (1972). The Phytoplankton of the Tampa Bay System, Florida.
# M.A. Thesis, Marine Science Institute, University of South Florida, St. Petersburg.
# Associated article: Turner, J. T. & Hopkins, T. L. (1974). Phytoplankton of
# the Tampa Bay System. Bull. Mar. Sci. 24(1):101-119.

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(readxl)
library(lubridate)
library(janitor)
library(DBI)
library(RSQLite)
library(this.path)

# LOAD PATHS ────────────────────────────────────────────────────────────────

script_dir <- this.dir()

setwd(script_dir)

fn <- file.path(script_dir, "original-data-Turner-1972.xlsx")

# IMPORT DATA ───────────────────────────────────────────────────────────────

d1 <- read_excel(fn, "station-environment-data")

d2 <- read_excel(fn, "count-data", guess_max = 5000)

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

station_metadata <- d1 |>
  select(
    year,
    month,
    day,
    season,
    latitude,
    longitude,
    depth_m,
    station_id,
    nitrate_umol_L,
    ammonia_umol_L,
    phosphate_umol_L,
    silicate_umol_L,
    temperature_C,
    salinity,
    site_impacted
  ) |>
  mutate(
    station_id = as.character(station_id),
    site_impacted = ifelse(site_impacted == "Y", TRUE, FALSE)
  ) |>
  clean_names()

# COUNT DATA ────────────────────────────────────────────────────────────────

date_station_taxa <- d2 |>
  select(!contains("biomass")) |>
  pivot_longer(
    cols = starts_with("station"),
    names_to = "station_id",
    values_to = "cells"
  ) |>
  mutate(
    abundance = cells * 1000000,
    station_id = str_remove(station_id, "^station-"),
    station_id = str_extract(station_id, "^[0-9]+")
  ) |>
  filter(!is.na(abundance), !is.na(station_id)) |>
  select(season, year, taxon, station_id, abundance) |>
  group_by(season, year, station_id, taxon) |>
  summarise(abundance = sum(abundance), .groups = "drop") |>
  pivot_wider(
    names_from = taxon,
    values_from = abundance,
    values_fill = 0
  ) |>
  pivot_longer(
    cols = -(season:station_id),
    names_to = "taxon",
    values_to = "abundance"
  )

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data <- date_station_taxa |>
  left_join(
    station_metadata |>
      select(year, month, day, season, latitude, longitude, station_id),
    by = c("year", "season", "station_id")
  ) |>
  filter(!is.na(latitude), !is.na(longitude)) |>
  mutate(
    date = make_date(year, month, day),
    datetime = format(
      as.POSIXct(date, tz = "UTC"),
      "%Y-%m-%dT00:00:00Z"
    ),
    taxon = taxon |>
      str_to_lower() |>
      str_replace_all("\\.", "") |>
      str_squish(),
    internal_id = paste(
      year,
      season,
      station_id,
      taxon,
      sep = "-"
    ) |>
      str_replace_all(" ", "_"),
    source = "turner-1972",
    depth = 0,
    depth_label = "exact",
    s_count_method = "microscopy",
    s_enumeration_method = "light microscopy",
    s_sampling_method = "surface water bottle",
    s_concentration_method = "settling",
    s_microscopy_method = "light microscope",
    s_preservation = "formalin-seawater 5%",
    s_volume_sampled = 0.3
  ) |>
  select(
    internal_id,
    source,
    year,
    month,
    day,
    datetime,
    latitude,
    longitude,
    depth,
    depth_label,
    taxon,
    abundance,
    s_count_method,
    s_enumeration_method,
    s_sampling_method,
    s_concentration_method,
    s_volume_sampled,
    s_microscopy_method,
    s_preservation
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
