# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: NPDC-NPOLAR
# Fernando Rayón Viña
# last update: 20260402T00:00 AST
#
# Norwegian Polar Institute — Marine Plankton Biodiversity
# https://npolar.no/
# https://data.npolar.no/plankton/
# - Events:      https://data.npolar.no/plankton/event
# - Expeditions: https://data.npolar.no/plankton/expedition
# - Occurrences: https://data.npolar.no/plankton/occurrence
# Repository: https://www.re3data.org/repository/r3d100012291

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(lubridate)
library(janitor)
library(DBI)
library(RSQLite)
library(this.path)
library(httr2)
library(fs)
library(ncdf4)

# LOAD PATHS ────────────────────────────────────────────────────────────────

script_dir <- this.dir()

setwd(script_dir)

raw_dir <- file.path(script_dir, "raw")

# IMPORT DATA ───────────────────────────────────────────────────────────────

events_raw <- read_csv(
  file.path(raw_dir, "npdc-plankton-events.csv"),
  show_col_types = FALSE
) |>
  clean_names()

occurrences_raw <- read_csv(
  file.path(raw_dir, "npdc-plankton-occurrences.csv"),
  show_col_types = FALSE
) |>
  clean_names()

# INTEGRATE DATA ────────────────────────────────────────────────────────────

events_supplement <- events_raw |>
  select(event_id, all_of(setdiff(names(events_raw), names(occurrences_raw))))

plankton <- occurrences_raw |>
  left_join(events_supplement, by = "event_id")

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data <- plankton |>
  filter(
    sample_type %in%
      c(
        "Phytoplankton taxonomy",
        "Microplankton taxonomy",
        "Ice algae taxonomy"
      ),
    gear_type != "Sediment trap",
    organism_quantity_type == "cells/l" | is.na(organism_quantity_type),
    sample_status != "missing" | is.na(sample_status),
    !is.na(organism_quantity),
    !(is.na(minimum_depth_in_meters) &
      is.na(maximum_depth_in_meters) &
      is.na(sea_ice_core_minimum_depth) &
      is.na(sea_ice_core_maximum_depth))
  ) |>
  mutate(
    date = as_date(event_date),
    year = year(date),
    month = month(date),
    day = day(date),
    taxon = paste(
      scientific_name,
      if_else(!is.na(identification_qualifier), identification_qualifier, "")
    ) |>
      str_squish() |>
      str_to_lower() |>
      str_replace_all("\\.", "") |>
      str_squish(),
    depth = case_when(
      sample_type == "Ice algae taxonomy" &
        !is.na(sea_ice_core_minimum_depth) &
        !is.na(sea_ice_core_maximum_depth) ~
        (sea_ice_core_minimum_depth + sea_ice_core_maximum_depth) / 2,
      sample_type == "Ice algae taxonomy" &
        !is.na(sea_ice_core_maximum_depth) ~ sea_ice_core_maximum_depth,
      sample_type == "Ice algae taxonomy" &
        !is.na(sea_ice_core_minimum_depth) ~ sea_ice_core_minimum_depth,
      !is.na(minimum_depth_in_meters) & !is.na(maximum_depth_in_meters) ~
        (minimum_depth_in_meters + maximum_depth_in_meters) / 2,
      !is.na(maximum_depth_in_meters) ~ maximum_depth_in_meters,
      !is.na(minimum_depth_in_meters) ~ minimum_depth_in_meters,
      TRUE ~ NA_real_
    ),
    depth_label = case_when(
      sample_type == "Ice algae taxonomy" &
        !is.na(sea_ice_core_minimum_depth) &
        !is.na(sea_ice_core_maximum_depth) &
        sea_ice_core_minimum_depth != sea_ice_core_maximum_depth ~ "range",
      !is.na(minimum_depth_in_meters) &
        !is.na(maximum_depth_in_meters) &
        minimum_depth_in_meters != maximum_depth_in_meters ~ "range",
      TRUE ~ "exact"
    ),
    depth_label_min = case_when(
      depth_label == "range" & sample_type == "Ice algae taxonomy" ~
        sea_ice_core_minimum_depth,
      depth_label == "range" ~ minimum_depth_in_meters,
      TRUE ~ NA_real_
    ),
    depth_label_max = case_when(
      depth_label == "range" & sample_type == "Ice algae taxonomy" ~
        sea_ice_core_maximum_depth,
      depth_label == "range" ~ maximum_depth_in_meters,
      TRUE ~ NA_real_
    ),
    internal_id = if_else(!is.na(field_number), field_number, occurrence_id),
    source = "npdc-npolar",
    subsource = paste(
      case_when(
        sample_type == "Phytoplankton taxonomy" ~ "phyto",
        sample_type == "Ice algae taxonomy" ~ "ice",
        sample_type == "Microplankton taxonomy" ~ "micro"
      ),
      str_replace_all(expedition, "[^a-zA-Z0-9]", "_"),
      sep = "_"
    ),
    s_count_method = "utermohl",
    s_enumeration_method = "light microscopy",
    s_sampling_method = gear_type,
    s_concentration_method = "sedimentation",
    s_microscopy_method = "light microscope",
    s_light_microscope = "inverted",
    s_preservation = fixative,
    s_counting_chamber = "utermohl"
  ) |>
  rename(
    latitude = decimal_latitude,
    longitude = decimal_longitude,
    abundance = organism_quantity
  ) |>
  select(
    internal_id,
    source,
    subsource,
    year,
    month,
    day,
    latitude,
    longitude,
    depth,
    depth_label,
    depth_label_min,
    depth_label_max,
    taxon,
    abundance,
    s_count_method,
    s_enumeration_method,
    s_sampling_method,
    s_concentration_method,
    s_microscopy_method,
    s_light_microscope,
    s_preservation,
    s_counting_chamber
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
