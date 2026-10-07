# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Kang-1992
# Fernando Rayón Viña
# last update: 20260416T00:00 ADT
#
# Kang, S.-H. (1992)
# Phytoplankton in the Antarctic Marginal Ice Zone.
# PhD thesis, Texas A&M University, December 1992.
# ProQuest: https://www.proquest.com/docview/304001198
# Data digitised from Appendix B (count-data-3): phytoplankton densities
# (cells·L⁻¹). Campaigns: AMERIEZ 86 (autumn 1986), AMERIEZ 88 (winter 1988),
# ODP Leg 119 (summer 1988), ICECOLORS 90 (spring 1990).
# Only full-cell counts are used as abundance; empty frustules discarded.
# Zero-filling: taxa absent at a station×depth are treated as true absences
# and reconstructed to 0 via double-pivot.

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

fn <- file.path(
  script_dir,
  "original-data-updated-Kang-Phd-thesis-1992.xlsx"
)

# IMPORT DATA ───────────────────────────────────────────────────────────────

station_raw <- read_excel(fn, sheet = "station-environment-data") |>
  clean_names() |>
  remove_empty(c("rows", "cols"))

count_raw <- read_excel(fn, sheet = "count-data-3", guess_max = 5000) |>
  clean_names() |>
  remove_empty(c("rows", "cols"))

station_metadata <- station_raw |>
  mutate(
    year = as.integer(year),
    month = as.integer(month),
    day = as.integer(day),
    date = make_date(year, month, day)
  )

# TRANSFORM: COUNT DATA ─────────────────────────────────────────────────────

key_cols <- c(
  "station_id",
  "transect",
  "station_number",
  "day",
  "month",
  "year",
  "depth_m",
  "taxon"
)

meta_cols <- c(
  "station_id",
  "transect",
  "station_number",
  "day",
  "month",
  "year",
  "depth_m"
)

count_summed <- count_raw |>
  group_by(across(all_of(key_cols))) |>
  summarise(
    abundance = sum(abundance_l_full, na.rm = TRUE),
    .groups = "drop"
  )

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data_raw <- count_summed |>
  pivot_wider(
    names_from = taxon,
    values_from = abundance,
    values_fill = 0
  ) |>
  pivot_longer(
    cols = -all_of(meta_cols),
    names_to = "taxon",
    values_to = "abundance"
  )

count_data <- count_data_raw |>
  left_join(
    station_metadata |>
      select(station_id, latitude, longitude) |>
      distinct(station_id, .keep_all = TRUE),
    by = "station_id"
  ) |>
  mutate(
    year = as.integer(year),
    month = as.integer(month),
    day_raw = day,
    day = if_else(is.na(day), 15L, as.integer(day)),
    date = make_date(year, month, day),
    datetime = format(
      as.POSIXct(paste(date, "00:00"), tz = "UTC"),
      "%Y-%m-%dT%H:%M:%SZ"
    ),
    taxon = taxon |>
      str_to_lower() |>
      str_replace_all("\\.", "") |>
      str_squish(),
    subsource = case_when(
      str_detect(station_id, "^AMERIEZ-83") ~ "ameriez-83",
      str_detect(station_id, "^AMERIEZ-86") ~ "ameriez-86",
      str_detect(station_id, "^AMERIEZ-88") ~ "ameriez-88",
      str_detect(station_id, "^ICECOLORS-90") ~ "icecolors-90",
      str_detect(station_id, "^ODP-Leg-119") ~ "odp-leg-119",
      TRUE ~ "kang-1992-other"
    ),
    internal_id = paste(station_id, year, month, depth_m, sep = "-"),
    source = "kang-1992",
    depth_label = "exact",
    s_count_method = "utermohl",
    s_enumeration_method = "light microscopy",
    s_sampling_method = "bottle",
    n_notes = if_else(
      is.na(day_raw),
      "Day set to 15 by convention (month and year only in source)",
      NA_character_
    ),
    s_concentration_method = case_when(
      str_detect(subsource, "icecolors") ~ "filtration (HPMA)",
      str_detect(subsource, "odp") ~ "utermohl and filtration (HPMA)",
      TRUE ~ "utermohl and filtration (HPMA)"
    ),
    s_microscopy_method = "light microscope",
    s_light_microscope = "inverted",
    s_preservation = case_when(
      str_detect(subsource, "ameriez-88") ~ "1% glutaraldehyde",
      TRUE ~ "1% buffered formalin (hexamine)"
    )
  ) |>
  select(
    internal_id,
    source,
    subsource,
    datetime,
    year,
    month,
    day,
    latitude,
    longitude,
    depth = depth_m,
    depth_label,
    taxon,
    abundance,
    n_notes,
    s_count_method,
    s_enumeration_method,
    s_sampling_method,
    s_concentration_method,
    s_microscopy_method,
    s_light_microscope,
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
