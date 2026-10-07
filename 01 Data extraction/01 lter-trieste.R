# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: LTER_Trieste
# Original author: Sadra Dehghani
# Fernando Rayón Viña
# last update: 20260312T15:00 AST
#
# Phytoplankton North Adriatic - Gulf of Trieste C1 LTER
# Cabrini M., Fonda-Umani S. (2021) - OGS, Italy

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(janitor)
library(this.path)
library(RSQLite)
library(DBI)
library(lubridate)

# LOAD PATHS ────────────────────────────────────────────────────────────────

script_dir <- this.dir()

fn_event <- file.path(script_dir, "event.txt")

fn_occ <- file.path(script_dir, "occurrence.txt")

fn_emof <- file.path(script_dir, "extendedmeasurementorfact.txt")

# EVENT TABLE ───────────────────────────────────────────────────────────────

d_event <- read_delim(fn_event, delim = "\t", show_col_types = FALSE) |>
  clean_names() |>
  filter(type == "sample") |>
  mutate(
    date = ymd(event_date),
    year = year(date),
    month = month(date),
    day = day(date),
    depth_m = case_when(
      !is.na(minimum_depth_in_meters) & !is.na(maximum_depth_in_meters) ~
        (minimum_depth_in_meters + maximum_depth_in_meters) / 2,
      !is.na(minimum_depth_in_meters) ~ minimum_depth_in_meters,
      !is.na(maximum_depth_in_meters) ~ maximum_depth_in_meters,
      TRUE ~ NA_real_
    ),
    latitude = 45.6976666,
    longitude = 13.7083333
  ) |>
  select(event_id, year, month, day, depth_m, latitude, longitude)

# OCCURRENCE TABLE ──────────────────────────────────────────────────────────

d_occ <- read_delim(fn_occ, delim = "\t", show_col_types = FALSE) |>
  clean_names() |>
  rename(
    taxon_original = scientific_name,
    worms_lsid = scientific_name_id,
    occurrence_status = occurrence_status
  ) |>
  select(id, occurrence_id, taxon_original, worms_lsid, occurrence_status)

# EXTENDED MEASUREMENT OR FACT TABLE ────────────────────────────────────────

d_emof <- read_delim(fn_emof, delim = "\t", show_col_types = FALSE) |>
  clean_names() |>
  filter(id != "root") |>
  rename(abundance_l = measurement_value) |>
  filter(
    measurement_type_id ==
      "http://vocab.nerc.ac.uk/collection/P01/current/SDBIOL01/"
  ) |>
  mutate(abundance_l = suppressWarnings(as.numeric(abundance_l))) |>
  filter(!is.na(abundance_l)) |>
  select(id, occurrence_id, abundance_l)

# MERGE AND TRANSFORM ───────────────────────────────────────────────────────

count_data <- d_occ |>
  inner_join(d_emof, by = c("id", "occurrence_id")) |>
  inner_join(d_event, by = c("id" = "event_id")) |>
  mutate(
    taxon = taxon_original |>
      str_to_lower() |>
      str_replace_all("\\.", "") |>
      str_squish(),
    source = "lter-trieste",
    s_count_method = "utermohl",
    s_enumeration_method = "light microscopy",
    s_sampling_method = "niskin bottles",
    s_concentration_method = "sedimentation",
    s_microscopy_method = "light microscope",
    s_light_microscope = "inverted, phase contrast",
    s_preservation = "formaldehyde, Ca(HCO3)2-buffered (0.8%)",
    s_magnification = "200-320-400",
    s_counting_chamber = "utermohl"
  ) |>
  select(
    id,
    latitude,
    longitude,
    depth_m,
    year,
    month,
    day,
    taxon,
    abundance_l,
    source,
    s_count_method,
    s_enumeration_method,
    s_sampling_method,
    s_concentration_method,
    s_microscopy_method,
    s_light_microscope,
    s_preservation,
    s_magnification,
    s_counting_chamber
  ) |>
  distinct()

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
