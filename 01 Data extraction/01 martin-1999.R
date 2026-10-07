# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Martin-1999
# Fernando Rayón Viña
# last update: 20260416T10:00AST
#
# Martin, J.L., LeGresley, M.M., Strain, P.M., Clement, P. (1999)
# Phytoplankton Monitoring in the Southwest Bay of Fundy During 1993-96.
# Can. Tech. Rep. Fish. Aquat. Sci. 2265: iv + 132 p.
# Fisheries and Oceans Canada — Biological Station, St. Andrews, NB

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
  "original-data-Martin-LeGresley-Strain-Clement-1999.xlsx"
)

# IMPORT DATA ───────────────────────────────────────────────────────────────

station_raw_aux_1 <- read_excel(fn, sheet = "station-environment-data_1") |>
  clean_names() |>
  remove_empty(c("rows", "cols"))

station_raw_aux_2 <- read_excel(fn, sheet = "station-environment-data_2") |>
  clean_names() |>
  remove_empty(c("rows", "cols"))

count_raw <- read_excel(fn, sheet = "count-data", guess_max = 5000) |>
  clean_names() |>
  remove_empty(c("rows", "cols"))

# STATION COORDINATES ───────────────────────────────────────────────────────

station_coords <- bind_rows(
  station_raw_aux_1 |> select(station_id, latitude, longitude) |> distinct(),
  station_raw_aux_2 |> select(station_id, latitude, longitude) |> distinct()
) |>
  mutate(station_id = as.integer(station_id)) |>
  distinct(station_id, .keep_all = TRUE)

# TRANSFORM: COUNT DATA ─────────────────────────────────────────────────────

station_cols <- setdiff(colnames(count_raw), c("date", "taxon"))

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data <- count_raw |>
  # Preserve NAs to distinguish true absences from structural zeros.
  # NAs = taxon not assessed on that date × station × depth combination.
  # Zeros = taxon assessed and not observed (confirmed absence).
  pivot_longer(
    cols = all_of(station_cols),
    names_to = "col_key",
    values_to = "abundance"
  ) |>
  # Fill structural zeros: taxon × date × station × depth combos not in source
  # exist as real absences and must be retained as 0.
  pivot_wider(
    names_from = col_key,
    values_from = abundance,
    values_fill = 0
  ) |>
  pivot_longer(
    cols = all_of(station_cols),
    names_to = "col_key",
    values_to = "abundance"
  ) |>
  filter(!is.na(abundance)) |>
  # col_key format after clean_names(): "station_3_0m", "station_16_50m"
  # Remove "station_" prefix, then split at last "_" before depth digits.
  mutate(col_key = str_remove(col_key, "^station_")) |>
  separate(
    col_key,
    into = c("station_id", "depth_str"),
    sep = "_(?=[0-9]+m?$)"
  ) |>
  mutate(
    station_id = as.integer(station_id),
    depth = as.numeric(str_remove(depth_str, "m$")),
    date = as.Date(date),
    year = year(date),
    month = month(date),
    day = day(date),
    datetime = format(
      as.POSIXct(paste(date, "00:00"), tz = "UTC"),
      "%Y-%m-%dT%H:%M:%SZ"
    ),
    taxon = taxon |>
      str_to_lower() |>
      str_replace_all("\\.", "") |>
      str_squish()
  ) |>
  left_join(station_coords, by = "station_id") |>
  mutate(
    internal_id = paste(
      station_id,
      format(date, "%Y%m%d"),
      depth,
      sep = "-"
    ),
    source = "martin-1999",
    depth_label = "exact",
    s_count_method = "utermohl",
    s_enumeration_method = "light microscopy",
    s_sampling_method = if_else(depth == 0, "bucket", "niskin bottle"),
    s_concentration_method = "sedimentation",
    s_microscopy_method = "light microscope",
    s_light_microscope = "inverted",
    s_preservation = "formalin:acetic acid (1:1)",
    s_volume_sampled = 0.25,
    s_volume_concentrated = 0.05,
    s_counting_chamber = "zeiss"
  ) |>
  select(
    internal_id,
    source,
    datetime,
    year,
    month,
    day,
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
    s_microscopy_method,
    s_light_microscope,
    s_preservation,
    s_volume_sampled,
    s_volume_concentrated,
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
