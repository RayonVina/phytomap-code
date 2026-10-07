# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Thompson_Bonham-2011
# Fernando Rayón Viña
# last update: 20260427T00:00 ADT
#
# Thompson, P. A. & Bonham, P. (2011)
# New insights into the Kimberley phytoplankton and their ecology.
# Journal of the Royal Society of Western Australia, 94: 161-169.
# Data from Table 1: mean phytoplankton abundances (cells L-1) from
# transects A and C (n = 12), Kimberley region, April-May 2010.
# Coordinates are approximate (midpoint of transects A-C).
# Depth set to 25 m by convention (water-column average: surface to ~50 m).
# s_volume_sampled set to 3 L (midpoint of reported 1-5 L range).

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

fn <- file.path(script_dir, "original-data-Thompson-Bonham-2011.xlsx")

# IMPORT DATA ───────────────────────────────────────────────────────────────

station_raw <- read_excel(fn, sheet = "station-environment-data") |>
  clean_names() |>
  remove_empty(c("rows", "cols"))

count_raw <- read_excel(fn, sheet = "count-data") |>
  clean_names() |>
  remove_empty(c("rows", "cols"))

# TRANSFORM ─────────────────────────────────────────────────────────────────

count_fixed <- count_raw |>
  mutate(
    taxon = case_when(
      taxon == "Chaetoceros spp." &
        abundance_cells_l == 1639 ~ "Chaetoceros spp. < 10 um",
      taxon == "Chaetoceros spp." &
        abundance_cells_l == 1461 ~ "Chaetoceros spp. > 10 um",
      TRUE ~ taxon
    )
  )

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data <- count_fixed |>
  left_join(station_raw, by = "station_id") |>
  mutate(
    year = as.integer(year),
    month = as.integer(month),
    day = as.integer(day),
    date = make_date(year, month, day),
    datetime = format(
      as.POSIXct(paste(date, "00:00"), tz = "UTC"),
      "%Y-%m-%dT%H:%M:%SZ"
    ),
    taxon = taxon |>
      str_to_lower() |>
      str_replace_all("\\.", "") |>
      str_squish(),
    internal_id = paste(station_id, taxon, sep = "-"),
    source = "thompson-bonham-2011",
    abundance = abundance_cells_l,
    depth = 25.0,
    depth_label = "wc",
    depth_label_min = 0.0,
    depth_label_max = 50.0,
    s_count_method = "sedgwick-rafter",
    s_enumeration_method = "light microscopy",
    s_sampling_method = "niskin bottle",
    s_concentration_method = "sedimentation",
    s_volume_sampled = 3.0,
    s_microscopy_method = "light microscope",
    s_light_microscope = "inverted",
    s_preservation = "acid lugols solution",
    s_volume_counted = 1.0,
    s_counting_chamber = "sedgwick-rafter"
  ) |>
  filter(!is.na(abundance)) |>
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
    depth_label_min,
    depth_label_max,
    abundance,
    taxon,
    s_count_method,
    s_enumeration_method,
    s_sampling_method,
    s_concentration_method,
    s_volume_sampled,
    s_microscopy_method,
    s_light_microscope,
    s_preservation,
    s_volume_counted,
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
