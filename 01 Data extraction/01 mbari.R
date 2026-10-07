# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: MBARI
# Fernando Rayón Viña
# last update: 20260429T00:00 ADT
#
# MBARI Phytoplankton Time Series (MBTS)
# Monterey Bay Aquarium Research Institute (MBARI)
# Monterey Bay and central California coast (USA), 1989–2014
# Source: CeNCOOS / CALOOS MBON GeoServer
# https://data.caloos.org/#search?type_group=all&query=mbari
# WFS: https://data.axds.co/gs/mbon/ows?service=WFS&version=1.0.0&request=GetFeature&outputFormat=csv&typeName=mbon:mbari_phytoplankton_json
# Methods: Epifluorescence microscopy on glutaraldehyde-preserved seawater (~60 mL)
# Ref: CeNCOOS IOOS RA, Appendix F1.4 – Phytoplankton Counts, MBARI
# NOTE: `vf` (volume factor, dimensionless) stored in s_volume_concentrated;
# it is a per-filter conversion factor, not a volume in litres — see metadata.
# NOTE: dataset includes heterotrophic taxa (e.g. Heterotrophic Flagellates);
# no autotroph-only filter applied.

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(lubridate)
library(janitor)
library(DBI)
library(RSQLite)
library(this.path)

# LOAD PATHS ────────────────────────────────────────────────────────────────

script_dir <- this.dir()

setwd(script_dir)

fn <- file.path(script_dir, "mbari_phytoplankton_json.csv")

# IMPORT DATA ───────────────────────────────────────────────────────────────

raw <- read_csv(fn, na = c("NaN", "NA", ""), show_col_types = FALSE) |>
  clean_names()

# TRANSFORM: COUNT DATA ─────────────────────────────────────────────────────

count_data <- raw |>
  mutate(
    dt = coalesce(
      suppressWarnings(ymd_hms(datetime_gmt, tz = "UTC")),
      suppressWarnings(dmy_hms(date_str, tz = "UTC"))
    ),
    datetime = paste0(format(dt, "%Y-%m-%dT%H:%M:%S"), "Z"),
    year = year(dt),
    month = month(dt),
    day = day(dt),
    internal_id = paste(
      str_remove(fid, "^mbari_phytoplankton_json\\.fid-"),
      cruise,
      station_name,
      groups,
      sep = "_"
    ),
    taxon = param_id |>
      str_to_lower() |>
      str_replace_all("\\.", "") |>
      str_squish(),
    abundance = replace_na(as.numeric(perliter), 0),
    source = "mbari",
    subsource = str_to_lower(project),
    depth_label = "exact",
    s_count_method = "epifluorescence microscopy",
    s_enumeration_method = "epifluorescence microscopy",
    s_sampling_method = "bottle",
    s_preservation = "glutaraldehyde",
    s_microscopy_method = "epifluorescence",
    s_volume_concentrated = as.numeric(vf)
  ) |>
  filter(
    !is.na(target_depth_m),
    !is.na(dec_lat),
    !is.na(dec_long),
    !is.na(datetime_gmt)
  ) |>
  rename(
    latitude = dec_lat,
    longitude = dec_long,
    depth = target_depth_m,
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
    depth,
    depth_label,
    taxon,
    abundance,
    s_count_method,
    s_enumeration_method,
    s_sampling_method,
    s_preservation,
    s_microscopy_method,
    s_volume_concentrated
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
