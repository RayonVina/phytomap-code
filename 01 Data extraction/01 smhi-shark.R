# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: SMHI-SHARK
# Fernando Rayón Viña
# last update: 20260406T15:00 AST
#
# SMHI SHARK – Swedish National Marine Monitoring Programme
# Phytoplankton (abundance, cells/L)
# SMHI / Havs- och vattenmyndigheten (HaV)
# https://www.smhi.se/en/services/open-data/national-archive-for-oceanographic-data
# Data portal : https://sharkweb.smhi.se | https://sharkdata.se
# R package   : SHARK4R (CRAN)   https://sharksmhi.github.io/SHARK4R/
# Report      : SMHI Report Oceanography No. 81 (2025)

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(lubridate)
library(janitor)
library(DBI)
library(RSQLite)
library(this.path)
library(SHARK4R)

# LOAD PATHS ────────────────────────────────────────────────────────────────

script_dir <- this.dir()

setwd(script_dir)

phyto_raw <- fetch_shark_type("Phytoplankton", "phyto_raw")

# PHYTOPLANKTON ─────────────────────────────────────────────────────────────

phyto_base <- phyto_raw |>
  filter(
    parameter == "Abundance",
    is.na(water_land_station_type_code) |
      water_land_station_type_code %in% c("MO", "C")
  ) |>
  mutate(
    date = sample_date,
    year = year(date),
    month = month(date),
    day = day(date),
    datetime = if_else(
      !is.na(sample_time),
      format(
        as.POSIXct(paste(date, format(sample_time, "%H:%M")), tz = "UTC"),
        "%Y-%m-%dT%H:%M:%SZ"
      ),
      NA_character_
    ),
    depth = case_when(
      !is.na(sample_min_depth_m) &
        !is.na(sample_max_depth_m) &
        sample_min_depth_m != sample_max_depth_m ~ (sample_min_depth_m +
        sample_max_depth_m) /
        2,
      !is.na(sample_min_depth_m) ~ sample_min_depth_m,
      !is.na(sample_max_depth_m) ~ sample_max_depth_m,
      TRUE ~ NA_real_
    ),
    depth_label = if_else(
      !is.na(sample_min_depth_m) &
        !is.na(sample_max_depth_m) &
        sample_min_depth_m != sample_max_depth_m,
      "range",
      "exact"
    ),
    depth_label_min = if_else(
      depth_label == "range",
      sample_min_depth_m,
      NA_real_
    ),
    depth_label_max = if_else(
      depth_label == "range",
      sample_max_depth_m,
      NA_real_
    ),
    abundance = value,
    source = "smhi-shark",
    subsource = str_remove(dataset_name, "^SHARK_Phytoplankton_"),
    taxon = scientific_name |>
      str_to_lower() |>
      str_replace_all("\\.", "") |>
      str_squish(),
    n_size_info = case_when(
      !is.na(size_min_um) & !is.na(size_max_um) ~ paste0(
        size_min_um,
        "-",
        size_max_um,
        "um"
      ),
      !is.na(size_class) ~ paste0("class_", size_class),
      TRUE ~ NA_character_
    ),
    n_marginalia = if_else(
      !is.na(species_flag_code) & species_flag_code != "",
      species_flag_code,
      NA_character_
    ),
    n_notes = paste0(
      if_else(!is.na(aphia_id), paste0("aphia_id=", aphia_id), ""),
      if_else(
        !is.na(reported_cell_volume_um3) & reported_cell_volume_um3 != "",
        paste0(";cell_vol_um3=", reported_cell_volume_um3),
        ""
      )
    ) |>
      na_if(""),
    s_count_method = "Utermohl",
    s_enumeration_method = "Microscopy",
    s_sampling_method = recode(
      sampler_type_code,
      "HOS" = "hose sampler",
      "HOSE" = "hose sampler",
      "NSK" = "niskin bottle",
      "ROS" = "rosette sampler",
      "WAT" = "water bottle",
      "FBW" = "fixed-depth bottle",
      "RU" = "ruttner sampler",
      "RB" = "ruttner bottle",
      "LW" = "large water sampler",
      "IND" = "individual sample",
      "QSM" = "quantitative sample",
      .default = sampler_type_code
    ),
    s_concentration_method = "sedimentation",
    s_volume_sampled = sampled_volume_l,
    s_microscopy_method = "light microscope",
    s_light_microscope = "inverted",
    s_preservation = recode(
      preservation_method_code,
      "CLU" = "Lugol's solution (acid)",
      "ALU" = "Lugol's solution (alkaline)",
      "LUG" = "Lugol's solution",
      .default = preservation_method_code
    ),
    s_volume_concentrated = sedimentation_volume_ml / 1000,
    s_magnification = as.character(magnification),
    s_counting_chamber = "utermohl"
  )

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data <- phyto_base |>
  filter(
    !is.na(abundance),
    !is.na(depth),
    !is.na(date),
    !is.na(sample_latitude_dd),
    !is.na(sample_longitude_dd),
    is.na(quality_flag) | quality_flag %in% c("", "B", " ", "S")
  ) |>
  mutate(
    internal_id = shark_sample_id_md5,
    datetime = case_when(
      !is.na(date) & !is.na(sample_time) ~
        format(
          as.POSIXct(paste(date, format(sample_time, "%H:%M")), tz = "UTC"),
          "%Y-%m-%dT%H:%M:%SZ"
        ),
      !is.na(date) ~
        format(
          as.POSIXct(paste(date, "00:00"), tz = "UTC"),
          "%Y-%m-%dT%H:%M:%SZ"
        ),
      TRUE ~ NA_character_
    ),
    n_uncertain = if_else(quality_flag == "S", 1L, 0L, missing = 0L)
  ) |>
  rename(latitude = sample_latitude_dd, longitude = sample_longitude_dd) |>
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
    depth_label_min,
    depth_label_max,
    taxon,
    abundance,
    n_size_info,
    s_count_method,
    s_enumeration_method,
    s_sampling_method,
    s_concentration_method,
    s_volume_sampled,
    s_microscopy_method,
    s_light_microscope,
    s_preservation,
    s_volume_concentrated,
    s_magnification,
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
