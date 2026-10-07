# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: DFO-BioChem-AZMP
# Fernando Rayón Viña
# last update: 20260330T00:00 AST
#
# Canada DFO BioChem - AZMP Phytoplankton 1996-2020
# BioChem Query 1942 / AZMP High-Frequency Stations

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(lubridate)
library(janitor)
library(rerddap)
library(DBI)
library(RSQLite)
library(this.path)

# LOAD PATHS ────────────────────────────────────────────────────────────────

script_dir <- this.dir()

setwd(script_dir)

fn_phyto <- file.path(script_dir, "phyto_AZMP", "BioChem_Query_1942_Phyto.csv")

# IMPORT DATA ───────────────────────────────────────────────────────────────

phyto_raw <- read_csv(fn_phyto, show_col_types = FALSE) |> clean_names()

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

phyto_base <- phyto_raw |>
  rename(
    latitude = start_lat,
    longitude = start_lon
  ) |>
  mutate(
    # --- date & time ---
    date = dmy(start_date),
    year = year(date),
    month = month(date),
    day = day(date),
    time_chr = as.character(start_time),
    hours = as.integer(substr(time_chr, 1, nchar(time_chr) - 2)),
    minutes = as.integer(substr(
      time_chr,
      nchar(time_chr) - 1,
      nchar(time_chr)
    )),
    minutes = if_else(minutes >= 60L, minutes - 40L, minutes),
    time = sprintf("%02d:%02d", hours, minutes),
    datetime = format(
      as.POSIXct(paste(date, time), format = "%Y-%m-%d %H:%M", tz = "UTC"),
      "%Y-%m-%dT%H:%M:%SZ"
    ),

    # --- depth: exact vs. range ---
    depth_label = if_else(start_depth == end_depth, "exact", "range"),
    depth = (start_depth + end_depth) / 2,
    depth_label_min = if_else(depth_label == "range", start_depth, NA_real_),
    depth_label_max = if_else(depth_label == "range", end_depth, NA_real_),

    # --- abundance ---
    abundance = c3 / 1000, # cells/m³ → cells/L
    presence_absence = as.integer(c3 > 0),

    # --- identifiers ---
    source = "dfo-biochem-azmp",
    subsource = area_name,
    taxon = taxonomic_name |>
      str_to_lower() |>
      str_replace_all("\\.", "") |>
      str_squish(),

    # --- methods ---
    s_count_method = "Utermohl",
    s_enumeration_method = "Microscopy",
    s_sampling_method = "niskin bottle",
    s_concentration_method = "sedimentation",
    s_volume_sampled = volume * 1000, # m³ → L
    s_microscopy_method = "light microscope",
    s_light_microscope = "inverted",
    s_preservation = preservation,
    s_counting_chamber = "utermohl"
  )

count_data <- phyto_base |>
  filter(!is.na(abundance)) |>
  mutate(
    internal_id = paste(
      "azmp",
      str_extract(descriptor_name, "^\\S+"),
      tsn,
      sep = "_"
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
