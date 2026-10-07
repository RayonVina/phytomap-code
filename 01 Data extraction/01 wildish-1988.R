# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Wildish-1988
# Fernando Rayón Viña
# last update: 20260416T11:51 ADT
#
# Wildish, D.J., Martin, J.L., Wilson, A.J., & DeCoste, A.M. (1988)
# Environmental Monitoring of the Bay of Fundy Salmonid Mariculture Industry
# During 1986 and 1987.
# Can. Tech. Rep. Fish. Aquat. Sci. 1648: iii + 44 p.
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
  "original-data-Wildish-Martin-Wilson-Decoste-1988.xlsx"
)

# IMPORT DATA ───────────────────────────────────────────────────────────────

station_raw <- read_excel(fn, sheet = "station-environment-data-1") |>
  clean_names() |>
  remove_empty(c("rows", "cols"))

count_raw <- read_excel(fn, sheet = "count-data") |>
  remove_empty(c("rows", "cols"))

# TRANSFORM: COUNT DATA ─────────────────────────────────────────────────────

station_cols <- setdiff(colnames(count_raw), c("taxon", "date", "depth_m"))

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data_raw <- count_raw |>
  mutate(
    depth_label = if_else(depth_m == "bottom", "bottom", NA_character_),
    depth_m = suppressWarnings(as.numeric(depth_m))
  ) |>
  pivot_longer(
    cols = all_of(station_cols),
    names_to = "station_id",
    values_to = "abundance"
  ) |>
  replace_na(list(abundance = Inf)) |> # NA = station not sampled
  pivot_wider(
    names_from = station_id,
    values_from = abundance,
    values_fill = 0 # reconstruct true-absence zeros
  ) |>
  pivot_longer(
    cols = all_of(station_cols),
    names_to = "station_id",
    values_to = "abundance"
  ) |>
  filter(is.finite(abundance))

count_data <- count_data_raw |>
  filter(!is.na(depth_m)) |>
  mutate(
    station_id = as.integer(str_remove(station_id, "station_")),
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
  left_join(
    station_raw |> mutate(station_id = as.integer(station_id)),
    by = "station_id"
  ) |>
  mutate(
    internal_id = paste(
      station_id,
      format(date, "%Y%m%d"),
      depth_m,
      sep = "-"
    ),
    source = "wildish-1988",
    depth_label = "exact",
    s_count_method = "utermohl",
    s_enumeration_method = "light microscopy",
    s_sampling_method = if_else(depth_m == 0, "bucket", "nansen bottle"),
    s_concentration_method = "sedimentation",
    s_microscopy_method = "light microscope",
    s_light_microscope = "Zeiss inverted microscope",
    s_preservation = "2.5% formalin:acetic acid mixture (1:1)",
    s_volume_sampled = 0.2,
    s_volume_concentrated = 0.05,
    s_magnification = "10x-40x objective, 12.5x eyepiece"
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
    depth = depth_m,
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
    s_magnification
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
