# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Garrison-et-al-1987
# Fernando Rayón Viña
#
# Garrison et al. 1987 — Algal assemblages in Antarctic pack ice
# and in ice-edge plankton. J. Phycol. 23:564-572.
# - count-data: average abundance (cells·L⁻¹) per species and sample
# type (ice / water_column), as per Table 2 of the paper.
# N ice samples = 33; N water column samples = 74.
# - Coordinates: group-level averages of individual station coords
# (see Fig. 1 in paper). No per-species coordinates available.
# - Date: all samples collected 23 Nov–2 Dec 1983; midpoint 28 Nov
# used as representative date (per Excel Readme sheet).
# - Depth: recorded as 50 m in the Excel (represents upper 50 m of
# water column / ice core depth). Treated as depth_label = "range"
# (0–50 m) for water_column; exact for ice core sections (~20 cm,
# coded here as the recorded 50 m average for consistency).
# - "+" values: < 10³ cells·L⁻¹ per paper notation — no exact
# numeric value available, set to NA and excluded.

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(readxl)
library(here)
library(lubridate)
library(janitor)
library(DBI)
library(RSQLite)
library(this.path)

# LOAD PATHS ────────────────────────────────────────────────────────────────

script_dir <- this.dir()

setwd(script_dir)

hereHelper <- function(dataset) {
  paths <- list.dirs(here::here(), recursive = TRUE, full.names = TRUE)
  paths <- paths[basename(paths) == dataset]
  if (length(paths) != 1L) {
    stop(sprintf("Expected one input directory named %s.", dataset))
  }
  paths
}

directory <- hereHelper("Garrison-et-al-1987")

fn <- here(directory, "original-data-Garrison-1987.xlsx")

# IMPORT DATA ───────────────────────────────────────────────────────────────

d1 <- read_excel(fn, sheet = "station-environment-data") |>
  clean_names() |>
  filter(station_id %in% c("Ice", "Water_column")) |>
  mutate(
    depth_label = "exact"
  )

d2 <- read_excel(fn, sheet = "count-data") |>
  clean_names()

# RESHAPE & CLEAN COUNT DATA ────────────────────────────────────────────────

count_data <- d2 |>
  pivot_longer(
    cols = c(ice, water_column),
    names_to = "station_id",
    values_to = "abundance_raw"
  ) |>
  mutate(
    station_id_join = case_when(
      station_id == "ice" ~ "Ice",
      station_id == "water_column" ~ "Water_column"
    ),
    abundance_raw_chr = as.character(abundance_raw),
    abundance = case_when(
      abundance_raw_chr == "+" ~ NA_real_,
      # Excel formula stored as string (e.g. "=6.4*10^6")
      str_detect(abundance_raw_chr, "^=") ~
        as.numeric(str_extract(abundance_raw_chr, "[0-9.]+")) *
        10^as.numeric(str_extract(abundance_raw_chr, "(?<=\\^)[0-9]+")),
      TRUE ~ suppressWarnings(as.numeric(abundance_raw_chr))
    )
  ) |>
  filter(!is.na(abundance)) |>
  left_join(
    d1 |> rename(station_id_join = station_id),
    by = "station_id_join"
  ) |>
  mutate(
    source = "garrison-et-al-1987",
    subsource = station_id,
    n_notes = if_else(
      subsource == "ice",
      "Depth recorded as 50 m representative value; original samples are ice core sections (~20 cm) melted and processed per Garrison & Buck (1986).",
      NA_character_
    )
  ) |>
  select(
    source,
    subsource,
    year,
    month,
    day,
    latitude,
    longitude,
    depth,
    depth_label,
    t_taxon = taxon,
    abundance,
    n_notes
  )

# SAMPLING METADATA ─────────────────────────────────────────────────────────

count_data <- count_data |>
  mutate(
    s_count_method = "utermohl",
    s_enumeration_method = "microscopy",
    s_sampling_method = if_else(subsource == "ice", "ice core", "bottles"),
    s_concentration_method = "sedimentation",
    s_preservation = "glutaraldehyde [1%] (+ paraformaldehyde [1%])",
    s_volume_sampled = if_else(subsource == "ice", "20 cm", "10 mL"),
    s_volume_concentrated = "200 mL",
    s_volume_counted = "10-100 mL",
    s_microscopy_method = "light",
    s_light_microscope = "inverted",
    s_magnification = NA_character_,
    s_counting_chamber = NA_character_
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
