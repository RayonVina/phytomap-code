# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Martin-2014
# Fernando Rayón Viña
# last update: 20260415T15:00 AST
#
# Martin, J.L. & LeGresley, M.M. (2014)
# Phytoplankton Monitoring in the Western Isles Region of the Bay of Fundy
# during 2003–2006
# Can. Tech. Rep. Fish. Aquat. Sci. 3100: v + 190 p.
# Fisheries and Oceans Canada – Biological Station, St. Andrews, NB
# Data digitised from Appendix 2 (phytoplankton densities, cells·L⁻¹, 2003–2006)
# 5 stations: Brandy Cove, Lime Kiln Bay, Deadmans Harbour,
# Wolves Islands (0/10/25/50 m), Mid Passamaquoddy Bay

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

fn <- file.path(script_dir, "original-data-Martin-LeGresley-2014.xlsx")

# IMPORT DATA ───────────────────────────────────────────────────────────────

station_raw <- read_excel(fn, sheet = "station-environment-data")

count_raw <- read_excel(fn, sheet = "count-data")

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

station_base <- station_raw |>
  mutate(
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude)
  )

# STATION METADATA ──────────────────────────────────────────────────────────

station_cols <- setdiff(colnames(count_raw), c("year", "month", "day", "taxon"))

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data <- count_raw |>
  pivot_longer(
    cols = all_of(station_cols),
    names_to = "station_depth",
    values_to = "abundance"
  ) |>
  separate(
    col = station_depth,
    into = c("station_ID", "depth"),
    sep = "_(?=\\d+$)"
  ) |>
  mutate(
    depth = as.numeric(depth),
    abundance = as.numeric(abundance)
  ) |>
  filter(!is.na(abundance)) |>
  left_join(station_base, by = "station_ID") |>
  mutate(
    date = make_date(year, month, day),
    datetime = format(
      as.POSIXct(paste(date, "00:00"), tz = "UTC"),
      "%Y-%m-%dT%H:%M:%SZ"
    ),
    internal_id = paste(
      str_replace_all(station_ID, "\\s+", "-"),
      format(date, "%Y%m%d"),
      depth,
      sep = "-"
    ),
    source = "martin-2014",
    taxon = taxon |>
      str_to_lower() |>
      str_replace_all("\\.", "") |>
      str_squish(),
    depth_label = "exact",
    s_count_method = "utermohl",
    s_enumeration_method = "light microscopy",
    s_sampling_method = if_else(depth == 0, "bucket", "niskin bottle"),
    s_concentration_method = "sedimentation",
    s_volume_sampled = 0.25,
    s_microscopy_method = "light microscope",
    s_light_microscope = "Nikon inverted microscope",
    s_preservation = "formalin:acetic acid (1:1)",
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
    s_volume_sampled,
    s_microscopy_method,
    s_light_microscope,
    s_preservation,
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
