# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Narragansett
# Fernando Rayón Viña
# last update: 20260320T15:30 AST
#
# NARRAGANSETT BAY LONG-TERM PLANKTON TIME SERIES
# Historical: Smayda 1959-1997 (URI GSO / NABATS)
# Modern: URI GSO 1999-present

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(readxl)
library(lubridate)
library(DBI)
library(this.path)
library(RSQLite)

# DATA SOURCES ──────────────────────────────────────────────────────────────

script_dir <- this.dir()

# LOAD PATHS ────────────────────────────────────────────────────────────────

setwd(script_dir)

fn_hist <- file.path(
  script_dir,
  "smayda-phytoplankton-narragansett-bay-station-2-surface-1959-to-1997.xlsx"
)

fn_mod <- file.path(script_dir, "Cell-Count-Master-Data.xlsx")

LAT <- 41.570

LON <- -71.390

# IMPORT ────────────────────────────────────────────────────────────────────

hist_raw <- read_excel(fn_hist, sheet = "Stn II surface phyto", skip = 9) |>
  slice(-c(1, 2))

# IMPORT DATA ───────────────────────────────────────────────────────────────

mod_raw <- read_excel(fn_mod, sheet = "Count data") |>
  slice(-1)

# TRANSFORM HISTORICAL ──────────────────────────────────────────────────────

count_hist <- hist_raw |>
  rename(date = `Sample Date`) |>
  filter(!is.na(date)) |>
  mutate(
    date = as.Date(date),
    year = year(date),
    month = month(date),
    day = day(date)
  ) |>
  select(-`Time Series Week`) |>
  pivot_longer(
    cols = -c(date, year, month, day),
    names_to = "taxon",
    values_to = "abundance"
  ) |>
  filter(!is.na(abundance)) |>
  mutate(
    taxon = taxon |>
      str_remove("\\.\\.\\.\\d+$") |>
      str_to_lower() |>
      str_replace_all("\\.", "") |>
      str_squish(),
    abundance = abundance * 1000, # cells/mL → cells/L
    internal_id = paste(
      "sta2",
      "historical",
      format(date, "%Y%m%d"),
      sep = "_"
    ),
    source = "narragansett",
    subsource = "historical",
    latitude = LAT,
    longitude = LON,
    depth = 0.5,
    depth_label = "range",
    depth_label_min = 0.0,
    depth_label_max = 1.0,
    s_count_method = "Sedgwick-Rafter",
    s_enumeration_method = "Microscopy",
    s_sampling_method = "water bottle",
    s_microscopy_method = "light microscope",
    s_light_microscope = "compound",
    s_magnification = "250x-500x",
    s_counting_chamber = "Sedgwick-Rafter cell",
    notes = "surface sample; depth nominal 0.5 m (±0.5 m); whole-water unconcentrated counts"
  )

# TRANSFORM MODERN ──────────────────────────────────────────────────────────

count_mod <- mod_raw |>
  rename(date = DATE, count_type = `COUNT TYPE`, location = LOCATION) |>
  filter(!is.na(date), !is.na(count_type)) |>
  mutate(
    date = as.Date(date),
    year = year(date),
    month = month(date),
    day = day(date),
    station = case_when(
      str_detect(location, regex("sta\\.?\\s*2", ignore_case = TRUE)) ~ "sta2",
      str_detect(location, regex("gso", ignore_case = TRUE)) ~ "gso-dock",
      str_detect(
        location,
        regex("jamestown", ignore_case = TRUE)
      ) ~ "jamestown-bridge",
      location == "0" ~ "unknown",
      TRUE ~ "unknown"
    ) |>
      str_remove_all("[^a-z0-9]") |>
      str_to_lower(),
    subsource = case_when(
      str_detect(
        count_type,
        regex("mixed", ignore_case = TRUE)
      ) ~ "modern-mixed",
      str_detect(
        count_type,
        regex("surface", ignore_case = TRUE)
      ) ~ "modern-surface",
      str_detect(
        count_type,
        regex("depth", ignore_case = TRUE)
      ) ~ "modern-depth",
      TRUE ~ "modern-other"
    ),
    depth = case_when(
      subsource %in% c("modern-surface", "modern-mixed") ~ 0.5,
      subsource == "modern-depth" ~ 8.0,
      TRUE ~ NA_real_
    ),
    depth_label = case_when(
      subsource %in% c("modern-surface", "modern-depth") ~ "range",
      subsource == "modern-mixed" ~ "wc",
      TRUE ~ NA_character_
    ),
    depth_label_min = case_when(
      subsource == "modern-surface" ~ 0.0,
      subsource == "modern-depth" ~ 7.5,
      subsource == "modern-mixed" ~ 0.0,
      TRUE ~ NA_real_
    ),
    depth_label_max = case_when(
      subsource == "modern-surface" ~ 1.0,
      subsource == "modern-depth" ~ 8.5,
      subsource == "modern-mixed" ~ 8.0,
      TRUE ~ NA_real_
    )
  ) |>
  select(-`Total abundance`, -location, -count_type) |>
  pivot_longer(
    cols = -c(
      date,
      year,
      month,
      day,
      subsource,
      station,
      depth,
      depth_label,
      depth_label_min,
      depth_label_max
    ),
    names_to = "taxon",
    values_to = "abundance"
  ) |>
  mutate(
    abundance = suppressWarnings(as.numeric(abundance)),
    taxon = taxon |>
      str_remove("\\.\\.\\.\\d+$") |>
      str_to_lower() |>
      str_replace_all("\\.", "") |>
      str_squish()
  ) |>
  filter(!is.na(abundance)) |>
  group_by(
    date,
    year,
    month,
    day,
    subsource,
    station,
    depth,
    depth_label,
    depth_label_min,
    depth_label_max,
    taxon
  ) |>
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") |>
  mutate(
    internal_id = paste(
      "narragansett",
      station,
      subsource,
      format(date, "%Y%m%d"),
      sep = "_"
    ),
    source = "narragansett",
    latitude = LAT,
    longitude = LON,
    s_count_method = "Sedgwick-Rafter",
    s_enumeration_method = "Microscopy",
    s_sampling_method = "water bottle",
    s_concentration_method = NA_character_,
    s_microscopy_method = "light microscope",
    s_light_microscope = "compound",
    s_preservation = "1% Lugol's iodine",
    s_magnification = NA_character_,
    s_counting_chamber = "Sedgwick-Rafter cell",
    notes = case_when(
      subsource == "modern-surface" ~
        "surface sample counted separately; depth nominal 0.5 m (±0.5 m)",
      subsource == "modern-depth" ~
        "bottom sample counted separately; depth nominal 8 m (±0.5 m), based on NBNERR Tech Report 2009 and Durbin & Durbin 1998",
      subsource == "modern-mixed" ~
        "equal volumes surface and bottom combined before counting; depth nominal 4 m (±0.5 m, midpoint of ~8 m water column)",
      TRUE ~ "count type unspecified"
    )
  )

# COMBINE & SELECT ──────────────────────────────────────────────────────────

cols <- c(
  "internal_id",
  "source",
  "subsource",
  "year",
  "month",
  "day",
  "latitude",
  "longitude",
  "depth",
  "depth_label",
  "depth_label_min",
  "depth_label_max",
  "taxon",
  "abundance",
  "s_count_method",
  "s_enumeration_method",
  "s_sampling_method",
  "s_microscopy_method",
  "s_light_microscope",
  "s_magnification",
  "s_counting_chamber",
  "notes"
)

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data <- bind_rows(
  count_hist |> select(all_of(cols)),
  count_mod |> select(all_of(cols))
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
