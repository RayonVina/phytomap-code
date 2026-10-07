# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Stoecker-1984
# Fernando Rayón Viña
# last update: 20260427T00:00 ADT
#
# Stoecker, D.K., Davis, L.H. & Anderson, D.M. (1984)
# Fine scale spatial correlations between planktonic ciliates and dinoflagellates
# Journal of Plankton Research, 6(5): 829-842
# https://doi.org/10.1093/plankt/6.5.829
# Diel vertical distribution studies (Study No. 1: Sep 28-29; No. 2: Oct 5-6,
# 1981) at Station A, Perch Pond, Falmouth, MA; and horizontal distribution
# study (Apr 19, 1982) at 17 locations (~50 m spacing) in Perch Pond.
# Abundances are water-column averages (0-4 m, cells L-1).
# Note: Table I of the paper heads the horizontal study date "April 9, 1982";
# the Methods section reads "April 19, 1982". April 19 is used here.

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

fn <- file.path(script_dir, "original-data-Stoecker-et-al-1984.xlsx")

# IMPORT DATA ───────────────────────────────────────────────────────────────

station_raw_aux <- read_excel(fn, sheet = "station-environment-data") |>
  clean_names() |>
  remove_empty(c("rows", "cols"))

count1_raw <- read_excel(fn, sheet = "count-data1") |>
  clean_names() |>
  remove_empty(c("rows", "cols"))

count2_raw <- read_excel(fn, sheet = "count-data2") |>
  clean_names() |>
  remove_empty(c("rows", "cols"))

# HELPERS ───────────────────────────────────────────────────────────────────

coords <- station_raw_aux |>
  select(station_id, latitude, longitude) |>
  distinct(station_id, .keep_all = TRUE)

iso_datetime <- function(year, month, day, time) {
  h <- as.integer(time) %/% 100
  m <- as.integer(time) %% 100
  base <- make_date(year, month, day)
  h_adj <- if_else(h == 24L, 0L, h)
  dt <- as.POSIXct(
    paste0(
      if_else(h == 24, base + days(1), base),
      " ",
      sprintf("%02d:%02d", h_adj, m)
    ),
    tz = "UTC"
  )
  format(dt, "%Y-%m-%dT%H:%M:%SZ")
}

meta_cols <- c("station_id", "year", "month", "day", "time")

s2_cols <- names(count1_raw)[6:9]

s1_cols <- names(count1_raw)[10:14]

taxon_map <- setNames(
  c(
    "Favella sp.",
    "Prorocentrum redfieldi",
    "Prorocentrum minimum",
    "Dinoflagellata indeterminata",
    "Favella sp.",
    "Cochlodinium helicoides",
    "Prorocentrum minimum",
    "Prorocentrum redfieldi",
    "Dinoflagellata indeterminata"
  ),
  c(s2_cols, s1_cols)
)

diel_s2 <- count1_raw |>
  filter(month == 10) |>
  select(all_of(c(meta_cols, s2_cols))) |>
  pivot_longer(all_of(s2_cols), names_to = "taxon_raw", values_to = "abundance")

diel_s1 <- count1_raw |>
  filter(month == 9) |>
  select(all_of(c(meta_cols, s1_cols))) |>
  pivot_longer(all_of(s1_cols), names_to = "taxon_raw", values_to = "abundance")

diel_all <- bind_rows(diel_s2, diel_s1) |>
  filter(!is.na(abundance)) |>
  mutate(taxon = recode(taxon_raw, !!!taxon_map)) |>
  select(-taxon_raw)

horiz <- count2_raw |>
  mutate(
    day = if_else(
      as.integer(day) == 9L &
        as.integer(month) == 4L &
        as.integer(year) == 1982L,
      19L,
      as.integer(day)
    )
  ) |>
  rename(abundance = abundance_l) |>
  select(station_id, year, month, day, time, taxon, abundance)

count_diel <- diel_all |>
  left_join(coords, by = "station_id") |>
  mutate(
    year = as.integer(year),
    month = as.integer(month),
    day = as.integer(day),
    datetime = iso_datetime(year, month, day, time),
    source = "stoecker-1984",
    subsource = "stoecker-1984-diel",
    internal_id = paste(
      station_id,
      sprintf("%04d%02d%02d", year, month, day),
      sprintf("%04d", as.integer(time)),
      str_replace_all(str_to_lower(taxon), "[^a-z0-9]+", "-"),
      sep = "-"
    ),
    depth = 2.0,
    depth_label = "range",
    depth_label_min = 0.0,
    depth_label_max = 4.0,
    n_notes = if_else(
      taxon == "Favella sp.",
      "Favella counted using dissecting microscope after concentration of 1 L subsample through 30 um mesh",
      NA_character_
    ),
    s_count_method = "utermohl",
    s_enumeration_method = "light microscopy",
    s_sampling_method = "bottle",
    s_concentration_method = if_else(
      taxon == "Favella sp.",
      "filtration (30 um mesh)",
      "sedimentation"
    ),
    s_volume_sampled = 2000.0,
    s_volume_counted = if_else(taxon == "Favella sp.", NA_real_, 25.0),
    s_light_microscope = if_else(
      taxon == "Favella sp.",
      "dissecting",
      "inverted"
    ),
    s_preservation = "formalin"
  )

count_horiz <- horiz |>
  left_join(coords, by = "station_id") |>
  mutate(
    year = as.integer(year),
    month = as.integer(month),
    day = as.integer(day),
    datetime = iso_datetime(year, month, day, time),
    source = "stoecker-1984",
    subsource = "stoecker-1984-horizontal",
    internal_id = paste(
      station_id,
      sprintf("%04d%02d%02d", year, month, day),
      sprintf("%04d", as.integer(time)),
      str_replace_all(str_to_lower(taxon), "[^a-z0-9]+", "-"),
      sep = "-"
    ),
    depth = 2.0,
    depth_label = "range",
    depth_label_min = 0.0,
    depth_label_max = 4.0,
    n_notes = "Ciliates counted over entire 100 mL top cylinder; dinoflagellates counted in crossed diameter transects equivalent to 2.3 mL. Station coordinates are a centroid approximation of 17 sampling locations (~50 m spacing); individual positions unavailable.",
    s_count_method = "utermohl",
    s_enumeration_method = "light microscopy",
    s_sampling_method = "integrated water column (PVC pipe)",
    s_concentration_method = "sedimentation",
    s_light_microscope = "inverted",
    s_preservation = "1% buffered formaldehyde"
  )

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data <- bind_rows(count_diel, count_horiz) |>
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
    n_notes,
    s_count_method,
    s_enumeration_method,
    s_sampling_method,
    s_concentration_method,
    s_volume_sampled,
    s_volume_counted,
    s_light_microscope,
    s_preservation
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
