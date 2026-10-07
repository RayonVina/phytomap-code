# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Capriulo-Carpent-1983
# Original author: Amelie Frappier
# Fernando Rayón Viña
# last update: 20260330T15:30 ADT
# Original version: 2022-07-08 Amelie Frappier

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(readxl)
library(lubridate)
library(DBI)
library(RSQLite)
library(this.path)

# LOAD PATHS ────────────────────────────────────────────────────────────────

script_dir <- this.dir()

setwd(script_dir)

fn <- file.path(script_dir, "original-data-Capriulo-Cerpenter-1983.xlsx")

# IMPORT DATA ───────────────────────────────────────────────────────────────

station_raw <- read_excel(fn, sheet = "station-environment-data")

count_raw <- read_excel(fn, sheet = "count-data")

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

station_base <- station_raw |>
  mutate(
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude),
    source = "capriulo-1983"
  )

date_cols <- setdiff(colnames(count_raw), c("station_ID", "taxon"))

totals <- tibble(
  date_depth = date_cols,
  total_cells_l = c(
    4220000, # 7-23-1979_1
    1030000, # 8-9-1979_1
    438000, # 8-22-1979_1
    481000, # 9-5-1979_0
    116000, # 10-4-1979_1
    153000, # 11-7-1979_1
    158000, # 12-12-1979_1
    79200, # 1-8-1980_1
    354000, # 1-22-1980_1
    444000, # 2-11-1980_1
    251000, # 3-20-1980_1
    296000, # 4-15-1980_1
    301000, # 5-9-1980_1
    380000, # 5-9-1980_5
    272000, # 6-25-1980_1
    259000, # 6-25-1980_5
    290000, # 7-24-1980_1
    364000, # 7-24-1980_5
    417000, # 8-7-1980_0
    554000, # 8-7-1980_1
    465000, # 8-7-1980_5
    1150000, # 9-5-1980_1
    1620000, # 9-5-1980_5
    1460000, # 9-24-1980_0
    750000, # 10-7-1980_0
    391000 # 10-17-1980_0
  )
)

# COUNT DATA ────────────────────────────────────────────────────────────────

count_data <- count_raw |>
  pivot_longer(
    cols = all_of(date_cols),
    names_to = "date_depth",
    values_to = "percent"
  ) |>
  mutate(percent = as.numeric(percent)) |>
  left_join(totals, by = "date_depth") |>
  mutate(abundance = percent / 100 * total_cells_l) |>
  filter(!is.na(abundance), abundance > 0) |>
  separate(date_depth, into = c("date_str", "depth"), sep = "_") |>
  mutate(
    date = mdy(date_str),
    year = year(date),
    month = month(date),
    day = day(date),
    depth = as.numeric(depth)
  ) |>
  left_join(station_base, by = "station_ID") |>
  mutate(
    internal_id = paste(
      "capriulo-1983",
      station_ID,
      format(date, "%Y%m%d"),
      depth,
      sep = "-"
    ),
    source = "capriulo-1983",
    taxon = taxon |>
      str_to_lower() |>
      str_replace_all("\\.", "") |>
      str_squish(),
    s_count_method = "sedgwick-rafter",
    s_enumeration_method = "Microscopy",
    s_sampling_method = "niskin bottle, bucket",
    s_concentration_method = "centrifuge",
    s_microscopy_method = "light microscope",
    s_preservation = "lugols",
    s_magnification = "500",
    s_counting_chamber = "sedgwick-rafter"
  ) |>
  select(
    internal_id,
    source,
    year,
    month,
    day,
    latitude,
    longitude,
    depth,
    taxon,
    abundance,
    s_count_method,
    s_enumeration_method,
    s_sampling_method,
    s_concentration_method,
    s_microscopy_method,
    s_preservation,
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
