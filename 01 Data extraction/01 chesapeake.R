# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Chesapeake
# Original author: Julianne Jager; Sadra Dehghani
# Fernando Rayón Viña
# last update: 20250710T14:10 AST
# Original version: 2022-12-12 Julianne Jager
# Original version: 2023-05-31 Sadra Dehghani
#
# CHESAPEAKE BAY
# Modified to add simplified depth information and correct LayerSampled recoding

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(readr)
library(dplyr)
library(lubridate)
library(hms)
library(stringr)
library(janitor)
library(this.path)
library(DBI)
library(RSQLite)

# LOAD PATHS ────────────────────────────────────────────────────────────────

script_dir <- this.dir()

setwd(script_dir)

# SIMPLIFIED FUNCTION TO ADD DEPTH INFORMATION ──────────────────────────────

add_simplified_depth_info <- function(data) {
  data |>
    mutate(
      depth_label = case_when(
        layer_sampled == "integrated" ~ "wc",
        TRUE ~ "range"
      ),
      depth_label_min = case_when(
        layer_sampled == "above_pycnocline" ~ 0.5,
        layer_sampled == "below_pycnocline" ~ 4.0,
        layer_sampled == "surface" ~ 0.5,
        layer_sampled == "bottom" ~ 8.0,
        layer_sampled == "integrated" ~ 0.5,
        layer_sampled == "middle" ~ 3.0,
        TRUE ~ NA_real_
      ),
      depth_label_max = case_when(
        layer_sampled == "above_pycnocline" ~ 3.5,
        layer_sampled == "below_pycnocline" ~ 9.5,
        layer_sampled == "surface" ~ 1.0,
        layer_sampled == "bottom" ~ 10.0,
        layer_sampled == "integrated" ~ 10.0,
        layer_sampled == "middle" ~ 6.0,
        TRUE ~ NA_real_
      ),
      depth = (depth_label_min + depth_label_max) / 2
    )
}

# SPECIES DATA: IMPORT AND PROCESS ──────────────────────────────────────────

b <- file.path(script_dir, "PHYTP-station-1986-2025.csv")

cols_to_use <- c(
  "Station",
  "FieldActivityId",
  "BiologicalFieldActivityId",
  "Project",
  "DataProvider",
  "DataCollector",
  "Latitude",
  "Longitude",
  "SampleDate",
  "SampleTime",
  "SampleReplicate",
  "LayerSampled",
  "SampleType",
  "SamplingEquipment",
  "SampleMeasurement",
  "SampleVolume",
  "SampleVolumeUnits",
  "BiologicalAnalyticalMethodCode",
  "SalinityRegime",
  "TaxonomicSerialNumber_BAY",
  "ScientificName",
  "TaxonSynonymStatus",
  "TaxonAttribute",
  "ReportingParameter",
  "ReportingValue",
  "ReportingUnits"
)

col_types <- cols(
  Station = col_character(),
  FieldActivityId = col_double(),
  BiologicalFieldActivityId = col_double(),
  Project = col_character(),
  DataProvider = col_character(),
  DataCollector = col_character(),
  Latitude = col_double(),
  Longitude = col_double(),
  SampleDate = col_date(format = "%m/%d/%Y"),
  SampleTime = col_time("%H:%M:%S"),
  SampleReplicate = col_character(),
  LayerSampled = col_character(),
  SampleType = col_character(),
  SamplingEquipment = col_character(),
  SampleMeasurement = col_character(),
  SampleVolume = col_double(),
  SampleVolumeUnits = col_character(),
  BiologicalAnalyticalMethodCode = col_character(),
  SalinityRegime = col_character(),
  TaxonomicSerialNumber_BAY = col_double(),
  ScientificName = col_character(),
  TaxonSynonymStatus = col_character(),
  TaxonAttribute = col_character(),
  ReportingParameter = col_character(),
  ReportingValue = col_double(),
  ReportingUnits = col_character()
)

# IMPORT DATA ───────────────────────────────────────────────────────────────

fn_b <- read_csv(b, col_types = col_types, col_select = all_of(cols_to_use))

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data <- fn_b |>
  mutate(
    SampleDate_fmt = format(SampleDate, "%Y-%m-%d"),
    datetime_utc = as.POSIXct(
      paste(SampleDate_fmt, "T", as.character(SampleTime), sep = ""),
      format = "%Y-%m-%dT%H:%M:%S",
      tz = "UTC"
    ),
    day = day(datetime_utc),
    month = month(datetime_utc),
    year = year(datetime_utc),
    source = "chesapeake",
    count_method = "utermohl",
    taxon_clean = ScientificName |>
      str_replace_all("\\.", "") |>
      str_replace_all("_", " ") |>
      str_squish(),
    layer_sampled = recode(
      LayerSampled,
      "AP" = "above_pycnocline",
      "BP" = "below_pycnocline",
      "WC" = "integrated",
      "S" = "surface",
      "SUR" = "surface",
      "BOT" = "bottom",
      "INT" = "integrated",
      "MID" = "middle",
      .default = NA_character_
    )
  ) |>
  separate(
    taxon_clean,
    into = c("tax_genus", "tax_species", "tax_epithet"),
    sep = " ",
    extra = "merge",
    fill = "right",
    remove = FALSE
  ) |>
  rename(
    station = Station,
    internal_source = DataCollector,
    internal_sample_replicate = SampleReplicate,
    taxonomic_serial_number = TaxonomicSerialNumber_BAY,
    abundance_l = ReportingValue
  ) |>
  select(
    station,
    Latitude,
    Longitude,
    layer_sampled,
    datetime_utc,
    day,
    month,
    year,
    t_taxon = ScientificName,
    taxon_clean,
    tax_genus,
    tax_species,
    tax_epithet,
    taxonomic_serial_number,
    abundance_l,
    count_method,
    source,
    internal_source,
    internal_sample_replicate
  ) |>
  janitor::clean_names() |>
  add_simplified_depth_info()

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
