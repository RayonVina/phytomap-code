# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Australian-reference-stations
# Original author: Sadra Dehghani; Julianne Jager
# Fernando Rayón Viña
# last update: 20250616T16:10 AST
# Original version: 2022-12-12 Julianne Jager
#
# Australian reference stations

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(readxl)
library(this.path)
library(DBI)
library(RSQLite)
library(lubridate)

# LOAD PATHS ────────────────────────────────────────────────────────────────

script_dir <- this.dir()

fn1 <- file.path(
  script_dir,
  "IMOS_-_Phytoplankton_Abundance_and_Biovolume_(reference_stations)-abundance_-_raw_data.csv"
)

# IMPORT DATA ───────────────────────────────────────────────────────────────

d <- read_csv(fn1, show_col_types = FALSE)

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data <- d |>
  pivot_longer(
    cols = c(`Acantharia`:`unid square cells`),
    names_to = "taxon",
    values_to = "abundance_l"
  ) |>
  select(
    datetime_utc = SampleTime_UTC,
    latitude = Latitude,
    longitude = Longitude,
    depth_m = SampleDepth_m,
    station_id = StationCode,
    taxon,
    abundance_l,
    count_method = Method
  ) |>
  filter(!is.na(abundance_l)) |>
  mutate(
    year = year(datetime_utc),
    month = month(datetime_utc),
    day = day(datetime_utc),
    source = "aus-ref",
    count_method = recode(count_method, "LM" = "light microscopy"),
    taxon_original = taxon,
    taxon_clean = taxon |>
      str_to_lower() |>
      str_replace_all("\\.", "") |>
      str_squish()
  ) |>
  separate(
    taxon_clean,
    into = c("tax_genus", "tax_species", "tax_epithet"),
    sep = " ",
    extra = "merge",
    fill = "right"
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
