# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: SRN-English-Channel
# Fernando Rayón Viña
# last update: 20260318T15:00 AST
#
# PHYTOBS / SRN - French coastal phytoplankton monitoring
# SRN-English-Channel

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(data.table)
library(lubridate)
library(janitor)
library(DBI)
library(this.path)
library(RSQLite)

# LOAD PATHS ────────────────────────────────────────────────────────────────

script_dir <- this.dir()

setwd(script_dir)

fn_phytobs <- file.path(script_dir, "phytobs-data", "combined.csv")

fn_srn_counts <- file.path(script_dir, "seanoe-data", "118659.csv")

# IMPORT DATA ───────────────────────────────────────────────────────────────

phytobs_raw <- read_delim(fn_phytobs, delim = ";", show_col_types = FALSE)

srn_counts_raw <- read_delim(
  fn_srn_counts,
  delim = ";",
  show_col_types = FALSE,
  locale = locale(encoding = "latin1")
)

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

phytobs_base <- phytobs_raw |>
  select(-`...49`) |>
  filter(!is.na(scientific_name) & scientific_name != "") |>
  mutate(
    date = as.Date(sampling_date),
    year = year(date),
    month = month(date),
    day = day(date),
    depth = max_depth,
    subsource = str_to_lower(biotic_dataset_program),
    source = "srn-english-channel"
  )

srn_counts_base <- srn_counts_raw |>
  rename(
    passage_id = `Passage : Identifiant interne`,
    year = `Passage : Année`,
    month = `Passage : Mois`,
    time = `Passage : Heure`,
    latitude = `Coordonnées passage : Coordonnées minx`,
    longitude = `Coordonnées passage : Coordonnées miny`,
    depth = `Prélèvement : Immersion`,
    device = `Libellé de l'engin de prélévement`,
    taxon = `Résultat : Nom du taxon référent`,
    abundance = `Résultat : Valeur de la mesure`
  ) |>
  mutate(
    date = dmy(`Passage : Date`),
    year = as.integer(year),
    month = as.integer(month),
    day = day(date),
    subsource = "srn",
    source = "srn-english-channel"
  )

# COUNT DATA ────────────────────────────────────────────────────────────────

count_data_phytobs <- phytobs_base |>
  filter(!is.na(depth)) |>
  mutate(
    taxon = scientific_name |>
      str_to_lower() |>
      str_replace_all("\\.", "") |>
      str_squish(),
    internal_id = paste(
      subsource,
      sample_name,
      sep = "-"
    ),
    abundance = `single_taxon_count (cells/l)`,
    s_count_method = "Utermohl",
    s_enumeration_method = "Microscopy",
    s_sampling_method = "niskin bottle",
    s_concentration_method = "sedimentation",
    s_microscopy_method = "light microscope",
    s_light_microscope = "inverted",
    s_preservation = if_else(
      preservation_solution == "not available",
      NA_character_,
      preservation_solution
    ),
    s_counting_chamber = "utermohl"
  ) |>
  filter(!is.na(abundance)) |>
  select(
    internal_id,
    source,
    subsource,
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
    s_light_microscope,
    s_preservation,
    s_counting_chamber
  )

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data_srn <- srn_counts_base |>
  mutate(
    taxon = taxon |>
      str_to_lower() |>
      str_replace_all("\\.", "") |>
      str_squish(),
    internal_id = paste(
      subsource,
      passage_id,
      sep = "-"
    ),
    s_count_method = "Utermohl",
    s_enumeration_method = "Microscopy",
    s_sampling_method = "niskin bottle",
    s_concentration_method = "sedimentation",
    s_microscopy_method = "light microscope",
    s_light_microscope = "inverted",
    s_counting_chamber = "utermohl"
  ) |>
  filter(!is.na(abundance)) |>
  select(
    internal_id,
    source,
    subsource,
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
    s_light_microscope,
    s_counting_chamber
  )

count_data <- bind_rows(count_data_phytobs, count_data_srn)

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
