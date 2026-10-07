# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: MAREDAT
# Fernando Rayón Viña
# last update: 20260521T00:00 ADT
#
# MAREDAT – Global phytoplankton compilation
# Sources:
# Diatoms          : PANGAEA.777384  – Leblanc et al. 2012, ESSD 4:149–165
# Coccolithophores : PANGAEA.785092  – O'Brien et al. 2013, ESSD 5:259–275

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(readxl)
library(lubridate)
library(janitor)
library(stringr)
library(this.path)
library(DBI)
library(RSQLite)

# LOAD PATHS ────────────────────────────────────────────────────────────────

script_dir <- this.dir()

setwd(script_dir)

# ── SECTION 1: DIATOMS (PANGAEA.777384) ────────────────────────────────────

fn3 <- "./10.1594_PANGAEA.777384/diatom database-corr-06-12-2016.xlsx"

# IMPORT DATA ───────────────────────────────────────────────────────────────

first_sheet3 <- read_excel(fn3, sheet = 1, guess_max = 100000) |>
  rename_with(~ gsub("\u00b5", "u", .x, fixed = TRUE)) |>
  clean_names() |>
  mutate(
    s_sampling_method = case_when(
      str_detect(str_to_lower(gear_used), "niskin") ~ "niskin bottle",
      str_detect(
        str_to_lower(gear_used),
        "go.?flo|teflon go"
      ) ~ "go-flo bottle",
      str_detect(str_to_lower(gear_used), "nansen") ~ "nansen bottle",
      str_detect(str_to_lower(gear_used), "van dor[mn]") ~ "van dorn bottle",
      str_detect(str_to_lower(gear_used), "ctd") ~ "ctd rosette",
      str_detect(str_to_lower(gear_used), "gpump") ~ "in-situ pump",
      str_detect(
        str_to_lower(gear_used),
        "scanfish"
      ) ~ "continuous sampler (ScanFish)",
      str_detect(str_to_lower(gear_used), "bucket") ~ "bucket",
      str_detect(str_to_lower(gear_used), "bottle") ~ "bottle",
      gear_used == "?" ~ NA_character_,
      !is.na(gear_used) & gear_used != "" ~ gear_used,
      TRUE ~ NA_character_
    ),
    s_count_method = case_when(
      str_detect(
        str_to_lower(couting_method),
        "utermöhl|utermohl|utermol|uthermöhl|sedimentation|settlement|inverted microscop"
      ) ~ "Utermohl",
      str_detect(
        str_to_lower(couting_method),
        "epifluor"
      ) ~ "epifluorescence microscopy",
      str_detect(str_to_lower(couting_method), "flow cytom") ~ "flow cytometry",
      str_to_lower(couting_method) %in%
        c(
          "microscope (taxon counting)",
          "microscope",
          "light microscopy",
          "acid lugol stain"
        ) ~ "Utermohl",
      couting_method == "?" ~ NA_character_,
      !is.na(couting_method) & couting_method != "" ~ couting_method,
      TRUE ~ NA_character_
    ),
    s_enumeration_method = case_when(
      str_detect(
        str_to_lower(couting_method),
        "utermöhl|utermohl|utermol|uthermöhl|sedimentation|settlement|inverted|microscop|light micros|acid lugol stain"
      ) ~ "light microscopy",
      str_detect(
        str_to_lower(couting_method),
        "epifluor"
      ) ~ "epifluorescence microscopy",
      str_detect(str_to_lower(couting_method), "flow cytom") ~ "flow cytometry",
      TRUE ~ NA_character_
    ),
    s_concentration_method = case_when(
      str_detect(
        str_to_lower(couting_method),
        "utermöhl|utermohl|utermol|uthermöhl|sedimentation|settlement|inverted"
      ) ~ "sedimentation",
      TRUE ~ NA_character_
    ),
    s_microscopy_method = case_when(
      str_detect(
        str_to_lower(couting_method),
        "utermöhl|utermohl|utermol|uthermöhl|sedimentation|settlement|inverted|microscop|light micros|acid lugol stain"
      ) ~ "light microscope",
      str_detect(str_to_lower(couting_method), "epifluor") ~ "epifluorescence",
      TRUE ~ NA_character_
    ),
    s_light_microscope = case_when(
      str_detect(
        str_to_lower(couting_method),
        "utermöhl|utermohl|utermol|uthermöhl|sedimentation|settlement|inverted"
      ) ~ "inverted",
      TRUE ~ NA_character_
    ),
    s_counting_chamber = case_when(
      str_detect(
        str_to_lower(couting_method),
        "utermöhl|utermohl|utermol|uthermöhl|sedimentation|settlement|inverted"
      ) ~ "Utermohl",
      TRUE ~ NA_character_
    ),
    s_preservation = case_when(
      str_detect(
        str_to_lower(preservative),
        "lugol.*(glut|fix)|glut.*lugol"
      ) ~ "Lugol + glutaraldehyde",
      str_detect(str_to_lower(preservative), "lugol") ~ "Lugol",
      str_detect(
        str_to_lower(preservative),
        "formalin|formol|formaldehyde"
      ) ~ "formalin",
      str_detect(
        str_to_lower(preservative),
        "glutaraldehyde"
      ) ~ "glutaraldehyde",
      str_detect(str_to_lower(preservative), "bouin") ~ "Bouin",
      !is.na(preservative) & preservative != "" ~ preservative,
      TRUE ~ NA_character_
    )
  ) |>
  select(
    internal_id = primary_key,
    subsource = origin_database,
    year,
    month,
    day,
    date,
    longitude,
    latitude,
    depth,
    taxon = corrected_name_entry,
    abundance = abundance_cells_l_1,
    s_sampling_method,
    s_count_method,
    s_enumeration_method,
    s_concentration_method,
    s_microscopy_method,
    s_light_microscope,
    s_counting_chamber,
    s_preservation
  )

# ── SECTION 2: COCCOLITHOPHORES (PANGAEA.785092) ───────────────────────────

fn5 <- "./10.1594_PANGAEA.785092/ESSD_Coccolithophores_revised_2013- copied.xlsx"

# IMPORT DATA ───────────────────────────────────────────────────────────────

first_sheet5 <- read_excel(fn5, sheet = 1, guess_max = 10000) |>
  rename_with(~ gsub("\u00b5", "u", .x, fixed = TRUE)) |>
  clean_names() |>
  mutate(
    # ── s_sampling_method ─────────────────────────────────────────────────
    s_sampling_method = case_when(
      str_detect(str_to_lower(sampling_method), "niskin") ~ "niskin bottle",
      str_detect(str_to_lower(sampling_method), "nansen") ~ "nansen bottle",
      str_detect(str_to_lower(sampling_method), "rosette") ~ "ctd rosette",
      str_detect(str_to_lower(sampling_method), "ctd") ~ "ctd rosette",
      str_detect(str_to_lower(sampling_method), "gpump") ~ "in-situ pump",
      str_detect(str_to_lower(sampling_method), "multi sonde") ~ "multi sonde",
      str_detect(str_to_lower(sampling_method), "heron net") ~ "net",
      str_detect(str_to_lower(sampling_method), "bottle") ~ "bottle",
      !is.na(sampling_method) & sampling_method != "" ~ sampling_method,
      TRUE ~ NA_character_
    ),

    # ── s_count_method and related fields derived from instrument_method ────────────────

    s_count_method = case_when(
      str_detect(
        str_to_lower(instrument_method),
        "utermohl|utermöhl|sedimentation|inverted"
      ) ~ "Utermohl",
      str_detect(
        str_to_lower(instrument_method),
        "epifluor"
      ) ~ "epifluorescence microscopy",
      str_detect(
        str_to_lower(instrument_method),
        "flow cytom"
      ) ~ "flow cytometry",
      str_detect(str_to_lower(instrument_method), "sem") ~ "SEM",
      str_detect(
        str_to_lower(instrument_method),
        "phase.contrast"
      ) ~ "phase-contrast microscopy",
      str_to_lower(instrument_method) %in%
        c(
          "light microscopy",
          "light microscopy?",
          "lm",
          "light microscope (x400)"
        ) ~ "Utermohl",
      !is.na(instrument_method) &
        !instrument_method %in% c("GF10", "Z", "NI10") ~ instrument_method,
      TRUE ~ NA_character_
    ),
    s_enumeration_method = case_when(
      str_detect(
        str_to_lower(instrument_method),
        "utermohl|utermöhl|sedimentation|inverted|light micros|lm"
      ) ~ "light microscopy",
      str_detect(
        str_to_lower(instrument_method),
        "epifluor"
      ) ~ "epifluorescence microscopy",
      str_detect(
        str_to_lower(instrument_method),
        "flow cytom"
      ) ~ "flow cytometry",
      str_detect(str_to_lower(instrument_method), "sem") ~ "SEM",
      str_detect(
        str_to_lower(instrument_method),
        "phase.contrast"
      ) ~ "light microscopy",
      TRUE ~ NA_character_
    ),
    s_concentration_method = case_when(
      str_detect(
        str_to_lower(instrument_method),
        "utermohl|utermöhl|sedimentation|inverted"
      ) ~ "sedimentation",
      TRUE ~ NA_character_
    ),
    s_microscopy_method = case_when(
      str_detect(
        str_to_lower(instrument_method),
        "utermohl|utermöhl|sedimentation|inverted|light micros|lm"
      ) ~ "light microscope",
      str_detect(
        str_to_lower(instrument_method),
        "epifluor"
      ) ~ "epifluorescence",
      str_detect(
        str_to_lower(instrument_method),
        "phase.contrast"
      ) ~ "light microscope",
      str_detect(str_to_lower(instrument_method), "sem") ~ "SEM",
      TRUE ~ NA_character_
    ),
    s_light_microscope = case_when(
      str_detect(
        str_to_lower(instrument_method),
        "utermohl|utermöhl|sedimentation|inverted"
      ) ~ "inverted",
      str_detect(str_to_lower(instrument_method), "phase.contrast") ~ "upright",
      TRUE ~ NA_character_
    ),
    s_counting_chamber = case_when(
      str_detect(
        str_to_lower(instrument_method),
        "utermohl|utermöhl|sedimentation|inverted"
      ) ~ "Utermohl",
      TRUE ~ NA_character_
    ),

    # ── s_preservation ────────────────────────────────────────────────────

    s_preservation = case_when(
      str_detect(
        str_to_lower(preservative),
        "lugol.*(glut|fix)|glut.*lugol"
      ) ~ "Lugol + glutaraldehyde",
      str_detect(str_to_lower(preservative), "acid lugol") ~ "Lugol (acid)",
      str_detect(str_to_lower(preservative), "lugol") ~ "Lugol",
      str_detect(
        str_to_lower(preservative),
        "formalin|formol|formaldehyde"
      ) ~ "formalin",
      str_detect(
        str_to_lower(preservative),
        "glutaraldehyde"
      ) ~ "glutaraldehyde",
      preservative %in% c("GF10", "NI10") ~ NA_character_,
      !is.na(preservative) & preservative != "" ~ preservative,
      TRUE ~ NA_character_
    )
  ) |>
  select(
    subsource = database,
    year,
    month,
    day,
    longitude = longitutude,
    latitude,
    depth = depth_m,
    taxon = species_lowest_classification,
    abundance = cells_l,
    s_sampling_method,
    s_count_method,
    s_enumeration_method,
    s_concentration_method,
    s_microscopy_method,
    s_light_microscope,
    s_counting_chamber,
    s_preservation
  )

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data <- bind_rows(first_sheet3, first_sheet5) |>
  distinct() |>
  filter(!is.na(abundance), !is.na(depth)) |>
  mutate(
    year = if_else(is.na(year) & !is.na(date), year(date), year),
    month = if_else(is.na(month) & !is.na(date), month(date), month),
  ) |>
  mutate(
    datetime = sprintf("%04d-%02d-%02dT00:00:00Z", year, month, day),
    source = "maredat",
    depth = as.integer(depth),
    across(c(source, subsource, taxon), ~ replace_na(.x, ""))
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
    taxon,
    abundance,
    s_sampling_method,
    s_count_method,
    s_enumeration_method,
    s_concentration_method,
    s_microscopy_method,
    s_light_microscope,
    s_counting_chamber,
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
