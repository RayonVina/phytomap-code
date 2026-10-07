# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Australian-Phytoplankton-Database-AODN
# Original author: Sadra Dehghani
# Fernando Rayón Viña
# last update: 2025-07-08T15:15 AST
#
# AODN - Australian phytoplankton database

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(stringr)
library(this.path)
library(DBI)
library(RSQLite)

# LOAD PATHS ────────────────────────────────────────────────────────────────

script_dir <- this.dir()

fn1 <- file.path(script_dir, "organized-data.csv")

# IMPORT DATA ───────────────────────────────────────────────────────────────

d <- read_csv(fn1, show_col_types = FALSE) |>
  rename(
    datetime_utc = sample_time_utc,
    taxon = taxon_name,
    sample_id = record_id,
    abundance_l = cells_l,
    presence = presence_absence
  ) |>
  mutate(abundance_l = coalesce(abundance_l, 0))

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data <- d |>
  select(
    subsource = sample_sn,
    datetime = datetime_utc,
    year,
    month,
    day,
    latitude,
    longitude,
    depth = depth_m,
    internal_id = sample_id,
    s_sampling_method = sample_method,
    s_microscopy_method = identification_method,
    s_preservation = preservative,
    t_taxon = taxon,
    abundance = abundance_l
  ) |>
  filter(!is.na(datetime), !is.na(latitude), !is.na(depth)) |>
  mutate(
    source = "aodn",
    subsource = subsource |>
      str_remove(
        "^aus_phyto_db_ongoing_data\\.fid--5c8f4e3e_"
      ) |>
      str_remove(".{5}$"),
    depth_label = "exact",
    s_preservation = case_when(
      str_to_lower(s_preservation) == "lugols" ~ "Lugol's",
      str_to_lower(s_preservation) == "none" ~ "none",
      TRUE ~ s_preservation
    ),
    s_microscopy_method = case_when(
      str_detect(
        str_to_lower(s_microscopy_method),
        "light microscop"
      ) ~ "light microscope",
      str_detect(
        str_to_lower(s_microscopy_method),
        "sem"
      ) ~ "SEM",
      TRUE ~ NA_character_
    ),
    s_enumeration_method = "Light microscopy",
    s_magnification = "640x",
    t_taxon_clean = t_taxon |>
      str_to_lower() |>
      str_replace_all("\\.", "") |>
      str_squish()
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
