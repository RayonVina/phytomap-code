# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: LTER-NAS-v4
# Original author: Sadra Dehghani
# Fernando Rayón Viña
# last update: 20260316T15:00 ADT
#
# LTER-NAS v4 - North Atlantic Study (NAS) Phytoplankton 1965-2015

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(janitor)
library(this.path)
library(RSQLite)
library(DBI)
library(lubridate)

# LOAD PATHS ────────────────────────────────────────────────────────────────

script_dir <- this.dir()

fn_data <- file.path(script_dir, "database_LTER_NAS_1965_2015_v4.csv")

# RAW DATA ──────────────────────────────────────────────────────────────────

data_raw <- read_csv(fn_data, show_col_types = FALSE) |>
  clean_names() |>
  select(-ends_with("_sensor_uri")) |> # drop instrument metadata columns
  rename(
    longitude = long,
    latitude = lat,
    depth_m = depth,
    bottom_depth_m = bot_depth_m,
    date = yyyy_mm_dd,
    time = hh_mm_ss,
    temperature_c = temp,
    salinity = sal,
    density = dens,
    oxygen_ml_l = oxyg_ml_l,
    oxygen_percent = ox_percent,
    ammonia_um = nh3_micro_mol,
    nitrite_um = no2_micro_mol, # NO2 = nitrite (not nitrate)
    nitrate_um = no3_micro_mol, # NO3 = nitrate
    din_um = din_micro_mol,
    phosphate_um = po4_micro_mol,
    silicate_um = si_micro_mol,
    alkalinity = alky,
    total_nitrogen_um = ntot,
    total_phytoplankton_ml = phyto_tot_cell_ml,
    zooplankton_ind_m3 = zoo_ind_m_3,
    diatom_ml = diato_cell_ml_urn_lsid_marinespecies_org_taxname_148899,
    dinoflagellate_ml = dino_cell_ml_urn_lsid_marinespecies_org_taxname_19542,
    flagellate_ml = flag_cell_ml,
    cocco_ml = cocco_cell_ml
  ) |>
  mutate(
    date = ymd(date),
    year = year(date),
    month = month(date),
    day = day(date)
  )

# COUNT DATA ────────────────────────────────────────────────────────────────

count_data <- data_raw |>
  select(
    longitude,
    latitude,
    station,
    cruise,
    ship,
    year,
    month,
    day,
    depth_m,
    diatom_ml,
    dinoflagellate_ml,
    cocco_ml
  ) |>
  group_by(station, cruise, ship) |>
  mutate(internal_id = paste(station, cruise, ship, sep = "-")) |>
  ungroup() |>
  mutate(across(
    c(diatom_ml, dinoflagellate_ml, cocco_ml),
    ~ suppressWarnings(as.numeric(.x)) * 1000
  )) |>
  pivot_longer(
    cols = c(diatom_ml, dinoflagellate_ml, cocco_ml),
    names_to = "taxon",
    values_to = "abundance_l"
  ) |>
  mutate(
    taxon = recode(
      taxon,
      "diatom_ml" = "Bacillariophyceae",
      "dinoflagellate_ml" = "Dinophyceae",
      "cocco_ml" = "Coccolithophyceae"
    )
  ) |>
  filter(!is.na(abundance_l)) |>
  mutate(source = "lter-nas") |>
  select(
    internal_id,
    latitude,
    longitude,
    depth_m,
    year,
    month,
    day,
    taxon,
    abundance_l,
    source
  ) |>
  distinct()

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
