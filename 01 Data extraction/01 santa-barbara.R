# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: santa-barbara
# Fernando Rayón Viña
# last update: 20260317T17:00 AST
#
# SBC LTER - Santa Barbara Channel Cross-Shelf Study 2008-2009
# knb-lter-sbc.45.8

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(data.table)
library(lubridate)
library(DBI)
library(this.path)
library(RSQLite)

# LOAD PATHS ────────────────────────────────────────────────────────────────

script_dir <- this.dir()

setwd(script_dir)

fn_ctd <- file.path(
  script_dir,
  "SBC_LTER_cross_shelf_study_2008_CTD_downcasts.txt"
)

fn_phyto <- file.path(
  script_dir,
  "SBC_LTER_cross_shelf_study_2008_5m_phyto_counts.csv"
)

fn_bottle <- file.path(
  script_dir,
  "SBC_LTER_cross_shelf_study_2008_profile_bottle_data.csv"
)

# IMPORT DATA ───────────────────────────────────────────────────────────────

ctd_raw <- read_delim(fn_ctd, show_col_types = FALSE)

phyto_raw <- read_csv(fn_phyto, show_col_types = FALSE)

bottle_raw <- fread(fn_bottle, encoding = "Latin-1")

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

ctd_base <- ctd_raw |>
  rename(
    cruise = `Cruise`,
    station = ` Station`,
    type = ` Type`,
    date = ` mon/day/yr`,
    time = ` hh:mm`,
    longitude = ` longitude [degrees_east]`,
    latitude = `latitude [degrees_north]`,
    bottom_depth_m = ` Bot. Depth [m]`,
    depth_m = ` Depth SW [m]`,
    pressure_db = ` Pressure [db]`,
    temp_C = ` potemp [C]`,
    salinity_psu = ` Salinity [PSU]`,
    sigma_theta = ` Sigma-theta`,
    fluorescence = ` Fluor`,
    transmission = ` Trans [%]`,
    beam_atten_m = ` Beam C [1/m]`
  ) |>
  mutate(
    date = mdy(date),
    year = year(date),
    month = month(date),
    day = day(date),
    source = "santa-barbara",
    across(where(is.numeric), ~ na_if(., 99999))
  )

bottle_base <- bottle_raw |>
  rename(
    seq_num = `Seq #`,
    cruise = Cruise,
    station = Station,
    time = `hh:mm`,
    latitude = `latitude [degrees_north]`,
    longitude = `longitude [degrees_west]`,
    depth_target_m = `Depth [target m]`,
    depth_actual_m = `PROF. Depth SW [m]`,
    bottom_depth_m = `Bot. Depth [m]`,
    prof_temp_C = `PROF. T00 [C]`,
    prof_salinity_psu = `PROF. Salinity [PSU]`,
    prof_sigma_theta = `PROF. Sigma-theta`,
    phosphate_uM = `PO4 [umol/l]`,
    silicate_uM = `Si [umol/l]`,
    nitrite_uM = `NO2 [umol/l]`,
    nitrate_nitrite_uM = `NO2+NO3 [umol/l]`,
    ammonium_uM = `NH4 [umol/l]`,
    din_uM = `DIN [µM N]`,
    tdn_uM = `TDN [uM N]`,
    don_uM = `DON [µM N]`,
    poc_uM = `POC [umol/l]`,
    pon_uM = `PON [umol/l]`,
    chla_ug_l = `Chl a [ug/L]`,
    phaeo_ug_l = `Phaeopigments [ug/l]`,
    chla_phaeo_ug_l = `Chl a + phaeo [ug/l]`,
    primary_prod_mgC_m3_d = `Primary Production [mg C/m3/d]`,
    er_mgC_m3_d = `ER [mgC/m3/d]`,
    per = PER,
    doc_uM = `DOC [uM C]`,
    leu_pmol_l_hr = `3H Leu [pmol/L/hr]`,
    bac_prod_mgC_m3_d = `Bacterial Production [mgC/m3/d]`,
    bac_abund_cells_ml = `Bac Abund [cells/ml]`,
    diatom = `Diatoms [cells/L]`,
    dinoflagellate = `Dinoflagellates [cells/L]`,
    silicoflagellate = `Silicoflagellates [cells/L]`,
    prymnesiophyte = `Prymnesiophytes [cells/L]`
  ) |>
  mutate(
    date = as.Date(Date),
    year = year(date),
    month = month(date),
    day = day(date),
    cruise = str_remove(as.character(cruise), "^20"),
    source = "santa-barbara",
    across(where(is.numeric), ~ na_if(., 99999)),
    across(where(is.numeric), ~ na_if(., 100000)),
    across(where(is.numeric), ~ if_else(. < 0, NA_real_, .))
  )

ctd_coords <- ctd_base |>
  mutate(station_num = str_extract(station, "\\d$")) |>
  group_by(station_num, date) |>
  summarise(
    latitude = first(na.omit(latitude)),
    longitude = first(na.omit(longitude)),
    .groups = "drop"
  )

# COUNT DATA ────────────────────────────────────────────────────────────────

count_data_phyto <- phyto_raw |>
  rename(cruise = Cruise, station = Station) |>
  mutate(
    date = dmy(date),
    year = year(date),
    month = month(date),
    day = day(date),
    station_num = as.character(station)
  ) |>
  pivot_longer(
    cols = `Pseudo-nitzschia`:Phaeocystis,
    names_to = "taxon",
    values_to = "abundance"
  ) |>
  left_join(ctd_coords, by = c("station_num", "date")) |>
  filter(!is.na(abundance)) |>
  mutate(
    internal_id = paste(
      cruise,
      station,
      sep = "-"
    ),
    source = "santa-barbara",
    depth = 5,
    taxon = taxon |>
      str_to_lower() |>
      str_replace_all("\\.", "") |>
      str_squish(),
    s_count_method = "Utermohl",
    s_enumeration_method = "Microscopy",
    s_sampling_method = "niskin bottle",
    s_concentration_method = "sedimentation",
    s_microscopy_method = "light microscope",
    s_light_microscope = "inverted",
    s_preservation = "formaldehyde 2%",
    s_magnification = "350",
    s_counting_chamber = "utermohl"
  )

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data_bottle <- bottle_base |>
  select(
    seq_num,
    cruise,
    station,
    source,
    date,
    year,
    month,
    day,
    latitude,
    longitude,
    depth_actual_m,
    diatom,
    dinoflagellate,
    silicoflagellate,
    prymnesiophyte
  ) |>
  pivot_longer(
    cols = c(diatom, dinoflagellate, silicoflagellate, prymnesiophyte),
    names_to = "taxon_raw",
    values_to = "abundance"
  ) |>
  filter(!is.na(abundance)) |>
  mutate(
    taxon = recode(
      taxon_raw,
      "diatom" = "Bacillariophyta",
      "dinoflagellate" = "Dinoflagellata",
      "silicoflagellate" = "Dictyochophyceae",
      "prymnesiophyte" = "Prymnesiophyceae"
    ),
    internal_id = paste(seq_num, cruise, station, sep = "-"),
    depth = depth_actual_m,
    taxon = taxon |> str_to_lower() |> str_squish(),
    s_count_method = NA_character_,
    s_enumeration_method = "Microscopy",
    s_sampling_method = "niskin bottle",
    s_concentration_method = NA_character_,
    s_microscopy_method = NA_character_,
    s_light_microscope = NA_character_,
    s_preservation = NA_character_,
    s_magnification = NA_character_,
    s_counting_chamber = NA_character_
  )

count_data <- bind_rows(
  count_data_phyto |>
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
      s_light_microscope,
      s_preservation,
      s_magnification,
      s_counting_chamber
    ),
  count_data_bottle |>
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
      s_light_microscope,
      s_preservation,
      s_magnification,
      s_counting_chamber
    )
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
