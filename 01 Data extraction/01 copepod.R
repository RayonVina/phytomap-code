# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: COPEPOD
# Original author: Sadra Dehghani
# Fernando Rayón Viña
# last update: 20250519T00:00 ADT
#
# COPEPOD (Coastal & Oceanic Plankton Ecology, Production & Observation Database)
# NOAA National Marine Fisheries Service (NMFS)
# Global phytoplankton compilation (taxa atlas 2000000)
# Source: https://www.st.nmfs.noaa.gov/copepod/atlas/html/taxatlas_2000000.html
# Metadata: https://www.st.nmfs.noaa.gov/copepod/documentation/short-format_description.html

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(lubridate)
library(janitor)
library(DBI)
library(RSQLite)
library(this.path)

# LOAD PATHS ────────────────────────────────────────────────────────────────

script_dir <- this.dir()

setwd(script_dir)

# READ DATASET FROM SOURCE ──────────────────────────────────────────────────

url <- "https://www.st.nmfs.noaa.gov/copepod/atlas/data_src/copepod__2000000-compilation.txt"

header_line <- read_lines(url, skip = 3, n_max = 1)

col_names <- strsplit(gsub("^#", "", header_line), split = ",")[[1]]

col_names <- trimws(col_names)

# IMPORT DATA ───────────────────────────────────────────────────────────────

fn_a <- read_delim(
  file = url,
  delim = ",",
  skip = 5,
  col_names = col_names,
  col_types = cols()
)

# CLEAN NA SENTINELS ────────────────────────────────────────────────────────

fn_a[] <- lapply(fn_a, function(x) {
  x <- str_trim(x)
  x[x %in% c("n/a", "null", "-9", "-999", "-9999", "")] <- NA
  x
})

# FILTER AND TRANSFORM ──────────────────────────────────────────────────────

count_data <- fn_a |>
  filter(`UNITS...25` %in% c("#/L", "#/ml")) |>
  mutate(
    # ── identifiers ───────────────────────────────────────────────────────────
    source = "copepod",
    subsource = `DATASET-ID` |>
      str_to_lower() |>
      str_replace_all("[^a-z0-9]+", "-") |>
      str_remove_all("^-|-$"),
    internal_id = `SHP-CRUISE`,

    # ── date / time ───────────────────────────────────────────────────────────
    year = as.integer(YEAR),
    month = as.integer(MON),
    day = as.integer(DAY),
    time_raw = as.numeric(TIMEgmt),
    time_valid = if_else(
      is.na(time_raw) | time_raw < 0 | time_raw >= 24 | time_raw == 99.990,
      0,
      time_raw
    ),
    datetime = format(
      as.POSIXct(
        paste(year, month, day, sep = "-"),
        format = "%Y-%m-%d",
        tz = "UTC"
      ) +
        hms::hms(seconds = round(time_valid * 3600)),
      "%Y-%m-%dT%H:%M:%SZ"
    ),

    # ── geography / depth ─────────────────────────────────────────────────────
    latitude = as.numeric(LATITUDE),
    longitude = as.numeric(LONGITDE),
    depth_upper = as.numeric(UPPER_Z),
    depth_lower = as.numeric(LOWER_Z),
    depth = (depth_upper + depth_lower) / 2,
    depth_label = case_when(
      depth_upper == depth_lower ~ "exact",
      !is.na(depth_upper) & !is.na(depth_lower) ~ "range",
      TRUE ~ NA_character_
    ),
    depth_label_min = if_else(depth_label == "range", depth_upper, NA_real_),
    depth_label_max = if_else(depth_label == "range", depth_lower, NA_real_),

    # ── abundance ─────────────────────────────────────────────────────────────

    value_per_volume = as.numeric(`VALUE-per-volu`),
    abundance = case_when(
      `UNITS...25` == "#/ml" ~ value_per_volume * 1000,
      TRUE ~ value_per_volume
    ),

    # ── taxon ─────────────────────────────────────────────────────────────────
    taxon = `SCIENTIFIC NAME -[ modifiers ]-` |>
      str_remove(" -\\[.*\\]-$") |>
      str_trim(),

    # ── sampling metadata ─────────────────────────────────────────────────────
    tow_code = `T`,
    gear_raw = GEAR,
    water_strained_m3 = suppressWarnings(parse_number(`Water Strained`)),
    s_volume_sampled = water_strained_m3 * 1000, # m³ → L

    flowmeter_flag = suppressWarnings(as.integer(`F1...26`)),
    n_notes = case_when(
      flowmeter_flag == 1L ~ "s_volume_sampled measured by flowmeter",
      flowmeter_flag == 0L ~ "s_volume_sampled not measured by flowmeter",
      TRUE ~ NA_character_
    ),

    s_count_method = "cell count (abundance)",
    s_enumeration_method = NA_character_,
    s_microscopy_method = NA_character_,
    s_light_microscope = NA_character_,
    s_preservation = NA_character_,
    s_concentration_method = NA_character_,
    s_volume_concentrated = NA_real_,
    s_volume_counted = NA_real_,
    s_magnification = NA_character_,
    s_counting_chamber = NA_character_
  ) |>
  # Refine s_sampling_method using tow_code and gear_raw ----------------------
  left_join(
    tribble(
      ~tow_code , ~s_sampling_method_base ,
      "B"       , "bottle"                ,
      "V"       , "net (vertical tow)"    ,
      "O"       , "net (oblique tow)"     ,
      "H"       , "net (horizontal tow)"  ,
      "S"       , "net (surface tow)"
    ),
    by = "tow_code"
  ) |>
  mutate(
    s_sampling_method = case_when(
      tow_code == "B" &
        str_detect(str_to_lower(gear_raw), "^ni") ~ "niskin bottle",
      tow_code == "B" &
        str_detect(str_to_lower(gear_raw), "^go") ~ "go-flo bottle",
      tow_code == "B" &
        str_detect(str_to_lower(gear_raw), "^ro") ~ "ctd rosette",
      tow_code == "B" ~ "bottle",
      TRUE ~ s_sampling_method_base
    )
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
    s_microscopy_method,
    s_light_microscope,
    s_preservation,
    s_volume_concentrated,
    s_volume_counted,
    s_magnification,
    s_counting_chamber
  )

sampling_methods_lookup <- tribble(
  ~subsource                                    , ~s_sampling_method_override          , ~s_preservation , ~s_enumeration_method , ~n_notes_extra                                                                        ,

  # Bottle samples
  "meteor-collection"                           , "bottle"                             , NA_character_   , "light microscopy"    , "Phytoplankton abundance, METEOR cruises 1929-1935 (de-04101)"                        ,
  "atlantis-i"                                  , "bottle"                             , NA_character_   , "light microscopy"    , "Total phytoplankton abundance (us-04201); no species-level data"                     ,
  "scor-discoverer-expedition-1970"             , "bottle"                             , NA_character_   , "light microscopy"    , "Grouped counts: diatoms / dinoflagellates / coccolithophores only (us-01041)"        ,

  # Continuous Plankton Recorder (CPR)
  "sahfos-cpr-atlantic-ocean"                   , "continuous plankton recorder (CPR)" , "formalin"      , "light microscopy"    , "CPR survey (SAHFOS); silk mesh ~270 um; preservation in formalin"                    ,
  "ecomon-soop-gulf-of-maine"                   , "continuous plankton recorder (CPR)" , "formalin"      , "light microscopy"    , "CPR survey (EcoMon-SOOP, us-05103); preservation in formalin"                        ,
  "ecomon-soop-mid-atlantic-bight"              , "continuous plankton recorder (CPR)" , "formalin"      , "light microscopy"    , "CPR survey (EcoMon-SOOP, us-05104); preservation in formalin"                        ,

  # Vertical net tows with species-level microscopy — IBSS/Pelagic Ecosystems series
  "pelagic-ecosystems-of-the-tropical-atlantic" , NA_character_                        , "formalin"      , "light microscopy"    , "IBSS Ukraine (ru-05105); phyto + zoo species-level net samples"                      ,
  "pelagic-ecosystems-of-the-indian-ocean"      , NA_character_                        , "formalin"      , "light microscopy"    , "IBSS Ukraine (ru-05106); phyto + zoo species-level net samples"                      ,
  "pelagic-ecosystems-of-the-mediterranean"     , NA_character_                        , "formalin"      , "light microscopy"    , "IBSS Ukraine (ru-05107); phyto + zoo species-level net samples"                      ,
  "atlantniro-plankton"                         , NA_character_                        , "formalin"      , "light microscopy"    , "AtlantNIRO (ru-05301); phyto + zoo species-level, Antarctic"                         ,
  "biological-atlas-of-the-arctic-seas-2000"    , NA_character_                        , "formalin"      , "light microscopy"    , "Barents & Kara Seas atlas (ru-03101); phyto species-level"                           ,
  "hakuho-maru-collection"                      , NA_character_                        , "formalin"      , "light microscopy"    , "HAKUHO MARU cruises (jp-04301); phyto species-level"                                 ,
  "koyo-maru-brazil"                            , NA_character_                        , "formalin"      , "light microscopy"    , "KOYO MARU off N. Brazil (jp-04101); phyto species-level"                             ,

  # Argentine coast / Drake Passage — settled volume + diatom spp.
  "tridente-iii"                                , NA_character_                        , "formalin"      , "light microscopy"    , "Argentine Coast / Drake Passage (ar-01001); settled volume + diatom spp."            ,

  # USCG CHELAN — phytoplankton functional groups
  "uscg-chelan"                                 , NA_character_                        , NA_character_   , "light microscopy"    , "USCG CHELAN (us-01055); copepod spp. + phyto functional groups, N. Pacific / Arctic"
)

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data <- count_data |>
  left_join(sampling_methods_lookup, by = "subsource") |>
  mutate(
    s_sampling_method = if_else(
      !is.na(s_sampling_method_override),
      s_sampling_method_override,
      s_sampling_method
    ),
    s_preservation = if_else(
      is.na(s_preservation.x),
      s_preservation.y,
      s_preservation.x
    ),
    s_enumeration_method = if_else(
      is.na(s_enumeration_method.x),
      s_enumeration_method.y,
      s_enumeration_method.x
    ),
    n_notes = case_when(
      !is.na(n_notes) & !is.na(n_notes_extra) ~ paste(
        n_notes,
        n_notes_extra,
        sep = "; "
      ),
      is.na(n_notes) ~ n_notes_extra,
      TRUE ~ n_notes
    )
  ) |>
  select(
    -s_sampling_method_override,
    -s_preservation.x,
    -s_preservation.y,
    -s_enumeration_method.x,
    -s_enumeration_method.y,
    -n_notes_extra
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
