# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Wilfred-Allen
# Fernando Rayón Viña
# last update: 20260429T00:00 ADT
#
# W. E. Allen's Phytoplankton Species Time Series
# North American Pacific Coast, 1917-1939
# Hewes, C.D. and Thomas, W.H. (2002). W. E. Allen's Phytoplankton Species
# Time Series and Environmental Factors from the North American Pacific Coast,
# With Related Illustrations from Cupp and Kofoid.
# Scripps Institution of Oceanography Technical Report. December 2002.
# https://escholarship.org/uc/item/39v2v8w5

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

data_dir <- file.path(script_dir, "AllenData")

# STATION METADATA ──────────────────────────────────────────────────────────

stations <- tribble(
  ~station_id , ~station_name                                  , ~latitude , ~longitude ,
  "SIO"       , "Scripps Institution of Oceanography Pier, CA" , 32.8667   , -117.2500  ,
  "OP"        , "Oceanside Pier, CA"                           , 33.1833   , -117.3833  ,
  "PH"        , "Point Hueneme Pier, CA"                       , 34.1500   , -119.2000  ,
  "PG"        , "Pacific Grove Pier, CA"                       , 36.6333   , -121.9167  ,
  "FI"        , "Farallon Islands, CA"                         , 37.4167   , -122.6000  ,
  "SC"        , "Scotch Cap, AK"                               , 54.4000   , -164.7500
)

# IMPORT DATA ───────────────────────────────────────────────────────────────

lf_phyto <- list.files(
  file.path(data_dir, "Totals"),
  pattern = "(diat|dino)1?\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

# STATION METADATA ──────────────────────────────────────────────────────────

station_order <- c("FI", "OP", "PG", "PH", "SC", "SIO")

# IMPORT DATA ───────────────────────────────────────────────────────────────

read_phyto_file <- function(path, station_id, taxon_group) {
  df <- read_csv(
    path,
    show_col_types = FALSE,
    locale = locale(encoding = "latin1")
  ) |>
    clean_names()
  names(df)[1:3] <- c("year", "month", "day")
  df |>
    select(-matches("^\\.\\.\\.\\d+$")) |>
    pivot_longer(
      cols = -c(year, month, day),
      names_to = "taxon_raw",
      values_to = "abundance"
    ) |>
    mutate(station_id = station_id, taxon_group = taxon_group)
}

phyto_raw <- map2_dfr(
  seq_along(station_order),
  station_order,
  function(i, sid) {
    bind_rows(
      read_phyto_file(lf_phyto[2 * i - 1], sid, "diatom"),
      read_phyto_file(lf_phyto[2 * i], sid, "dinoflagellate")
    )
  }
)

# TRANSFORM: ABUNDANCE ──────────────────────────────────────────────────────

count_data <- phyto_raw |>
  filter(!is.na(year), !is.na(month)) |>
  mutate(
    year = as.integer(year),
    month = as.integer(month),
    day = as.integer(day),
    date = make_date(year, month, day),
    datetime = format(
      as.POSIXct(paste(date, "00:00"), tz = "UTC"),
      "%Y-%m-%dT%H:%M:%SZ"
    ),
    taxon = taxon_raw |>
      str_replace_all("_", " ") |>
      str_to_lower() |>
      str_remove("\\s+species$") |>
      str_squish(),
    abundance = replace_na(abundance, 0),
    source = "wilfred-allen",
    subsource = case_when(
      station_id == "FI" ~ "farallon-islands",
      station_id == "OP" ~ "oceanside",
      station_id == "PG" ~ "pacific-grove",
      station_id == "PH" ~ "port-hueneme",
      station_id == "SC" ~ "scotch-cap",
      station_id == "SIO" ~ "scripps-pier"
    ),
    internal_id = paste(str_to_lower(station_id), year, month, day, sep = "-"),
    depth = 0,
    depth_label = "surface",
    s_count_method = "sedgwick-rafter",
    s_enumeration_method = "light microscopy",
    s_microscopy_method = "light microscope",
    s_counting_chamber = "sedgwick-rafter",
    s_sampling_method = "pier surface grab",
    s_preservation = "formalin",
    s_volume_counted = 1,
    s_concentration_method = case_when(
      station_id == "SIO" & year >= 1930 ~ "gravity settling (24h, 1L)",
      TRUE ~ "net filtration (40um, 25L)"
    ),
    s_volume_sampled = case_when(
      station_id == "SIO" & year >= 1930 ~ 1000,
      TRUE ~ 25000
    ),
    n_notes = if_else(
      taxon_group == "dinoflagellate" & station_id == "SIO" & year >= 1930,
      "SIO 1930-1939 gravity settling yields higher dinoflagellate counts than net filtration (Allen 1930)",
      NA_character_
    )
  ) |>
  left_join(
    stations |> select(station_id, latitude, longitude),
    by = "station_id"
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
    taxon,
    abundance,
    n_notes,
    s_count_method,
    s_enumeration_method,
    s_sampling_method,
    s_concentration_method,
    s_volume_sampled,
    s_microscopy_method,
    s_preservation,
    s_volume_counted,
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
