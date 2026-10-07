# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Australian-tasmania
# Original author: Sadra Dehghani
# Fernando Rayón Viña
# last update: 20250708T15:45 AST
#
# Australian Tasmania

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(janitor)
library(this.path)
library(DBI)
library(RSQLite)
library(lubridate)

# LOAD PATHS ────────────────────────────────────────────────────────────────

script_dir <- this.dir()

fn1 <- file.path(
  script_dir,
  "Phytoplankton_sampling_in_the_coastal_waters_of_south_eastern_Tasmania.csv"
)

# IMPORT DATA ───────────────────────────────────────────────────────────────

d <- read_csv(fn1, show_col_types = FALSE) |>
  mutate(
    year = year(DATE_TRIP),
    month = month(DATE_TRIP),
    day = day(DATE_TRIP),
    source = "aus-tas",
    count_method = "Sedgewick-Rafter"
  ) |>
  rename(
    latitude = LATITUDE,
    longitude = LONGITUDE,
    datetime_utc = DATE_TRIP,
    depth_m = DEPTH_CODE,
    wind_speed_m = WIND_SPEED,
    depth_secchi_m = DEPTH_SECCHI,
    cell_counted = NUMBER_CELLS_COUNTED,
    biovolume = BIOVOLUME, # unit : ul_cell
    biovolume_per_cell = PER_CELL_BIOVOLUME,
    initial_volume_ml = INITIAL_SAMPLE_VOLUME,
    final_volume_ml = FINAL_SAMPLE_VOLUME,
    taxon = TAXON,
    taxon_group = TAXON_GROUP,
    taxon_trophy = TAXON_TROPHY,
    taxon_clean = GENUS_SPECIES,
    sample_type = SAMPLE_TYPE,
    cell_density_ml = CORRECTED_CELL_CONCENTRATION_CELLS_PER_MILLILITRE,
    wind_direction = WIND_DIRECTION,
    site_code = SITE_CODE,
    trip_id = TRIP_ID,
    abundance_l = CORRECTED_CELL_CONCENTRATION_CELLS_PER_LITRE,
    fid = FID
  ) |>
  select(-c(NUMBER_SQUARES_COUNTED, GEOM))

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data <- d |>
  select(
    fid,
    trip_id,
    site_code,
    latitude,
    longitude,
    datetime_utc,
    year,
    month,
    day,
    taxon,
    taxon_clean,
    depth_m,
    depth_secchi_m,
    biovolume,
    biovolume_per_cell,
    abundance_l,
    count_method,
    source
  ) |>
  mutate(
    depth_m = as.numeric(gsub("m", "", depth_m)),
    biovolume_um3_l = biovolume_per_cell * abundance_l,
    taxon_clean = taxon_clean |>
      str_to_lower() |>
      str_replace_all("\\.", "") |>
      str_squish()
  ) |>
  separate(
    taxon_clean,
    into = c("tax_genus", "tax_species", "tax_epithet"),
    sep = " ",
    extra = "merge",
    fill = "right",
    remove = FALSE
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
