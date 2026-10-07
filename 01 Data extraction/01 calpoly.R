# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: CalPoly-HAB-MAP
# Original author: Julianne Jager
# Fernando Rayón Viña
# last update: 20250617T16:30 AST
# Original version: 2022-12-12 Julianne Jager
#
# CalPoly - HAB

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(readxl)
library(janitor)
library(this.path)
library(DBI)
library(RSQLite)
library(lubridate)

# LOAD PATHS ────────────────────────────────────────────────────────────────

script_dir <- this.dir()

setwd(script_dir)

fn1 <- file.path(script_dir, "HABs-CalPoly_d103_1542_881a.csv")

# PREPARATION ───────────────────────────────────────────────────────────────

d <- read_csv(fn1, show_col_types = FALSE) |>
  slice(-1) |>
  clean_names() |>
  mutate(
    time = str_remove(time, "Z"),
    datetime_utc = as.POSIXct(time, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC"),
    day = day(datetime_utc),
    month = month(datetime_utc),
    year = year(datetime_utc),
    source = "calpoly"
  ) |>
  clean_names() |>
  rename(
    temperature_C = temp,
    depth_m = depth,
    phosphate_uM = phosphate,
    silicate_uM = silicate,
    nitrite_uM = nitrite,
    nitrate_uM = nitrate
  ) |>
  select(
    -c(
      time,
      air_temp,
      salinity,
      location_code,
      cochlodinium,
      other_dinoflagellates,
      t_da,
      d_da
    )
  )

# COUNT DATA ────────────────────────────────────────────────────────────────

count_data <- d |>
  select(-c(temperature_C:volume_settled_for_counting)) |>
  pivot_longer(
    cols = akashiwo_sanguinea:total_phytoplankton,
    names_to = "taxon",
    values_to = "abundance_l"
  ) |>
  filter(!is.na(abundance_l)) |>
  mutate(
    taxon_clean = taxon |>
      str_replace_all("\\.", "") |>
      str_replace_all("_", " ") |>
      str_replace_all("pseudo nitzschia", "pseudo-nitzschia") |>
      str_squish(),
    count_method = "utermohl"
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
