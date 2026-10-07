# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Sal-2013
# Fernando Rayón Viña
# last update: 20260428T00:00 ADT
#
# Sal, S., López-Urrutia, Á., Irigoien, X., Harbour, D.S. and Harris, R.P. (2013)
# Marine microplankton diversity database. Ecology 94:1658.
# https://doi.org/10.1890/13-0236.1
# 788 stations, 736 unique taxa (1335 cols in Table3: some taxa repeated at
# different carbon content entries in Table2). Abundance in cells/mL → cells/L.
# Multiple oceanographic cruises, 1992–2002, global coverage.
# All counts by D.S. Harbour; Utermöhl inverted light microscopy.
# Zero abundances retained as true absences.

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

fn_stations <- file.path(script_dir, "Table1-stations.csv")

fn_taxa <- file.path(script_dir, "Table2-taxa.csv")

fn_counts <- file.path(script_dir, "Table3-counts.csv")

# IMPORT DATA ───────────────────────────────────────────────────────────────

stations_raw <- read.csv(fn_stations, na.strings = c("NA", "NaN", ""))

taxa_raw <- read.csv(fn_taxa, na.strings = c("NA", "NaN", ""))

counts_raw <- read.csv(
  fn_counts,
  na.strings = c("NA", "NaN", ""),
  check.names = FALSE
)

# TRANSFORM: STATIONS ───────────────────────────────────────────────────────

stations <- stations_raw |>
  rename(
    sample_id = SampleID,
    cruise = Cruise,
    latitude = Lat,
    longitude = Lon,
    depth = Depth,
    temperature_C = Temperature,
    chla_ug_l = Chl,
    nitrate_uM = Nitrate,
    nitrite_uM = Nitrite,
    ammonium_uM = Ammonium,
    phosphate_uM = Phosphate,
    silicate_uM = Silicate,
    mld_m = MLD,
    surface_par = SurfacePAR,
    parz = PARz,
    kd490 = Kd490,
    daylength_h = Daylength
  ) |>
  separate(Date, into = c("day", "month", "year"), sep = "-", convert = TRUE) |>
  mutate(
    date = make_date(year, month, day),
    datetime = format(
      as.POSIXct(paste(date, "00:00"), tz = "UTC"),
      "%Y-%m-%dT%H:%M:%SZ"
    )
  )

# TRANSFORM: TAXA ───────────────────────────────────────────────────────────

taxa <- taxa_raw |>
  clean_names() |>
  mutate(
    col_index = row_number(),
    taxon_canonical = case_when(
      !is.na(genus) & !is.na(species) ~ paste(genus, species),
      !is.na(genus) ~ genus,
      !is.na(author_comments) ~ paste(group, author_comments),
      TRUE ~ group
    ) |>
      str_to_lower() |>
      str_replace_all("\\.", "") |>
      str_squish()
  )

# TRANSFORM: COUNTS ─────────────────────────────────────────────────────────

colnames(counts_raw)[1] <- "sample_id"

n_taxa <- ncol(counts_raw) - 1

counts_indexed <- counts_raw

colnames(counts_indexed)[-1] <- as.character(seq_len(n_taxa))

count_long <- counts_indexed |>
  pivot_longer(
    cols = -sample_id,
    names_to = "col_index",
    values_to = "abundance_ml",
    names_transform = list(col_index = as.integer)
  ) |>
  left_join(
    taxa |> select(col_index, taxon_canonical),
    by = "col_index"
  ) |>
  # Sum abundances where same taxon appears in multiple columns (different carbon entries)
  group_by(sample_id, taxon_canonical) |>
  summarise(abundance_ml = sum(abundance_ml, na.rm = TRUE), .groups = "drop") |>
  rename(taxon = taxon_canonical)

time_location <- stations |>
  select(
    sample_id,
    cruise,
    latitude,
    longitude,
    year,
    month,
    day,
    datetime,
    depth
  )

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data <- count_long |>
  left_join(time_location, by = "sample_id") |>
  filter(!is.na(latitude), !is.na(longitude)) |>
  mutate(
    abundance = abundance_ml * 1000, # cells/mL → cells/L
    source = "sal-2013",
    subsource = str_to_lower(cruise),
    internal_id = paste(sample_id, sprintf("%06d", row_number()), sep = "-"),
    depth_label = "exact",
    s_count_method = "utermohl",
    s_enumeration_method = "light microscopy",
    s_sampling_method = "CTD niskin bottle",
    s_concentration_method = "sedimentation",
    s_microscopy_method = "light microscope",
    s_light_microscope = "inverted",
    s_preservation = "acid-lugol's iodine (1-5%) + formalin",
    s_counting_chamber = "utermohl"
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
    s_count_method,
    s_enumeration_method,
    s_sampling_method,
    s_concentration_method,
    s_microscopy_method,
    s_light_microscope,
    s_preservation,
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
