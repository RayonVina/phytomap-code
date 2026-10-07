# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Vendrick_2026_FRV
# Fernando Rayón Viña
# last update: 20260319T15:00 AST
#
# Venrick - California Current Phytoplankton 1996-2022
# CCE LTER dataset 254 / Venrick (2002, 2012, in prep.)

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(readxl)
library(lubridate)
library(DBI)
library(RSQLite)
library(this.path)

# DATA SOURCE ───────────────────────────────────────────────────────────────

script_dir <- this.dir()

# LOAD PATHS ────────────────────────────────────────────────────────────────

setwd(script_dir)

fn_def <- file.path(script_dir, "Species_abundances_definitions.xlsx")

fn_ab1 <- file.path(script_dir, "Species_Abundances_1996-2012.xlsx")

fn_ab2 <- file.path(script_dir, "Species_abundances_2012-2018.xlsx")

fn_ab3 <- file.path(script_dir, "Species-Abundances_2019-2022.xlsx")

# TAXONOMY ──────────────────────────────────────────────────────────────────

taxonomy <- read_excel(fn_def, sheet = "Species Codes") |>
  filter(!is.na(Code)) |>
  select(Code, Taxa, Species) |>
  rename(species_code = Code, taxa_group = Taxa, species_name = Species) |>
  mutate(
    taxon = if_else(
      !is.na(species_name) & species_name != "",
      species_name,
      taxa_group
    ) |>
      str_to_lower() |>
      str_replace_all("\\.", "") |>
      str_squish()
  )

# REGION COORDINATES ────────────────────────────────────────────────────────

region_coords <- tibble(
  region = c("NE", "SE", "Alley", "Offshore"),
  latitude = c(32.5, 31.0, 33.0, 33.3),
  longitude = c(-119.3, -119.1, -121.8, -124.0)
)

# IMPORT DATA ───────────────────────────────────────────────────────────────

read_venrick_sheet <- function(file, sheet_name) {
  raw <- read_excel(
    file,
    sheet = sheet_name,
    col_names = TRUE,
    col_types = "text"
  )

  r_names <- colnames(raw)
  regions <- as.character(raw[1, ])

  # Build colnames: fill cruise label forward across groups of 4 regions
  new_names <- character(length(r_names))
  new_names[1] <- "species_code"
  cruise_current <- NA_character_

  for (i in seq_along(r_names)[-1]) {
    if (!startsWith(r_names[i], "...")) {
      cruise_current <- r_names[i]
    }
    new_names[i] <- paste(cruise_current, regions[i], sep = "___")
  }

  colnames(raw) <- new_names
  raw[-1, ] |>
    pivot_longer(
      -species_code,
      names_to = "cruise_region",
      values_to = "abundance_raw"
    ) |>
    separate(cruise_region, into = c("cruise_label", "region"), sep = "___") |>
    mutate(
      species_code = suppressWarnings(as.numeric(species_code)),
      abundance_raw = suppressWarnings(as.numeric(abundance_raw))
    )
}

# READ ALL FILES ────────────────────────────────────────────────────────────

read_venrick_file <- function(file, exclude_sheets = NULL) {
  sheets <- excel_sheets(file)
  sheets <- sheets[!str_detect(str_to_lower(sheets), "^sheet")]
  if (!is.null(exclude_sheets)) {
    sheets <- sheets[!str_to_lower(sheets) %in% str_to_lower(exclude_sheets)]
  }
  bind_rows(lapply(sheets, \(s) read_venrick_sheet(file, s)))
}

ab_raw <- bind_rows(
  read_venrick_file(fn_ab1, exclude_sheets = "CC12DATA"),
  read_venrick_file(fn_ab2),
  read_venrick_file(fn_ab3)
)

# PARSE CRUISE LABEL → YEAR + MONTH ─────────────────────────────────────────

ab_parsed <- ab_raw |>
  mutate(
    yymm = str_extract(cruise_label, "\\d{4}"),
    yy = as.integer(str_sub(yymm, 1, 2)),
    mm = as.integer(str_sub(yymm, 3, 4)),
    year = if_else(yy >= 96L, 1900L + yy, 2000L + yy),
    month = mm
  ) |>
  mutate(month = if_else(cruise_label == "CalCOFI 0704...2", 1L, month))

# HANDLE SPECIAL VALUES     ─ ───────────────────────────────────────────────

ab_clean <- ab_parsed |>
  mutate(
    abundance = case_when(
      abundance_raw %in% c(-13, -1) ~ NA_real_,
      TRUE ~ abundance_raw
    )
  ) |>
  filter(!is.na(abundance))

# COUNT DATA ────────────────────────────────────────────────────────────────

count_data <- ab_clean |>
  left_join(taxonomy, by = "species_code") |>
  left_join(region_coords, by = "region") |>
  filter(!is.na(taxon)) |>
  mutate(
    abundance = abundance * 10, # cells/100mL → cells/L
    depth = 10, # nominal: second depth 5-15m (abstract)
    day = 15L,
    internal_id = paste(
      "venrick",
      str_replace_all(cruise_label, " ", "-"),
      region,
      species_code,
      sep = "_"
    ),
    notes = "day set to 15 (monthly midpoint convention)",
    source = "venrick-2018",
    s_count_method = "Utermohl",
    s_enumeration_method = "Microscopy",
    s_sampling_method = "niskin bottle",
    s_concentration_method = "sedimentation",
    s_microscopy_method = "light microscope",
    s_light_microscope = "inverted",
    s_preservation = "formalin 1%",
    s_magnification = "100x-250x",
    s_counting_chamber = "utermohl"
  ) |>
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
    s_counting_chamber,
    notes
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
