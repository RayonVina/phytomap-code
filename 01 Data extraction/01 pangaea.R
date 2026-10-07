# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Pangaea
# Fernando Rayón Viña
# last update: 20251127T11:50 CET
#
# PANGAEA

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(DBI)
library(RSQLite)
library(lubridate)
library(this.path)

# LOAD PATHS ────────────────────────────────────────────────────────────────

script_dir <- file.path(this.dir())

setwd(script_dir)

all_datasets <- list.files(
  path = "files",
  pattern = "standardized_data\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

# NORMALIZATION ─────────────────────────────────────────────────────────────

normalize_column <- function(col_name) {
  col <- str_trim(col_name)

  standardized <- case_when(
    # Event variations
    col %in% c("event", "Event", "Event ") ~ "event",

    # Comment variations
    col %in% c("comment", "Comment", "Comment 2", "coment") ~ "comment",

    # Temperature variations
    col %in% c("Temp [°C]", "Temp [C]", "temperature_C") ~ "temperature_c",

    # Depth variations (single depth)
    col %in% c("depth", "depth_m") ~ "depth_m",

    # Depth range variations - UNIFIED (top = min, bot = max)
    col %in%
      c(
        "Depth top [m]",
        "Depth top",
        "depth top [m]",
        "depth top",
        "Depth water [m] (min)",
        "depth water [m] (min)"
      ) ~ "depth_top_m",

    col %in%
      c(
        "Depth bot [m]",
        "Depth bot",
        "depth bot [m]",
        "depth bot",
        "Depth water [m] (max)",
        "depth water [m] (max)"
      ) ~ "depth_bot_m",

    # Abundance variations
    col %in% c("cells_per_litre", "abundance_l", "cells") ~ "abundance_l",

    # Elevation variations - UNIFY ALL
    col %in%
      c(
        "elevation",
        "Elevation",
        "Elevation [m]",
        "Elevation [m a.s.l.]"
      ) ~ "elevation",

    # Taxon variations
    col %in% c("original_taxon", "original_taxon ") ~ "original_taxon",
    col %in% c("simplified_taxon", "simpilifed_taxon") ~ "simplified_taxon",
    col == "taxon" ~ "taxon",
    col == "taxon_ex" ~ "taxon_ex",

    # Size fraction variations
    col %in%
      c("size_fraction", "size_fraction ", "fraciton_size") ~ "size_fraction",

    # Sampling variations
    col %in%
      c(
        "sampling_method",
        "Sampling_method",
        "sampling_metod"
      ) ~ "sampling_method",

    # Station variations
    col %in% c("station", "Station") ~ "station",

    # Date/Time variations
    col %in%
      c("Date/Time", "Date/Time ", "Date/Tme", "Date_time") ~ "date_time",

    # Keep as-is (lowercase)
    TRUE ~ tolower(col)
  )

  return(standardized)
}

# IMPORT DATA ───────────────────────────────────────────────────────────────

read_and_normalize <- function(file_path) {
  tryCatch(
    {
      df <- read.csv(file_path, header = FALSE, stringsAsFactors = FALSE)

      headers_raw <- as.character(df[1, ])
      headers_normalized <- normalize_column(headers_raw)

      # Remove columns with empty names
      keep_cols <- headers_normalized != ""
      df <- df[, keep_cols]
      headers_raw <- headers_raw[keep_cols]
      headers_normalized <- headers_normalized[keep_cols]

      # Handle duplicate normalized names
      if (any(duplicated(headers_normalized))) {
        # For abundance_l: prioritize cells_per_litre over cells
        if (
          "abundance_l" %in% headers_normalized[duplicated(headers_normalized)]
        ) {
          abundance_idx <- which(headers_normalized == "abundance_l")
          if (length(abundance_idx) > 1) {
            priority <- which(headers_raw == "cells_per_litre")
            if (length(priority) > 0) {
              keep <- priority[1]
            } else {
              keep <- abundance_idx[1]
            }
            df <- df[, -setdiff(abundance_idx, keep)]
            headers_normalized <- headers_normalized[
              -setdiff(abundance_idx, keep)
            ]
          }
        }

        # For comment: concatenate if both present
        if ("comment" %in% headers_normalized[duplicated(headers_normalized)]) {
          comment_idx <- which(headers_normalized == "comment")
          if (length(comment_idx) > 1) {
            df$comment_combined <- apply(df[, comment_idx], 1, function(x) {
              paste(x[x != ""], collapse = " | ")
            })
            df <- df[, -comment_idx]
            headers_normalized <- headers_normalized[-comment_idx]
            df$comment <- df$comment_combined
            df$comment_combined <- NULL
            headers_normalized <- c(headers_normalized, "comment")
          }
        }

        # For elevation: coalesce (take first non-empty)
        if (
          "elevation" %in% headers_normalized[duplicated(headers_normalized)]
        ) {
          elevation_idx <- which(headers_normalized == "elevation")
          if (length(elevation_idx) > 1) {
            df$elevation_unified <- apply(df[, elevation_idx], 1, function(x) {
              x <- x[x != "" & !is.na(x)]
              if (length(x) > 0) x[1] else NA_character_
            })
            df <- df[, -elevation_idx]
            headers_normalized <- headers_normalized[-elevation_idx]
            df$elevation <- df$elevation_unified
            df$elevation_unified <- NULL
            headers_normalized <- c(headers_normalized, "elevation")
          }
        }

        # For depth_top_m: coalesce (take first non-empty)
        if (
          "depth_top_m" %in% headers_normalized[duplicated(headers_normalized)]
        ) {
          depth_top_idx <- which(headers_normalized == "depth_top_m")
          if (length(depth_top_idx) > 1) {
            df$depth_top_unified <- apply(df[, depth_top_idx], 1, function(x) {
              x <- x[x != "" & !is.na(x)]
              if (length(x) > 0) x[1] else NA_character_
            })
            df <- df[, -depth_top_idx]
            headers_normalized <- headers_normalized[-depth_top_idx]
            df$depth_top_m <- df$depth_top_unified
            df$depth_top_unified <- NULL
            headers_normalized <- c(headers_normalized, "depth_top_m")
          }
        }

        # For depth_bot_m: coalesce (take first non-empty)
        if (
          "depth_bot_m" %in% headers_normalized[duplicated(headers_normalized)]
        ) {
          depth_bot_idx <- which(headers_normalized == "depth_bot_m")
          if (length(depth_bot_idx) > 1) {
            df$depth_bot_unified <- apply(df[, depth_bot_idx], 1, function(x) {
              x <- x[x != "" & !is.na(x)]
              if (length(x) > 0) x[1] else NA_character_
            })
            df <- df[, -depth_bot_idx]
            headers_normalized <- headers_normalized[-depth_bot_idx]
            df$depth_bot_m <- df$depth_bot_unified
            df$depth_bot_unified <- NULL
            headers_normalized <- c(headers_normalized, "depth_bot_m")
          }
        }
      }

      colnames(df) <- headers_normalized
      df <- df[-1, ]

      df$source_doi <- basename(dirname(file_path))

      return(df)
    },
    error = function(e) {
      warning(paste("Failed to read:", file_path, "-", e$message))
      return(NULL)
    }
  )
}

all_data_list <- map(all_datasets, read_and_normalize)

read_pangaea_628513 <- function(tab_path, doi_folder) {
  lines <- readLines(
    tab_path,
    encoding = "UTF-8"
  )

  metadata_end <- which(str_detect(lines, "^\\*/"))[1]

  metadata_lines <- lines[2:(metadata_end - 1)]

  header_line <- lines[metadata_end + 1]

  col_names <- str_split(header_line, "\t")[[1]]

  data_lines <- lines[(metadata_end + 2):length(lines)]

  data_wide <- data_lines |>
    map(~ str_split(.x, "\t")[[1]]) |>
    (function(x) do.call(rbind, x))() |>
    as_tibble(.name_repair = "minimal") |>
    setNames(col_names) |>
    mutate(across(everything(), ~ na_if(.x, ""))) |>
    type_convert(guess_integer = TRUE)

  metadata_text <- paste(metadata_lines, collapse = "\n")

  param_pattern <- "\\t([^\\[]+) \\[#/l\\] \\(([^\\)]+)\\)"

  matches <- str_match_all(metadata_text, param_pattern)[[1]]

  taxon_mapping <- tibble(
    taxon = str_trim(matches[, 3]), # Abbreviation
    taxon_ex = str_trim(matches[, 2]) # Full name
  ) |>
    filter(!is.na(taxon), !is.na(taxon_ex))

  data_long <- data_wide |>
    pivot_longer(
      cols = -c(Latitude, Longitude, No, `Bathy depth [m]`),
      names_to = "taxon",
      values_to = "abundance_l"
    ) |>
    mutate(abundance_l = as.numeric(abundance_l)) |>
    filter(!is.na(abundance_l), abundance_l > 0) |>
    rename(
      latitude = Latitude,
      longitude = Longitude,
      depth_m = `Bathy depth [m]`
    ) |>
    mutate(taxon = str_remove(taxon, " \\[#/l\\]$")) |>
    left_join(taxon_mapping, by = "taxon") |>
    mutate(
      taxon_ex = str_remove(taxon_ex, " sp\\.$"),
      original_taxon = coalesce(taxon_ex, taxon),
      simplified_taxon = coalesce(taxon_ex, taxon),
      event = "Sechura_2007-01",
      year = 2007L,
      month = 1L,
      day = NA_integer_,
      time = NA_character_,
      counting_method = "Utermöhl (1958)",
      sampling_method = "Bottom water sample",
      elevation = 3.0,
      temperature_C = NA_real_,
      presence_absence = NA_character_,
      size_fraction = NA_character_,
      heterotrophic = NA_character_,
      mixotrophic = NA_character_
    ) |>
    select(
      event,
      latitude,
      longitude,
      year,
      time,
      month,
      day,
      depth_m,
      taxon,
      taxon_ex,
      counting_method,
      abundance_l,
      presence_absence,
      size_fraction,
      sampling_method,
      elevation,
      temperature_C,
      heterotrophic,
      mixotrophic
    )

  data_long |>
    transmute(
      source_doi = doi_folder,
      event,
      latitude,
      longitude,
      year,
      month,
      day,
      time,
      depth_m,
      original_taxon = coalesce(taxon_ex, taxon),
      counting_method,
      abundance_l,
      presence_absence,
      sampling_method
    ) |>
    mutate(across(everything(), as.character))
}

read_pangaea_849248 <- function(tab_path, doi_folder) {
  lines <- readLines(
    tab_path,
    encoding = "UTF-8"
  )

  metadata_end <- which(str_detect(lines, "^\\*/"))[1]

  metadata_lines <- lines[2:(metadata_end - 1)]

  header_line <- lines[metadata_end + 1]

  col_names <- str_split(header_line, "\t")[[1]]

  data_lines <- lines[(metadata_end + 2):length(lines)]

  data_wide <- data_lines |>
    map(~ str_split(.x, "\t")[[1]]) |>
    (function(x) do.call(rbind, x))() |>
    as_tibble(.name_repair = "minimal") |>
    setNames(col_names)

  data_wide <- data_wide[, !duplicated(names(data_wide))]

  data_wide <- data_wide |>
    mutate(across(everything(), ~ na_if(.x, ""))) |>
    type_convert(guess_integer = TRUE)

  metadata_text <- paste(metadata_lines, collapse = "\n")

  param_pattern <- "\\t([^\\[]+) \\[#/l\\] \\(([^\\)]+)\\)"

  matches <- str_match_all(metadata_text, param_pattern)[[1]]

  taxon_mapping <- tibble(
    taxon = str_trim(matches[, 3]),
    taxon_ex = str_trim(matches[, 2])
  ) |>
    filter(!is.na(taxon), !is.na(taxon_ex)) |>
    mutate(
      taxon_ex = str_remove(taxon_ex, " spp\\.$"),
      taxon_ex = str_remove(taxon_ex, " sp\\.$"),
      taxon_ex = str_remove(taxon_ex, " var\\. .+$")
    )

  event_pattern <- "\\t([^\\(]+) \\([^\\)]+\\) \\* LATITUDE: ([\\-]?[\\d\\.]+) \\* LONGITUDE: ([\\-]?[\\d\\.]+)"

  event_matches <- str_match_all(metadata_text, event_pattern)[[1]]

  coord_mapping <- tibble(
    event = str_trim(event_matches[, 2]),
    latitude = as.numeric(event_matches[, 3]),
    longitude = as.numeric(event_matches[, 4])
  ) |>
    filter(!is.na(event), !is.na(latitude), !is.na(longitude))

  data_long <- data_wide |>
    pivot_longer(
      cols = -c(Event, `Event 2`, `Date/Time`, Longitude, `Depth water [m]`),
      names_to = "taxon",
      values_to = "abundance_l"
    ) |>
    mutate(abundance_l = as.numeric(abundance_l)) |>
    filter(!is.na(abundance_l), abundance_l > 0) |>
    rename(
      event = Event,
      datetime = `Date/Time`,
      depth_m = `Depth water [m]`
    ) |>
    select(-Longitude) |>
    mutate(
      taxon = str_remove(taxon, " \\[#/l\\]$"),
      datetime = as.POSIXct(datetime, format = "%Y-%m-%dT%H:%M"),
      year = year(datetime),
      month = month(datetime),
      day = day(datetime),
      time = format(datetime, "%H:%M:%S")
    ) |>
    left_join(coord_mapping, by = "event") |>
    left_join(taxon_mapping, by = "taxon") |>
    mutate(
      original_taxon = coalesce(taxon_ex, taxon),
      simplified_taxon = coalesce(taxon_ex, taxon),
      counting_method = "Counting",
      sampling_method = "Bottle, Niskin (NIS)",
      elevation = NA_real_,
      temperature_C = NA_real_,
      presence_absence = NA_character_,
      size_fraction = NA_character_,
      heterotrophic = NA_character_,
      mixotrophic = NA_character_
    ) |>
    select(
      event,
      latitude,
      longitude,
      year,
      time,
      month,
      day,
      depth_m,
      taxon_ex,
      counting_method,
      abundance_l,
      presence_absence,
      size_fraction,
      sampling_method
    )

  data_long |>
    transmute(
      source_doi = doi_folder,
      event,
      latitude,
      longitude,
      year,
      month,
      day,
      time,
      depth_m,
      original_taxon = taxon_ex,
      counting_method,
      abundance_l,
      presence_absence,
      sampling_method
    ) |>
    mutate(across(everything(), as.character))
}


raw_inputs <- tribble(
  ~doi_folder                                    , ~filename , ~reader ,
  "10.1594_PANGAEA.628513"                       ,
  "Sechura_phytoplankton_sea_bottom_2007-01.tab" ,
  list(read_pangaea_628513)                      ,
  "10.1594_PANGAEA.849248"                       ,
  "Bilim2_200803_phyto.tab"                      ,
  list(read_pangaea_849248)
)

for (i in seq_len(nrow(raw_inputs))) {
  doi_folder <- raw_inputs$doi_folder[i]
  candidates <- c(
    file.path(script_dir, "files", doi_folder, raw_inputs$filename[i]),
    file.path(script_dir, doi_folder, raw_inputs$filename[i])
  )
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates)) {
    keep <- basename(dirname(all_datasets)) != doi_folder
    all_datasets <- all_datasets[keep]
    all_data_list <- all_data_list[keep]
  }
}

for (i in seq_len(nrow(raw_inputs))) {
  doi_folder <- raw_inputs$doi_folder[i]
  candidates <- c(
    file.path(script_dir, "files", doi_folder, raw_inputs$filename[i]),
    file.path(script_dir, doi_folder, raw_inputs$filename[i])
  )
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates)) {
    reader <- raw_inputs$reader[[i]][[1]]
    all_data_list <- append(
      all_data_list,
      list(reader(candidates[1], doi_folder))
    )
  }
}

optional_input_fields <- c(
  "depth_m",
  "depth_top_m",
  "depth_bot_m",
  "abundance_l",
  "latitude",
  "longitude",
  "year",
  "month",
  "day",
  "presence_absence",
  "time",
  "simplified_taxon",
  "comment",
  "size_fraction",
  "elevation",
  "picopl prod c [µg/l/day]",
  "nanopl prod c [µg/l/day]",
  "phytopl prod c [µg/l/day]",
  "short name",
  "principal investigator",
  "name_comment",
  "temperature_c",
  "phytopl v [ml/m**3]",
  "heterotrophic",
  "taxon",
  "mixotrophic",
  "date/time end",
  "date/time start",
  "coeff",
  "bottle",
  "biomass_unit",
  "depth therm [m] (lower boundary)",
  "depth therm [m] (upper boundary)",
  "depth count [m]",
  "biomass",
  "biomass [µg/l]",
  "biomass (mg/m**3)",
  "autotrophic",
  "bathy depth [m]",
  "bottle [#]",
  "par [%]",
  "source_doi",
  "file name",
  "date_time",
  "station",
  "date",
  "samp vol [l]",
  "vol [m**3]",
  "samp vol [m**3] (computed/converted)",
  "samp vol [ml]",
  "samp vol [ml] (computed/converted)",
  "sampled_volume",
  "niskin volume",
  "sample id",
  "sample label",
  "event",
  "event 2",
  "campaign",
  "original_taxon",
  "taxon_ex",
  "sample_comment",
  "sample comment",
  "method/device",
  "counting_method",
  "sampling_method",
  "preservative"
)
all_data_list <- map(all_data_list, function(df) {
  if (is.null(df)) {
    return(NULL)
  }
  for (field in setdiff(optional_input_fields, names(df))) {
    df[[field]] <- rep(NA_character_, nrow(df))
  }
  df
})


all_data_list <- compact(all_data_list)

combined_data <- bind_rows(all_data_list)

problematic_cols <- combined_data |>
  select(any_of(c("", "na"))) |>
  colnames()

if (length(problematic_cols) > 0) {
  combined_data <- combined_data |> select(-any_of(problematic_cols))
}

# DATA TRANSFORMATION ───────────────────────────────────────────────────────

count_data <- combined_data |>
  mutate(
    depth_m = as.numeric(depth_m),
    depth_top_m = as.numeric(depth_top_m),
    depth_bot_m = as.numeric(depth_bot_m),
    abundance_l = as.numeric(abundance_l),
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude),
    year = as.integer(year),
    month = as.integer(month),
    day = as.integer(day),
    source = "pangaea",
    subsource = source_doi,
    presence_absence = as.logical(as.integer(presence_absence)),
    across(where(is.character), ~ na_if(.x, "")),
    time = case_when(
      is.na(time) ~ NA_character_,
      str_detect(time, "^\\d{1,2}:\\d{2}$") ~ paste0(time, ":00"),
      TRUE ~ time
    )
  ) |>
  select(
    -c(
      simplified_taxon,
      comment,
      size_fraction,
      elevation,
      `picopl prod c [µg/l/day]`,
      `nanopl prod c [µg/l/day]`,
      `phytopl prod c [µg/l/day]`,
      `short name`,
      `principal investigator`,
      name_comment,
      temperature_c,
      `phytopl v [ml/m**3]`,
      heterotrophic,
      taxon,
      mixotrophic,
      `date/time end`,
      `date/time start`,
      coeff,
      bottle,
      biomass_unit,
      `depth therm [m] (lower boundary)`,
      `depth therm [m] (upper boundary)`,
      `depth count [m]`,
      biomass,
      `biomass [µg/l]`,
      `biomass (mg/m**3)`,
      autotrophic,
      `bathy depth [m]`,
      `bottle [#]`,
      `par [%]`,
      source_doi,
      `file name`
    )
  )

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data <- count_data |>
  mutate(
    day_imputed = is.na(day) & is.na(date_time) & !is.na(year) & !is.na(month),
    day = if_else(day_imputed, 15L, day),
    notes = if_else(
      day_imputed,
      "Day set to 15 by convention",
      NA_character_
    ),
    day_imputed = NULL,
    datetime_posix = coalesce(
      ymd_hms(date_time, quiet = TRUE),
      ymd_hms(
        paste(
          paste(year, month, day, sep = "-"),
          coalesce(time, "00:00:00")
        ),
        quiet = TRUE
      )
    ),
    time = NULL,
    datetime = if_else(
      is.na(datetime_posix),
      NA_character_,
      format(datetime_posix, "%Y-%m-%dT%H:%M:%S")
    ),
    datetime_posix = NULL,
    station_tag = if_else(
      is.na(station),
      NA_character_,
      paste0("st: ", station)
    ),
    date_time = NULL,
    date = NULL,
    station = NULL,
    sample_volume_ml = coalesce(
      as.numeric(`samp vol [l]`) * 1000,
      as.numeric(`vol [m**3]`) * 1e6,
      as.numeric(`samp vol [m**3] (computed/converted)`) * 1e6,
      as.numeric(`samp vol [ml]`),
      as.numeric(`samp vol [ml] (computed/converted)`),
      as.numeric(sampled_volume) * 1000
    ),
    niskin_volume_ml = as.numeric(`niskin volume`) * 1000,
    depth_label = if_else(
      !is.na(depth_top_m) | !is.na(depth_bot_m),
      "range",
      NA_character_
    ),
    `samp vol [l]` = NULL,
    `vol [m**3]` = NULL,
    `samp vol [m**3] (computed/converted)` = NULL,
    `samp vol [ml]` = NULL,
    `samp vol [ml] (computed/converted)` = NULL,
    sampled_volume = NULL,
    `niskin volume` = NULL
  ) |>
  unite(
    "internal_id",
    `sample id`,
    `sample label`,
    event,
    `event 2`,
    campaign,
    station_tag,
    sep = " | ",
    na.rm = TRUE,
    remove = TRUE
  ) |>
  unite(
    "taxon",
    original_taxon,
    taxon_ex,
    na.rm = TRUE,
    remove = TRUE
  ) |>
  unite(
    "notes",
    notes,
    sample_comment,
    `sample comment`,
    sep = " | ",
    na.rm = TRUE,
    remove = TRUE
  ) |>
  mutate(
    notes = na_if(notes, "")
  )

count_data <- count_data |>
  filter(!is.na(depth_m))

# PROCESS ABUNDANCE AND PRESENCE/ABSENCE WRONG DATA ─────────────────────────

count_data <- count_data |>
  filter(!(abundance_l == 0 & !is.na(abundance_l) & presence_absence == TRUE))

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data <- count_data |>
  mutate(
    presence_absence = case_when(
      # a) If abundance > 0 → presence must be TRUE (covers FALSE and NA)
      !is.na(abundance_l) & abundance_l > 0 ~ TRUE,

      # b) If abundance == 0 and presence is NA → set to FALSE
      !is.na(abundance_l) & abundance_l == 0 & is.na(presence_absence) ~ FALSE,

      # c) All other cases (including abundance = NA with P/A TRUE or FALSE) stay as they are
      TRUE ~ presence_absence
    )
  )

count_data_db <- count_data |>
  mutate(
    datetime = if_else(
      is.na(datetime),
      NA_character_,
      paste0(datetime, "Z")
    ),
    depth = as.integer(round(depth_m)),
    s_count_method = coalesce(`method/device`, counting_method),
    s_volume_concentrated = as.numeric(niskin_volume_ml)
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
    depth_label = depth_label,
    depth_label_min = depth_top_m,
    depth_label_max = depth_bot_m,
    abundance = abundance_l,
    taxon,
    notes,
    s_count_method,
    s_sampling_method = sampling_method,
    s_preservation = preservative,
    s_volume_sampled = sample_volume_ml,
    s_volume_concentrated
  )

# PREPARE COUNT DATA ────────────────────────────────────────────────────────

count_data <- count_data_db

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
