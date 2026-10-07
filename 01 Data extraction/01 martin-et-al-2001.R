# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Martin-LeGresley-Strain-2001-report-2349
# Original author: Julianne Jager
# Fernando Rayón Viña
# Original version: 2022-08-04 Julianne Jager

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(readxl)
library(here)
library(lubridate)
library(janitor)
library(dplyr)
library(phytaxr)
library(DBI)
library(RSQLite)

# LOAD PATHS ────────────────────────────────────────────────────────────────

hereHelper <- function(dataset) {
  paths <- list.dirs(here::here(), recursive = TRUE, full.names = TRUE)
  paths <- paths[basename(paths) == dataset]
  if (length(paths) != 1L) {
    stop(sprintf("Expected one input directory named %s.", dataset))
  }
  paths
}

directory <- hereHelper("Martin-LeGresley-Strain-2001-report-2349")

fn <- here(directory, "original-data-Martin-et-al-2001.xlsx")

sheet_names <- excel_sheets(fn)

# IMPORT DATA ───────────────────────────────────────────────────────────────

d1a <- read_excel(fn, sheet_names[4])

d1b <- read_excel(fn, sheet_names[5])

d2 <- read_excel(fn, sheet_names[6])

station_metadata_a <- d1a |>
  mutate(
    date = mdy(date),
    day = day(date),
    month = month(date),
    year = year(date)
  ) |>
  select(!contains("_bottom") & !date) |>
  mutate(depth_m = 0)

station_metadata_b <- d1b |>
  mutate(
    date = mdy(date),
    day = day(date),
    month = month(date),
    year = year(date)
  ) |>
  mutate(depth_m = str_remove(depth_m, " m")) |>
  select(!date)

station_metadata <- station_metadata_a |>
  rbind(station_metadata_b)

# STATION METADATA ──────────────────────────────────────────────────────────

station_coords <- station_metadata |>
  select(station_id, latitude, longitude) |>
  distinct()

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data <- d2 |>
  pivot_longer(
    starts_with("station_"),
    names_to = "station_depth",
    values_to = "abundance_l"
  ) |>
  separate(col = station_depth, into = c("station_id", "depth_m"), sep = ":") |>
  mutate(
    station_id = str_remove(station_id, "station_"),
    depth_m = str_remove(depth_m, "depth_"),
    station_id = as.numeric(station_id)
  ) |>
  mutate(presence = ifelse(abundance_l == 0, FALSE, TRUE)) |>
  left_join(station_coords, by = "station_id") |>
  filter(
    !is.na(taxon),
    !is.na(abundance_l),
    !is.na(year),
    !is.na(month),
    !is.na(day)
  ) |>
  mutate(
    datetime = ymd(paste(year, month, day, sep = "-")),
    source = "martin-et-al-2001"
  ) |>
  select(
    year,
    month,
    day,
    datetime,
    depth = depth_m,
    latitude,
    longitude,
    t_taxon = taxon,
    abundance = abundance_l,
    source
  ) |>
  clean_names()

df <- count_data |>
  select(taxon = t_taxon) |>
  distinct(taxon)

df_clean <- df |>
  normalize_characters() |> # encoding, diacritics, invisible spaces
  process_taxonomic_prefixes() |> # O./C./F./P. rank prefixes
  process_incertae_entries() |> # cf., aff., incertae sedis, s.l., etc.
  process_sp_entries() |> # sp., spp., sp1, etc.
  process_bracket_entries() |> # bracket-delimited qualifiers []
  move_size_to_epithet() |> # size annotations (µm, mm, ranges)
  process_epithet_entries() |> # parenthesised epithets
  normalize_infraspecific_ranks() |> # var., subsp., f., ssp., cv.
  remove_dots() |> # stray dots (protects rank dots)
  move_reproductive_structures() |> # cysts, spores, filaments, …
  move_uncertainty_descriptors() |> # unidentified, unknown, unclassified, …
  move_morphological_descriptors() |> # centric, pennate, fusiform, …
  move_with_descriptors() |>
  move_formia_to_epithet() |> # forma designations
  move_commas_to_epithet() |> # post-comma strings
  move_authors_to_epithet() |> # authorships and dates
  remove_sp_tokens() |> # residual sp/spp/ssp in epithet
  split_separator_entries() |> # unfold /, &, +, or entries
  process_generic_taxa() |> # vernacular → scientific name
  clean_trailing_hyphens() |> # trailing/isolated hyphens
  remove_short_interstitial_tokens() |> # 1–2 char genus abbreviation artefacts
  apply_vernacular_corrections()

df_res <- run_resolution_pipeline(df_clean)

genus_vocab <- build_genus_vocabulary(df_res)

epithet_vocab <- build_epithet_vocabulary(df_res)

df_final <- process_fuzzy_batch(
  df = df_res,
  genus_vocab = genus_vocab,
  epithet_vocab = epithet_vocab,
  taxon_col = "taxon",
  taxon_clean_col = "taxon_clean",
  batch_size = 10,
  checkpoint_file = "phytaxr_step3_checkpoint.rds",
  min_similarity = 0.85,
  max_suggestions = 15,
  edit_max_dist = 3,
  timeout_sec = 15
)

# SETUP ─────────────────────────────────────────────────────────────────────

df_final <- df_final |>
  mutate(
    .epi = tax_epithet,
    .res = resolution_notes
  )

combine_sources <- function(a, b, sep = " ") {
  case_when(
    !is.na(a) & !is.na(b) ~ paste(a, sep, b),
    !is.na(a) ~ a,
    !is.na(b) ~ b,
    TRUE ~ NA_character_
  )
}

# LTR FLAG ──────────────────────────────────────────────────────────────────

df_final <- df_final |>
  mutate(
    ltr = str_detect(coalesce(.epi, ""), regex("\\bLTR\\b")) |
      str_detect(coalesce(.res, ""), regex("\\bLTR\\b"))
  )

# AGGREGATES ────────────────────────────────────────────────────────────────

extract_aggregates <- function(x) {
  ifelse(
    str_detect(coalesce(x, ""), regex("\\bagg\\.?\\b", ignore_case = TRUE)),
    str_squish(x),
    NA_character_
  )
}

df_final <- df_final |>
  mutate(
    .aggepi = extract_aggregates(.epi),
    .aggres = extract_aggregates(.res),
    aggregates = combine_sources(.aggepi, .aggres),
    .epi = if_else(
      !is.na(.aggepi),
      NA_character_,
      .epi
    ),
    .res = if_else(
      !is.na(.aggres),
      NA_character_,
      .res
    )
  ) |>
  select(-starts_with(".agg"))

# SIZE INFO ─────────────────────────────────────────────────────────────────

.SIZENUMERIC_PAT <- paste0(
  "(?:",
  "(?:length|width|radius|height|diam(?:eter)?)\\s*=\\s*",
  "\\d+(?:[.,]\\d+)?(?:\\s*(?:µm|um|μm|nm|mm))?",
  "|",
  "(?:[<>≤≥~]|<=|>=)?\\s*\\d+(?:[.,]\\d+)?",
  "(?:\\s*[-–]\\s*\\d+(?:[.,]\\d+)?)?",
  "\\s*(?:µm|um|μm|nm|mm)",
  "(?:\\s*(?:cell|width|length|diameter|height|diam))?",
  ")"
)

.SIZEQUALITATIVE_PAT <- "\\b(?:small|medium|large|tiny|huge|minor|major)\\b"

extract_sizeinfo <- function(epi, res) {
  num_epi <- str_extract_all(
    coalesce(epi, ""),
    regex(.SIZENUMERIC_PAT, ignore_case = TRUE)
  )
  num_res <- str_extract_all(
    coalesce(res, ""),
    regex(.SIZENUMERIC_PAT, ignore_case = TRUE)
  )
  qual <- str_extract_all(
    coalesce(epi, ""),
    regex(.SIZEQUALITATIVE_PAT, ignore_case = TRUE)
  )
  pmap_chr(
    list(num_epi, num_res, qual),
    ~ {
      vals <- unique(str_squish(c(..1, ..2, ..3)))
      vals <- vals[vals != ""]
      if (length(vals) == 0) NA_character_ else paste(vals, collapse = " | ")
    }
  )
}

df_final <- df_final |>
  mutate(
    size_info = extract_sizeinfo(.epi, .res),
    .epi = na_if(
      str_squish(
        str_remove_all(.epi, regex(.SIZENUMERIC_PAT, ignore_case = TRUE)) |>
          str_remove_all(regex(.SIZEQUALITATIVE_PAT, ignore_case = TRUE))
      ),
      ""
    ),
    .res = na_if(
      str_squish(
        str_remove_all(.res, regex(.SIZENUMERIC_PAT, ignore_case = TRUE))
      ),
      ""
    )
  )

# MORPHOTYPES ───────────────────────────────────────────────────────────────

.MORPHO_PAT <- paste0(
  "(?:^|\\b)(?:no\\s+|with(?:out)?\\s+)?",
  "(?:centric|pennate?|pennales|naked|thecate|thekat[ae]|athecate|armou?red|unarmou?red|",
  "gymnodinoid|peridinoid|gonyaulacoid|dinophysoid|prorocentroid|suessioid|",
  "holococcolit[ha]s?|heterococcolit[ha]s?|holococolit[ha]s?|heterococolit[ha]s?|",
  "round(?:ed)?|spherical|sphere|oval|ellipsoid(?:al)?|elliptic|fusiform|",
  "elongat(?:ed?|um)|cylindr(?:ical|ial|er)|conical|cone|prismatic|prism|",
  "stellate|star(?:-shaped)?|s-shaped|crescent|curved|spiral|helical|",
  "bacillus|rod-shaped|prolate|boat|drum|ribbon|fan|zigzag|",
  "rhombi[cd]|rhomboid|discoid|lanceolate|sigmoid|lunate|bilateral|radial|",
  "net(?:like)?|filament(?:s|ous)?|chain(?:s)?|cadenas?|",
  "coloni(?:al|y|es|as)|solitary|unicellular|single|coccoid|",
  "cluster(?:s)?|aggregate[d]?|agglutinated|clump(?:s)?|agglomerate[d]?|",
  "coiled|twisted|branched|arborescent|",
  "flagellate[ds]?|flagellated|non-flagellated|ciliate[ds]?|",
  "loricae?|empty|full|fusing|diatom|",
  "benthic|epiphytic|planktonic|pelagic|neritic|oceanic|",
  "heterotrophic|autotrophic|mixotrophic|phototrophic|",
  "olive-green|blue-?green|olive|",
  "naviculoid|cylindrotheca|double|truncated)(?=$|\\b|[,;\\s])"
)

extract_morphotypes <- function(epi, res) {
  from_epi <- str_extract_all(
    coalesce(epi, ""),
    regex(.MORPHO_PAT, ignore_case = TRUE)
  )
  from_res <- str_extract_all(
    coalesce(res, ""),
    regex(.MORPHO_PAT, ignore_case = TRUE)
  )
  map2_chr(
    from_epi,
    from_res,
    ~ {
      vals <- unique(tolower(str_squish(c(.x, .y))))
      vals <- vals[vals != ""]
      if (length(vals) == 0) NA_character_ else paste(vals, collapse = " | ")
    }
  )
}

df_final <- df_final |>
  mutate(
    morphotypes = extract_morphotypes(.epi, .res),
    .epi = na_if(
      str_squish(
        str_remove_all(.epi, regex(.MORPHO_PAT, ignore_case = TRUE))
      ),
      ""
    ),
    .res = na_if(
      str_squish(
        str_remove_all(.res, regex(.MORPHO_PAT, ignore_case = TRUE))
      ),
      ""
    )
  )

# LIFE-STAGES ───────────────────────────────────────────────────────────────

.LIFESTAGE_PAT <- paste0(
  "(?:^|\\b)",
  "(?:resting\\s+(?:spores?|cysts?)|",
  "spores?|cysts?|eggs?|ova|gametes?|",
  "planozygotes?|zygotes?|auxospores?|akinetes?|",
  "duplets?|triplets?|quadruplets?|",
  "larvae?|juveniles?|adults?|",
  "copepodites?|naupli(?:us|i)|",
  "empty\\s+loricae?|loricae?|",
  "thecae?|",
  "vegetative\\s+cells?|cells?\\s+stage|",
  "life\\s+stage|",
  "statospores?|cyst\\s+stage)",
  "(?=$|\\b|[,;\\s])"
)

extract_lifestages <- function(epi, res) {
  from_epi <- str_extract_all(
    coalesce(epi, ""),
    regex(.LIFESTAGE_PAT, ignore_case = TRUE)
  )
  from_res <- str_extract_all(
    coalesce(res, ""),
    regex(.LIFESTAGE_PAT, ignore_case = TRUE)
  )
  map2_chr(
    from_epi,
    from_res,
    ~ {
      vals <- unique(tolower(str_squish(c(.x, .y))))
      vals <- vals[vals != ""]
      if (length(vals) == 0) NA_character_ else paste(vals, collapse = " | ")
    }
  )
}

df_final <- df_final |>
  mutate(
    life_stages = extract_lifestages(.epi, .res),
    .epi = na_if(
      str_squish(
        str_remove_all(.epi, regex(.LIFESTAGE_PAT, ignore_case = TRUE))
      ),
      ""
    ),
    .res = na_if(
      str_squish(
        str_remove_all(.res, regex(.LIFESTAGE_PAT, ignore_case = TRUE))
      ),
      ""
    )
  )

# INFRASPECIFIC TAXA ────────────────────────────────────────────────────────

.INFRASP_PAT <- paste0(
  "(?:",
  "(?:[a-z]+-?[a-z]+\\s+)?(?:var\\.?|v\\.(?!\\s*[A-Z]))\\s*[a-z][a-z-]+|",
  "(?:f\\.?|forma)\\s+[a-z][a-z-]+|",
  "(?:subsp\\.?|ssp\\.?|subspecies)\\s+[a-z][a-z-]+|",
  "infraspecific\\s+epithet:\\s*['\"]?[a-z][a-z.\\s-]+['\"]?|",
  "(?:cv\\.?|cultivar)\\s+[a-z][a-z-]+",
  ")"
)

extract_infraspecific <- function(epi, res) {
  from_epi <- str_extract_all(
    coalesce(epi, ""),
    regex(.INFRASP_PAT, ignore_case = FALSE)
  )
  from_res <- str_extract_all(
    coalesce(res, ""),
    regex(.INFRASP_PAT, ignore_case = FALSE)
  )
  map2_chr(
    from_epi,
    from_res,
    ~ {
      vals <- unique(str_squish(c(.x, .y)))
      vals <- vals[vals != ""]
      if (length(vals) == 0) NA_character_ else paste(vals, collapse = " | ")
    }
  )
}

df_final <- df_final |>
  mutate(
    infraspecific_taxa = extract_infraspecific(.epi, .res),
    .epi = na_if(
      str_squish(
        str_remove_all(.epi, regex(.INFRASP_PAT, ignore_case = FALSE))
      ),
      ""
    ),
    .res = na_if(
      str_squish(
        str_remove_all(.res, regex(.INFRASP_PAT, ignore_case = FALSE))
      ),
      ""
    )
  )

# INCERTAE ──────────────────────────────────────────────────────────────────

.INCERTAEPAT <- paste0(
  "(?i)",
  "\\bcf\\.?\\b|\\baff\\.?\\b|\\bnr\\.?\\b|",
  "\\bincertae\\b|incertae\\s*sedis|\\bsedis\\b|",
  "\\bindet\\.?\\b|\\bindeterminat[ae]+\\b|",
  "\\bundetermined\\b|\\bundifferentiated\\b|",
  "\\bunidentified\\b|\\bunid\\.?\\b|\\bundet\\.?\\b|",
  "\\bunclassified\\b|\\bnot\\s+classified\\b|",
  "\\buncertain\\b|\\bunknown\\b|\\bunresolved\\b|",
  "\\?"
)

df_final <- df_final |>
  dplyr::mutate(
    uncertain = uncertain |
      stringr::str_detect(
        coalesce(.epi, ""),
        regex(.INCERTAEPAT, ignore_case = TRUE)
      ) |
      stringr::str_detect(
        coalesce(.res, ""),
        regex(.INCERTAEPAT, ignore_case = TRUE)
      ) |
      stringr::str_detect(
        coalesce(taxon, ""),
        regex(.INCERTAEPAT, ignore_case = TRUE)
      )
  )

# MARGINALIA ────────────────────────────────────────────────────────────────

.MARGINALIA_PAT <- paste0(
  "(?:^|\\b)(?:",

  # Identification uncertainty qualifiers
  "unidentified|unid\\.?|undet\\.?|undetermined|",
  "unclassified|not\\s+classified|",
  "indeterminat[ae]|indet\\.?|undifferentiated|",
  "incertae(?:\\s+sedis)?|",
  "order\\s+uncertain|uncertain|unknown|unresolved|",

  # Taxonomic grouping modifiers
  "group|complex|aggregate[d]?|sensu\\s+(?:lato|stricto)|s\\.l\\.|s\\.s\\.|",
  "clade|lineage|radiation|assemblage|",

  # Data quality / observational notes
  "fractionated|mixed|other|",
  "misc(?:ellaneous)?|various|varia|",
  "undescribed|informal|provisional|",

  # Morphotype/strain codes (letter or letter+number, clearly non-taxonomic)
  "(?:morph|strain|type|form(?!a)|variant|clone)\\s*[A-Za-z]?\\d*|",

  # Contextual modifiers
  "like(?:\\s+incertae)?|",
  "sensu\\s+[A-Z][a-z]+(?:\\s+et\\s+al\\.?)?(?:\\s+\\d{4})?",

  ")(?=$|\\b|[,;\\s\"'])"
)

extract_marginalia <- function(epi, res) {
  from_epi <- str_extract_all(
    coalesce(epi, ""),
    regex(.MARGINALIA_PAT, ignore_case = TRUE)
  )
  from_res <- str_extract_all(
    coalesce(res, ""),
    regex(.MARGINALIA_PAT, ignore_case = TRUE)
  )
  map2_chr(
    from_epi,
    from_res,
    ~ {
      vals <- unique(tolower(str_squish(c(.x, .y))))
      vals <- vals[vals != ""]
      if (length(vals) == 0) NA_character_ else paste(vals, collapse = " | ")
    }
  )
}

df_final <- df_final |>
  mutate(
    marginalia = extract_marginalia(.epi, .res),
    .epi = na_if(
      str_squish(
        str_remove_all(.epi, regex(.MARGINALIA_PAT, ignore_case = TRUE))
      ),
      ""
    ),
    .res = na_if(
      str_squish(
        str_remove_all(.res, regex(.MARGINALIA_PAT, ignore_case = TRUE))
      ),
      ""
    )
  )

# CLEANUP ───────────────────────────────────────────────────────────────────

df_final <- df_final |>
  mutate(
    notes = pmap_chr(
      list(.epi, .res, marginalia),
      ~ {
        parts <- discard(c(..1, ..2, ..3), is.na)
        parts <- parts[str_squish(parts) != ""]
        if (length(parts) == 0) {
          NA_character_
        } else {
          paste(unique(parts), collapse = " | ")
        }
      }
    )
  ) |>
  select(-c(.epi, .res, marginalia, tax_epithet, resolution_notes))

# CONSOLIDATE AGGREGATE DUPLICATES ──────────────────────────────────────────

consolidate_aggregates <- function(df) {
  agg_taxa <- df |>
    group_by(taxon) |>
    filter(n() > 1) |>
    ungroup() |>
    pull(taxon) |>
    unique()

  if (length(agg_taxa) == 0) {
    message("No aggregate duplicates found.")
    return(df)
  }

  partner_notes <- df |>
    filter(taxon %in% agg_taxa) |>
    group_by(taxon) |>
    slice(-1) |>
    summarise(
      .partner_note = paste(
        sprintf(
          "agg. partner: %s | matched: %s | aphiaID: %s | accepted: %s",
          taxon_clean,
          matched_name,
          aphiaid,
          accepted_name
        ),
        collapse = " | "
      ),
      .groups = "drop"
    )

  df |>
    group_by(taxon) |>
    slice(1) |>
    ungroup() |>
    left_join(partner_notes, by = "taxon") |>
    mutate(
      notes = case_when(
        !is.na(.partner_note) & !is.na(notes) ~
          paste(notes, .partner_note, sep = " | "),
        !is.na(.partner_note) ~ .partner_note,
        TRUE ~ notes
      )
    ) |>
    select(-.partner_note)
}

df_final <- consolidate_aggregates(df_final)

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data_full <- count_data |>
  left_join(df_final, by = c("t_taxon" = "taxon"))

abundance_raw_insert <- count_data_full |>
  mutate(
    datetime = as.character(datetime),
    depth = as.numeric(depth),
    n_uncertain = as.integer(uncertain),
    n_ltr = as.integer(ltr),
    n_flag_for_removal = as.integer(flag_for_removal),
    depth_label = "exact"
  ) |>
  select(
    source,
    datetime,
    year,
    month,
    day,
    latitude,
    longitude,
    depth,
    depth_label,
    abundance,
    t_taxon,
    t_taxon_clean = taxon_clean,
    n_uncertain,
    t_matched_name = matched_name,
    t_matched_aphiaid = matched_aphiaid,
    t_accepted_name = accepted_name,
    t_accepted_aphiaid = accepted_aphiaid,
    t_taxonomic_status = taxonomic_status,
    t_resolution_method = resolution_method,
    t_aphiaid = aphiaid,
    t_kingdom = kingdom,
    t_subkingdom = subkingdom,
    t_infrakingdom = infrakingdom,
    t_phylum = phylum,
    t_subphylum = subphylum,
    t_infraphylum = infraphylum,
    t_parvphylum = parvphylum,
    t_gigaclass = gigaclass,
    t_superclass = superclass,
    t_class = class,
    t_subclass = subclass,
    t_infraclass = infraclass,
    t_subterclass = subterclass,
    t_superorder = superorder,
    t_order = order,
    t_suborder = suborder,
    t_infraorder = infraorder,
    t_parvorder = parvorder,
    t_superfamily = superfamily,
    t_family = family,
    t_subfamily = subfamily,
    t_tribe = tribe,
    t_genus = genus,
    t_subgenus = subgenus,
    t_section = section,
    t_subsection = subsection,
    t_species = species,
    t_subspecies = subspecies,
    t_variety = variety,
    t_forma = forma,
    t_rank = rank,
    n_ltr,
    n_aggregates = aggregates,
    n_size_info = size_info,
    n_morphotypes = morphotypes,
    n_life_stages = life_stages,
    n_infraspecific_taxa = infraspecific_taxa,
    n_flag_for_removal,
    n_notes = notes
  )

# PREPARE COUNT DATA ────────────────────────────────────────────────────────

count_data <- abundance_raw_insert

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
