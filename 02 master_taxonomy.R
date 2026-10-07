#       PHYTOPLANKTON TAXONOMY RESOLVER - UNIFIED WORKFLOW
# -----------------------------------------------------------------------------
# WORKFLOW:
# 0. Environment setup
# 1. Data extraction and dictionary loading
# 2. Taxon cleaning and epithet separation
# 3. Dictionary-based resolution (exact match)
# 4. WoRMS/GBIF automatic resolution (exact + taxamatch)
# 5. Semi-automatic resolution (fuzzy suggestions + manual confirmation)
# 6. Export results and update dictionary
# -----------------------------------------------------------------------------
# 20260513T12:00ADT
# Compatible with phytaxr
# Fernando Rayón Viña

# =========================================================================== #
# 0. ENVIRONMENT SETUP --------------------------------------------------------
# =========================================================================== #
rm(list = ls())

# Install/update phytaxr from GitHub
# remove.packages("phytaxr")
remotes::install_github(
  "RayonVina/phytaxr",
  ref = "main",
  force = TRUE,
  build = TRUE,
  upgrade = "never",
  quiet = FALSE
)
packageDescription("phytaxr")$RemoteSha

packages <- c(
  "tidyverse",
  "DBI",
  "RSQLite",
  "this.path",
  "furrr",
  "parallel",
  "phytaxr"
)

invisible(lapply(packages, require, character.only = TRUE))

# Parallelization
n_cores <- max(1, min(round(parallel::detectCores() * 0.9), 16))

# Paths
script_dir <- this.path::this.dir()
db_path <- Sys.getenv("PHYTOMAP_DB")

setwd(script_dir)

# =========================================================================== #
# 1. DATABASE CONNECTION AND DATA GATHERING -----------------------------------
# =========================================================================== #

db_path <- Sys.getenv("PHYTOMAP_DB")
con <- dbConnect(RSQLite::SQLite(), db_path)

# Extract unique taxon
unique_taxon <- tbl(con, "abundance_raw") |>
  select(taxon = t_taxon) |>
  filter(!is.na(taxon), taxon != "") |>
  distinct(taxon) |>
  collect() |>
  arrange(taxon)

# Load dictionary from database
tax_sinom <- tbl(con, "tax_sinom") |> collect()
tax_dict <- tbl(con, "tax_dict") |> collect()

# Create lookups
lookup_taxon_verbatim <- setNames(
  split(tax_sinom, seq_len(nrow(tax_sinom))),
  tax_sinom$taxon
)

lookup_taxon_clean <- tax_sinom |>
  filter(!is.na(taxon_clean)) |>
  distinct(taxon_clean, .keep_all = TRUE)

lookup_taxon_clean <- setNames(
  split(lookup_taxon_clean, seq_len(nrow(lookup_taxon_clean))),
  lookup_taxon_clean$taxon_clean
)

lookup_aphiaid <- setNames(
  split(tax_dict, seq_len(nrow(tax_dict))),
  as.character(tax_dict$aphiaid)
)

dbDisconnect(con)

# Memory purge
rm(tax_dict, tax_sinom, packages)

# ============================================================================ #
# 2. NAME CLEANING -------------------------------------------------------------
# ============================================================================ #
plan(multisession, workers = n_cores)
chunk_size <- ceiling(nrow(unique_taxon) / n_cores)
taxon_chunks <- split(
  unique_taxon,
  rep(1:n_cores, each = chunk_size, length.out = nrow(unique_taxon))
)

unique_taxon_clean <- future_map_dfr(
  taxon_chunks,
  \(chunk) {
    chunk |>
      phytaxr::normalize_characters() |>
      phytaxr::process_taxonomic_prefixes() |>
      phytaxr::process_incertae_entries() |>
      phytaxr::process_sp_entries() |>
      phytaxr::process_bracket_entries() |>
      phytaxr::move_size_to_epithet() |>
      phytaxr::process_epithet_entries() |>
      phytaxr::normalize_infraspecific_ranks() |>
      phytaxr::remove_dots() |>
      phytaxr::move_reproductive_structures() |>
      phytaxr::move_uncertainty_descriptors() |>
      phytaxr::move_morphological_descriptors() |>
      phytaxr::move_with_descriptors() |>
      phytaxr::move_formia_to_epithet() |>
      phytaxr::move_commas_to_epithet() |>
      phytaxr::move_authors_to_epithet() |>
      phytaxr::remove_sp_tokens() |>
      phytaxr::split_separator_entries() |>
      phytaxr::process_generic_taxa() |>
      phytaxr::clean_trailing_hyphens() |>
      phytaxr::remove_short_interstitial_tokens()
  },
  .options = furrr_options(seed = TRUE)
)

plan(sequential)

# Fallback: avoid empty taxon_clean after aggressive cleaning
unique_taxon_clean <- unique_taxon_clean |>
  mutate(
    taxon_clean = str_squish(taxon_clean),
    taxon_clean = na_if(taxon_clean, ""),
    needed_fallback = is.na(taxon_clean),
    taxon_clean = if_else(is.na(taxon_clean), taxon, taxon_clean),
    uncertain = if_else(needed_fallback, TRUE, uncertain)
  ) |>
  select(-needed_fallback)

# =========================================================================== #
# 3. DICTIONARY RESOLUTION (EXACT MATCH) --------------------------------------
# =========================================================================== #
## 3.0 OPTIONAL: EMPTY DICTIONARY (for further research) ----------------------
## Uncomment these lines to skip dictionary resolution and test pure API pipeline
# lookup_taxon_verbatim <- list()
# lookup_taxon_clean    <- list()
# lookup_aphiaid        <- list()

## 3.1 RESOLVE FROM DICTIONARY (PARALLELIZED) ---------------------------------
# Ensure internal structure (PENDINGD EXPORT)
unique_taxon_clean <- phytaxr::ensure_resolution_schema(unique_taxon_clean)

plan(multisession, workers = n_cores)

chunk_size <- ceiling(nrow(unique_taxon_clean) / n_cores)
taxon_chunks <- split(
  unique_taxon_clean,
  rep(1:n_cores, each = chunk_size, length.out = nrow(unique_taxon_clean))
)

unique_taxon_clean <- future_map_dfr(
  taxon_chunks,
  function(chunk) {
    chunk |>
      mutate(
        # Step 1: Search original taxon verbatim
        .temp = map(
          taxon,
          ~ phytaxr::lookup_taxon_info(
            .x,
            lookup_taxon_verbatim,
            include_notes = TRUE
          )
        )
      ) |>
      unnest_wider(.temp, names_sep = "_") |>
      mutate(
        matched_aphiaid = coalesce(.temp_matched_aphiaid, matched_aphiaid),
        matched_name = coalesce(.temp_matched_name, matched_name),
        accepted_name = coalesce(.temp_accepted_name, accepted_name),
        accepted_aphiaid = coalesce(.temp_accepted_aphiaid, accepted_aphiaid),
        taxonomic_status = coalesce(.temp_taxonomic_status, taxonomic_status),
        resolution_method = ifelse(
          .temp_found,
          .temp_resolution_method,
          resolution_method
        ),
        aphiaid = coalesce(.temp_aphiaid, aphiaid),
        marginalia = coalesce(.temp_marginalia, marginalia)
      ) |>
      select(-starts_with(".temp")) |>
      # Step 2: Search taxon_clean for unresolved
      mutate(
        .temp = if_else(
          is.na(aphiaid),
          map(
            taxon_clean,
            ~ phytaxr::lookup_taxon_info(
              .x,
              lookup_taxon_clean,
              include_notes = FALSE
            )
          ),
          map(
            seq_len(n()),
            ~ tibble(
              matched_aphiaid = NA_integer_,
              matched_name = NA_character_,
              accepted_name = NA_character_,
              accepted_aphiaid = NA_integer_,
              taxonomic_status = NA_character_,
              resolution_method = NA_character_,
              aphiaid = NA_integer_,
              marginalia = NA_character_,
              found = FALSE
            )
          )
        )
      ) |>
      unnest_wider(.temp, names_sep = "_") |>
      mutate(
        matched_aphiaid = coalesce(.temp_matched_aphiaid, matched_aphiaid),
        matched_name = coalesce(.temp_matched_name, matched_name),
        accepted_name = coalesce(.temp_accepted_name, accepted_name),
        accepted_aphiaid = coalesce(.temp_accepted_aphiaid, accepted_aphiaid),
        taxonomic_status = coalesce(.temp_taxonomic_status, taxonomic_status),
        resolution_method = ifelse(
          .temp_found,
          .temp_resolution_method,
          resolution_method
        ),
        aphiaid = coalesce(.temp_aphiaid, aphiaid)
      ) |>
      select(-starts_with(".temp")) |>
      # Step 3: Enrich with full taxonomy
      mutate(
        .temp = map(
          aphiaid,
          ~ phytaxr::lookup_taxonomy_info(.x, lookup_aphiaid)
        )
      ) |>
      unnest_wider(.temp, names_sep = "_") |>
      mutate(
        kingdom = coalesce(.temp_kingdom, kingdom),
        subkingdom = coalesce(.temp_subkingdom, subkingdom),
        infrakingdom = coalesce(.temp_infrakingdom, infrakingdom),
        phylum = coalesce(.temp_phylum, phylum),
        subphylum = coalesce(.temp_subphylum, subphylum),
        infraphylum = coalesce(.temp_infraphylum, infraphylum),
        parvphylum = coalesce(.temp_parvphylum, parvphylum),
        gigaclass = coalesce(.temp_gigaclass, gigaclass),
        superclass = coalesce(.temp_superclass, superclass),
        class = coalesce(.temp_class, class),
        subclass = coalesce(.temp_subclass, subclass),
        infraclass = coalesce(.temp_infraclass, infraclass),
        subterclass = coalesce(.temp_subterclass, subterclass),
        superorder = coalesce(.temp_superorder, superorder),
        order = coalesce(.temp_order, order),
        suborder = coalesce(.temp_suborder, suborder),
        infraorder = coalesce(.temp_infraorder, infraorder),
        parvorder = coalesce(.temp_parvorder, parvorder),
        superfamily = coalesce(.temp_superfamily, superfamily),
        family = coalesce(.temp_family, family),
        subfamily = coalesce(.temp_subfamily, subfamily),
        tribe = coalesce(.temp_tribe, tribe),
        genus = coalesce(.temp_genus, genus),
        subgenus = coalesce(.temp_subgenus, subgenus),
        section = coalesce(.temp_section, section),
        subsection = coalesce(.temp_subsection, subsection),
        species = coalesce(.temp_species, species),
        subspecies = coalesce(.temp_subspecies, subspecies),
        variety = coalesce(.temp_variety, variety),
        forma = coalesce(.temp_forma, forma),
        rank = coalesce(.temp_rank, rank),
        taxonomic_status = coalesce(.temp_taxonomic_status, taxonomic_status)
      ) |>
      select(-starts_with(".temp"))
  },
  .options = furrr_options(seed = TRUE)
)

plan(sequential)

# Memory purge
rm(
  lookup_taxon_verbatim,
  lookup_taxon_clean,
  lookup_aphiaid,
  taxon_chunks,
  chunk_size
)

# =========================================================================== #
# 4. EXTERNAL RESOLUTION (WoRMS / GBIF) ---------------------------------------
# =========================================================================== #

FUZZY_CONFIG <- list(
  max_levenshtein_distance = 3,
  min_similarity_score = 0.85,
  require_genus_match = TRUE
)

cat("\n=== STEP 4: EXTERNAL TAXONOMIC RESOLUTION ===\n")
cat(sprintf("Total entries to process: %d\n", nrow(unique_taxon_clean)))
cat(sprintf(
  "Unresolved at start: %d\n\n",
  sum(is.na(unique_taxon_clean$matched_aphiaid))
))

unique_taxon_clean <- phytaxr::search_worms_priority(
  unique_taxon_clean,
  col = "taxon_clean"
)
unique_taxon_clean <- phytaxr::search_worms_taxamatch(
  unique_taxon_clean,
  col = "taxon_clean"
)
unique_taxon_clean <- phytaxr::search_gbif_strict(
  unique_taxon_clean,
  col = "taxon_clean"
)
unique_taxon_clean <- phytaxr::resolve_taxonomic_status(
  unique_taxon_clean,
  col = "taxon_clean"
)
unique_taxon_clean <- phytaxr::search_worms_fuzzy_minor(
  unique_taxon_clean,
  col = "taxon_clean",
  fuzzy_config = FUZZY_CONFIG,
  ncores = n_cores
)
unique_taxon_clean <- phytaxr::get_taxonomy(
  unique_taxon_clean,
  col = "taxon_clean"
)

# Summary
cat("=== STEP 4 SUMMARY ===\n")
cat(sprintf("Total entries: %d\n", nrow(unique_taxon_clean)))
cat(sprintf(
  "Resolved: %d (%.1f%%)\n",
  sum(!is.na(unique_taxon_clean$matched_aphiaid)),
  sum(!is.na(unique_taxon_clean$matched_aphiaid)) /
    nrow(unique_taxon_clean) *
    100
))
cat("\nResolution methods:\n")
print(dplyr::count(unique_taxon_clean, resolution_method, sort = TRUE))

# Checkpoint save
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
save(unique_taxon_clean, file = paste0("taxonomy_step4_", timestamp, ".RData"))

rm(FUZZY_CONFIG, timestamp)

# =========================================================================== #
# 5. SEMI-AUTOMATIC FUZZY RESOLUTION ----------------------------------------
# =========================================================================== #
cat("\n=== STEP 5: SEMI-AUTOMATIC FUZZY RESOLUTION ===\n")
cat(sprintf("Total entries:       %d\n", nrow(unique_taxon_clean)))
cat(sprintf(
  "Unresolved at start: %d\n\n",
  sum(is.na(unique_taxon_clean$matched_aphiaid))
))

# Build vocabularies from already-resolved entries
# genus_vocab <- phytaxr::build_genus_vocabulary(unique_taxon_clean)
# epithet_vocab <- phytaxr::build_epithet_vocabulary(unique_taxon_clean)

# Build vocabularies (from checkpoint if available, else from data)
if (file.exists("phytaxr_step5_checkpoint.rds")) {
  ckpt_preview <- readRDS("phytaxr_step5_checkpoint.rds")
  genus_vocab <- ckpt_preview$genus_vocab %||%
    phytaxr::build_genus_vocabulary(unique_taxon_clean)
  epithet_vocab <- ckpt_preview$epithet_vocab %||%
    phytaxr::build_epithet_vocabulary(unique_taxon_clean)
  rm(ckpt_preview)
} else {
  genus_vocab <- phytaxr::build_genus_vocabulary(unique_taxon_clean)
  epithet_vocab <- phytaxr::build_epithet_vocabulary(unique_taxon_clean)
}

cat(sprintf("Genus vocabulary:    %d unique genera\n", length(genus_vocab)))
cat(sprintf(
  "Epithet vocabulary:  %d unique epithets\n\n",
  length(epithet_vocab)
))

# Run interactive fuzzy resolution loop.
# On first run, creates checkpoint_file after each accepted entry.
# On subsequent runs (new session), resumes automatically from the checkpoint.
unique_taxon_clean <- phytaxr::process_fuzzy_batch(
  df = unique_taxon_clean,
  genus_vocab = genus_vocab,
  epithet_vocab = epithet_vocab,
  batch_size = 10,
  checkpoint_file = "phytaxr_step5_checkpoint.rds",
  min_similarity = 0.85,
  max_suggestions = 15,
  max_longshot_suggestions = 10,
  longshot_threshold = 0.30,
  edit_max_dist = 3,
  edit_max_candidates = 15,
  timeout_sec = 15
)

# Coalesce aphiaids: prefer accepted (synonym-resolved) over matched
unique_taxon_clean <- unique_taxon_clean |>
  dplyr::mutate(aphiaid = dplyr::coalesce(accepted_aphiaid, matched_aphiaid))

# Checkpoint save
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
save(unique_taxon_clean, file = paste0("taxonomy_step5_", timestamp, ".RData"))

# =========================================================================== #
# 6. SORTING, ARRANGING AND EPITHET MANAGEMENT ------------------------------
# =========================================================================== #
## 6.1 SETUP ---------------------------------------------------------------
#Create target columns
unique_taxon_clean <- unique_taxon_clean |>
  mutate(
    .epi = tax_epithet,
    .res = resolution_notes
  )

#Combine sources function
combine_sources <- function(a, b, sep = " ") {
  case_when(
    !is.na(a) & !is.na(b) ~ paste(a, sep, b),
    !is.na(a) ~ a,
    !is.na(b) ~ b,
    TRUE ~ NA_character_
  )
}

## 6.2 PARSING NOTES EXTRACTION ---------------------------------------------
### 6.2.1 LTR FLAG ----------------------------------------------------------
unique_taxon_clean <- unique_taxon_clean |>
  mutate(
    ltr = str_detect(coalesce(.epi, ""), regex("\\bLTR\\b")) |
      str_detect(coalesce(.res, ""), regex("\\bLTR\\b"))
  )

### 6.2.2 AGGREGATES --------------------------------------------------------
extract_aggregates <- function(x) {
  ifelse(
    str_detect(coalesce(x, ""), regex("\\bagg\\.?\\b", ignore_case = TRUE)),
    str_squish(x),
    NA_character_
  )
}

unique_taxon_clean <- unique_taxon_clean |>
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

### 6.2.3 SIZE INFO ---------------------------------------------------------
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

unique_taxon_clean <- unique_taxon_clean |>
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

### 6.2.4 MORPHOTYPES -------------------------------------------------------
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

unique_taxon_clean <- unique_taxon_clean |>
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

### 6.2.5 LIFE-STAGES -------------------------------------------------------
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

unique_taxon_clean <- unique_taxon_clean |>
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

### 6.2.6 INFRASPECIFIC TAXA ------------------------------------------------
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

unique_taxon_clean <- unique_taxon_clean |>
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

### 6.2.7 INCERTAE ----------------------------------------------------------
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

unique_taxon_clean <- unique_taxon_clean |>
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

### 6.2.8 MARGINALIA --------------------------------------------------------
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

unique_taxon_clean <- unique_taxon_clean |>
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

## 6.3 CLEANUP --------------------------------------------------------------
unique_taxon_clean <- unique_taxon_clean |>
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

# 7. FINAL STEPS: CONSOLIDATION AND EXPORT ----------------------------------
## 7.1 CONSOLIDATE AGGREGATE DUPLICATES -------------------------------------
# Identify taxa split by split_separator_entries()
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

unique_taxon_clean <- consolidate_aggregates(unique_taxon_clean)
rm(consolidate_aggregates)

## 7.2 DATABASE EXPORTATION --------------------------------------------------
script_dir <- this.path::this.dir()

db_path <- Sys.getenv("PHYTOMAP_DB")
con <- dbConnect(RSQLite::SQLite(), db_path)

if (dbExistsTable(con, "master_taxonomy")) {
  dbRemoveTable(con, "master_taxonomy")
}

dbWriteTable(con, "master_taxonomy", unique_taxon_clean, row.names = FALSE)
cat("✓ master_taxonomy exported:", nrow(unique_taxon_clean), "rows\n")

dbDisconnect(con)
cat("✓ Connection closed.\n")
