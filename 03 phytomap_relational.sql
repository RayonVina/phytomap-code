--  PHYTOMAP RELATIONAL STRUCTURE CONSTRUCTION
--------------------------------------------------------------------------------
-- WORKFLOW:
-- 1. Tables and relations creation
-- 2. Aditional author-curated inputs
-- 3. Filter and curation of data
-- 4. Deterministic relational IDs
-- 5. Load dimensions and dependent tables
-- 6. Denormalized view creation
-- 7. Validate counts and reconstructed values
-- 8. Remove useless tables
-- 9. Verify the release object and database integrity
--------------------------------------------------------------------------------
-- 20261006T16:00ADT
-- Compatible with phytaxr
-- Fernando Rayón Viña

.bail on
PRAGMA foreign_keys = ON;
PRAGMA recursive_triggers = OFF;
BEGIN IMMEDIATE;

---------------------------------------------------------------
-- 1. Tables, constraints, indexes and modification triggers --
---------------------------------------------------------------
CREATE TABLE 
  coordinates (
    id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, 
    latitude REAL NOT NULL CHECK (
      latitude >= -90.0 
      AND latitude <= 90.0
    ), 
    longitude REAL NOT NULL CHECK (
      longitude >= -180.0 
      AND longitude <= 180.0
    ), 
    coord_precision TEXT, 
    notes TEXT, 
    created TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP, 
    modified TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
  ); 
CREATE TABLE 
  event_datetime (
    id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, 
    year INTEGER NOT NULL CHECK (year >= 1), 
    month INTEGER NOT NULL CHECK (
      month BETWEEN 1 
      AND 12
    ), 
    day INTEGER NOT NULL, 
    time_value TEXT, 
    date_precision TEXT, 
    day_is_imputed INTEGER NOT NULL DEFAULT 0 CHECK (
      day_is_imputed IN (0, 1)
    ), 
    notes TEXT, 
    created TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP, 
    modified TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP, 
    CONSTRAINT ck_event_datetime_valid_calendar_day CHECK (
      day >= 1 
      AND day <= CASE 
        month WHEN 2 THEN CASE 
          WHEN (
            year % 4 = 0 
            AND year % 100 != 0
          ) 
          OR (year % 400 = 0) THEN 29 ELSE 28
        END WHEN 4  THEN 30 
            WHEN 6  THEN 30 
            WHEN 9  THEN 30 
            WHEN 11 THEN 30 
            ELSE 31
      END
    ), 
    CONSTRAINT ck_event_datetime_time_value CHECK (
      time_value IS NULL 
      OR (
        length(time_value) = 5 
        AND time_value GLOB '[0-2][0-9]:[0-5][0-9]' 
        AND CAST(
          substr(time_value, 1, 2) AS INTEGER
        ) BETWEEN 0 
        AND 23
      ) 
      OR (
        length(time_value) = 8 
        AND time_value GLOB '[0-2][0-9]:[0-5][0-9]:[0-5][0-9]' 
        AND CAST(
          substr(time_value, 1, 2) AS INTEGER
        ) BETWEEN 0 
        AND 23
      ) 
      OR (
        length(time_value) = 9 
        AND time_value GLOB '[0-2][0-9]:[0-5][0-9]:[0-5][0-9]Z' 
        AND CAST(
          substr(time_value, 1, 2) AS INTEGER
        ) BETWEEN 0 
        AND 23
      )
    )
  ); 
CREATE TABLE 
  depth (
    id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, 
    depth REAL NOT NULL CHECK (depth >= 0.0), 
    depth_label TEXT, 
    depth_label_min REAL CHECK (
      depth_label_min IS NULL 
      OR depth_label_min >= 0.0
    ), 
    depth_label_max REAL CHECK (
      depth_label_max IS NULL 
      OR depth_label_max >= 0.0
    ), 
    notes TEXT, 
    created TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP, 
    modified TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP, 
    CONSTRAINT ck_depth_label_bounds CHECK (
      depth_label_min IS NULL 
      OR depth_label_max IS NULL 
      OR depth_label_min <= depth_label_max
    )
  ); 
CREATE TABLE 
  taxonomy (
    id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, 
    aphiaid INTEGER NOT NULL UNIQUE, kingdom TEXT, 
    subkingdom TEXT, infrakingdom TEXT, 
    phylum TEXT, subphylum TEXT, infraphylum TEXT, 
    parvphylum TEXT, gigaclass TEXT, superclass TEXT, 
    class TEXT, subclass TEXT, infraclass TEXT, 
    subterclass TEXT, superorder TEXT, 
    taxon_order TEXT, suborder TEXT, infraorder TEXT, 
    parvorder TEXT, superfamily TEXT, 
    family TEXT, subfamily TEXT, tribe TEXT, 
    genus TEXT, subgenus TEXT, section TEXT, 
    subsection TEXT, species TEXT, subspecies TEXT, 
    variety TEXT, forma TEXT, rank TEXT, 
    taxonomic_status TEXT, notes TEXT, 
    created TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP, 
    modified TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
  ); 
CREATE TABLE 
  sampling_protocol (
    id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, 
    count_method TEXT, enumeration_method TEXT, 
    sampling_method TEXT, sampling_method_notes TEXT, 
    concentration_method TEXT, concentration_method_notes TEXT, 
    volume_sampled TEXT, microscopy_method TEXT, 
    microscopy_notes TEXT, preservation TEXT, 
    volume_concentrated TEXT, volume_counted TEXT, 
    magnification TEXT, counting_chamber TEXT, 
    notes TEXT, created TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP, 
    modified TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
  ); 
CREATE TABLE 
  source (
    id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, 
    source TEXT NOT NULL, subsource TEXT, 
    metadata TEXT, url TEXT, notes TEXT, 
    created TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP, 
    modified TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
  ); 
CREATE TABLE 
  event (
    id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, 
    event_datetime_id INTEGER NOT NULL, 
    coordinates_id INTEGER NOT NULL, 
    source_id INTEGER NOT NULL, 
    notes TEXT, 
    created TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP, 
    modified TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP, 
    CONSTRAINT uq_event_datetime_coordinates_source UNIQUE (
      event_datetime_id, coordinates_id, 
      source_id
    ), 
    CONSTRAINT fk_event_event_datetime FOREIGN KEY (event_datetime_id) REFERENCES event_datetime(id) ON UPDATE CASCADE ON DELETE RESTRICT, 
    CONSTRAINT fk_event_coordinates FOREIGN KEY (coordinates_id) REFERENCES coordinates(id) ON UPDATE CASCADE ON DELETE RESTRICT, 
    CONSTRAINT fk_event_source FOREIGN KEY (source_id) REFERENCES source(id) ON UPDATE CASCADE ON DELETE RESTRICT
  ); 
CREATE TABLE 
  sample (
    id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, 
    event_id INTEGER NOT NULL, 
    depth_id INTEGER NOT NULL, 
    notes TEXT, 
    created TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP, 
    modified TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP, 
    CONSTRAINT uq_sample_event_depth UNIQUE (event_id, depth_id), 
    CONSTRAINT fk_sample_event FOREIGN KEY (event_id) REFERENCES event(id) ON UPDATE CASCADE ON DELETE RESTRICT, 
    CONSTRAINT fk_sample_depth FOREIGN KEY (depth_id) REFERENCES depth(id) ON UPDATE CASCADE ON DELETE RESTRICT
  ); 
CREATE TABLE 
  taxon_record (
    id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, 
    taxonomy_id INTEGER NOT NULL, 
    taxon_raw TEXT, 
    taxon_clean TEXT, 
    taxon_matched TEXT, 
    aphiaid_matched INTEGER, 
    taxon_accepted TEXT, 
    aphiaid_accepted INTEGER, 
    resolution_method TEXT, 
    uncertain INTEGER CHECK (
      uncertain IS NULL 
      OR uncertain IN (0, 1)
    ), 
    ltr INTEGER CHECK (
      ltr IS NULL 
      OR ltr IN (0, 1)
    ), 
    aggregates TEXT, 
    size_info TEXT, 
    morphotype TEXT, 
    life_stage TEXT, 
    infraspecific_taxa TEXT, 
    trophic TEXT, 
    local_authority TEXT, 
    notes TEXT, 
    created TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP, 
    modified TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP, 
    CONSTRAINT fk_taxon_record_taxonomy FOREIGN KEY (taxonomy_id) REFERENCES taxonomy(id) ON UPDATE CASCADE ON DELETE RESTRICT
  ); 
CREATE TABLE 
  abundance (
    id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, 
    sample_id INTEGER NOT NULL, 
    taxon_record_id INTEGER NOT NULL, 
    sampling_protocol_id INTEGER NOT NULL, 
    abundance REAL NOT NULL CHECK (abundance >= 0.0), 
    notes TEXT, 
    created TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP, 
    modified TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP, 
    CONSTRAINT uq_abundance_sample_taxon_protocol_value UNIQUE (
      sample_id, taxon_record_id, sampling_protocol_id, 
      abundance
    ), 
    CONSTRAINT fk_abundance_sample FOREIGN KEY (sample_id) REFERENCES sample(id) ON UPDATE CASCADE ON DELETE RESTRICT, 
    CONSTRAINT fk_abundance_taxon_record FOREIGN KEY (taxon_record_id) REFERENCES taxon_record(id) ON UPDATE CASCADE ON DELETE RESTRICT, 
    CONSTRAINT fk_abundance_sampling_protocol FOREIGN KEY (sampling_protocol_id) REFERENCES sampling_protocol(id) ON UPDATE CASCADE ON DELETE RESTRICT
  ); 
CREATE TABLE 
  IF NOT EXISTS metadata_dictionary (
    sort_order INTEGER, 
    table_name TEXT NOT NULL, 
    variable TEXT NOT NULL, 
    description TEXT, 
    data_type TEXT, 
    possible_values TEXT, 
    missing_value_interpretation TEXT, 
    notes TEXT, 
    CONSTRAINT METADATA_DICTIONARY_PK PRIMARY KEY (table_name, variable)
  ); 
CREATE INDEX 
  ix_abundance_sampling_protocol_id ON abundance(sampling_protocol_id); 
CREATE INDEX 
  ix_abundance_taxon_record_id ON abundance(taxon_record_id); 
CREATE INDEX 
  ix_event_coordinates_id ON event(coordinates_id); 
CREATE INDEX 
  ix_event_event_datetime_id ON event(event_datetime_id); 
CREATE INDEX 
  ix_event_source_id ON event(source_id); 
CREATE INDEX 
  ix_sample_depth_id ON sample(depth_id); 
CREATE INDEX 
  ix_taxon_record_aphiaid_accepted ON taxon_record(aphiaid_accepted); 
CREATE INDEX 
  ix_taxon_record_aphiaid_matched ON taxon_record(aphiaid_matched); 
CREATE INDEX 
  ix_taxon_record_taxonomy_id ON taxon_record(taxonomy_id); 
CREATE INDEX 
  ix_taxonomy_aphiaid ON taxonomy(aphiaid); CREATE UNIQUE INDEX uq_coordinates_natural_key ON coordinates (
  latitude, 
  longitude, 
  COALESCE(coord_precision, '')
); CREATE UNIQUE INDEX uq_depth_natural_key ON depth (
  depth, 
  COALESCE(depth_label, ''), 
  COALESCE(depth_label_min, -1.0), 
  COALESCE(depth_label_max, -1.0)
); CREATE UNIQUE INDEX uq_event_datetime_natural_key ON event_datetime (
  year, 
  month, 
  day, 
  COALESCE(time_value, ''), 
  COALESCE(date_precision, ''), 
  day_is_imputed
); CREATE UNIQUE INDEX uq_sampling_protocol_natural_key ON sampling_protocol (
  COALESCE(count_method, ''), 
  COALESCE(enumeration_method, ''), 
  COALESCE(sampling_method, ''), 
  COALESCE(sampling_method_notes, ''), 
  COALESCE(concentration_method, ''), 
  COALESCE(concentration_method_notes, ''), 
  COALESCE(volume_sampled, ''), 
  COALESCE(microscopy_method, ''), 
  COALESCE(microscopy_notes, ''), 
  COALESCE(preservation, ''), 
  COALESCE(volume_concentrated, ''), 
  COALESCE(volume_counted, ''), 
  COALESCE(magnification, ''), 
  COALESCE(counting_chamber, '')
); CREATE UNIQUE INDEX uq_source_natural_key ON source (
  source, 
  COALESCE(subsource, '')
); CREATE UNIQUE INDEX uq_taxon_record_natural_key ON taxon_record (
  taxonomy_id, 
  COALESCE(taxon_raw, ''), 
  COALESCE(taxon_clean, ''), 
  COALESCE(taxon_matched, ''), 
  COALESCE(aphiaid_matched, -1), 
  COALESCE(taxon_accepted, ''), 
  COALESCE(aphiaid_accepted, -1), 
  COALESCE(resolution_method, ''), 
  COALESCE(uncertain, -1), 
  COALESCE(ltr, -1), 
  COALESCE(aggregates, ''), 
  COALESCE(size_info, ''), 
  COALESCE(morphotype, ''), 
  COALESCE(life_stage, ''), 
  COALESCE(infraspecific_taxa, ''), 
  COALESCE(trophic, ''), 
  COALESCE(local_authority, '')
); 
CREATE TRIGGER 
  trg_abundance_modified 
AFTER 
UPDATE 
  OF sample_id, 
  taxon_record_id, 
  sampling_protocol_id, 
  abundance, 
  notes ON abundance FOR EACH ROW BEGIN 
    UPDATE 
      abundance 
    SET 
      modified = strftime('%Y-%m-%dT%H:%M:%fZ', 'now') 
    WHERE 
      id = NEW.id;
END; 
CREATE TRIGGER 
  trg_coordinates_modified 
AFTER 
UPDATE 
  OF latitude, 
  longitude, 
  coord_precision, 
  notes ON coordinates FOR EACH ROW BEGIN 
    UPDATE 
      coordinates 
    SET 
      modified = strftime('%Y-%m-%dT%H:%M:%fZ', 'now') 
    WHERE 
      id = NEW.id;
END; 
CREATE TRIGGER 
  trg_depth_modified 
AFTER 
UPDATE 
  OF depth, 
  depth_label, 
  depth_label_min, 
  depth_label_max, 
  notes ON depth FOR EACH ROW BEGIN 
    UPDATE 
      depth 
    SET 
      modified = strftime('%Y-%m-%dT%H:%M:%fZ', 'now') 
    WHERE 
      id = NEW.id;
END; 
CREATE TRIGGER 
  trg_event_datetime_modified 
AFTER 
UPDATE 
  OF year, 
  month, 
  day, 
  time_value, 
  date_precision, 
  day_is_imputed, 
  notes ON event_datetime FOR EACH ROW BEGIN 
    UPDATE 
      event_datetime 
    SET 
      modified = strftime('%Y-%m-%dT%H:%M:%fZ', 'now') 
    WHERE 
      id = NEW.id;
END; 
CREATE TRIGGER 
  trg_event_modified 
AFTER 
UPDATE 
  OF event_datetime_id, 
  coordinates_id, 
  source_id, 
  notes ON event FOR EACH ROW BEGIN 
    UPDATE 
      event 
    SET 
      modified = strftime('%Y-%m-%dT%H:%M:%fZ', 'now') 
    WHERE 
      id = NEW.id;
END; 
CREATE TRIGGER 
  trg_sample_modified 
AFTER 
UPDATE 
  OF event_id, 
  depth_id, 
  notes ON sample FOR EACH ROW BEGIN 
    UPDATE 
      sample 
    SET 
      modified = strftime('%Y-%m-%dT%H:%M:%fZ', 'now') 
    WHERE 
      id = NEW.id;
END; 
CREATE TRIGGER 
  trg_sampling_protocol_modified 
AFTER 
UPDATE 
  OF count_method, 
  enumeration_method, 
  sampling_method, 
  sampling_method_notes, 
  concentration_method, 
  concentration_method_notes, 
  volume_sampled, 
  microscopy_method, 
  microscopy_notes, 
  preservation, 
  volume_concentrated, 
  volume_counted, 
  magnification, 
  counting_chamber, 
  notes ON sampling_protocol FOR EACH ROW BEGIN 
    UPDATE 
      sampling_protocol 
    SET 
      modified = strftime('%Y-%m-%dT%H:%M:%fZ', 'now') 
    WHERE 
      id = NEW.id;
END; 
CREATE TRIGGER 
  trg_source_modified 
AFTER 
UPDATE 
  OF source, 
  subsource, 
  metadata, 
  url, 
  notes ON source FOR EACH ROW BEGIN 
    UPDATE 
      source 
    SET 
      modified = strftime('%Y-%m-%dT%H:%M:%fZ', 'now') 
    WHERE 
      id = NEW.id;
END; 
CREATE TRIGGER 
  trg_taxon_record_modified 
AFTER 
UPDATE 
  OF taxonomy_id, 
  taxon_raw, 
  taxon_clean, 
  taxon_matched, 
  aphiaid_matched, 
  taxon_accepted, 
  aphiaid_accepted, 
  resolution_method, 
  uncertain, 
  ltr, 
  aggregates, 
  size_info, 
  morphotype, 
  life_stage, 
  infraspecific_taxa, 
  trophic, 
  local_authority, 
  notes ON taxon_record FOR EACH ROW BEGIN 
    UPDATE 
      taxon_record 
    SET 
      modified = strftime('%Y-%m-%dT%H:%M:%fZ', 'now') 
    WHERE 
      id = NEW.id;
END; 
CREATE TRIGGER 
  trg_taxonomy_modified 
AFTER 
UPDATE 
  OF aphiaid, 
  kingdom, 
  subkingdom, 
  infrakingdom, 
  phylum, 
  subphylum, 
  infraphylum, 
  parvphylum, 
  gigaclass, 
  superclass, 
  class, 
  subclass, 
  infraclass, 
  subterclass, 
  superorder, 
  taxon_order, 
  suborder, 
  infraorder, 
  parvorder, 
  superfamily, 
  family, 
  subfamily, 
  tribe, 
  genus, 
  subgenus, 
  section, 
  subsection, 
  species, 
  subspecies, 
  variety, 
  forma, 
  rank, 
  taxonomic_status, 
  notes ON taxonomy FOR EACH ROW BEGIN 
    UPDATE 
      taxonomy 
    SET 
      modified = strftime('%Y-%m-%dT%H:%M:%fZ', 'now') 
    WHERE 
      id = NEW.id;
END; 
---------------------------------------
-- 2. Optional author-curated inputs --
---------------------------------------
-- CSVs are UTF-8 with a header in the exact column order below.
-- Blank nullable CSV fields are interpreted as SQL NULL.
-- A sampling row replaces the complete protocol for that original observation.
-- An uncertainty row overrides only the corresponding uncertainty flag.
CREATE TEMP TABLE curated_source (
  source TEXT NOT NULL, subsource TEXT, 
  metadata TEXT, url TEXT, notes TEXT
); CREATE TEMP TABLE curated_sampling (
  original_id INTEGER PRIMARY KEY NOT NULL, 
  count_method TEXT, enumeration_method TEXT, 
  sampling_method TEXT, sampling_method_notes TEXT, 
  concentration_method TEXT, concentration_method_notes TEXT, 
  volume_sampled TEXT, microscopy_method TEXT, 
  microscopy_notes TEXT, preservation TEXT, 
  volume_concentrated TEXT, volume_counted TEXT, 
  magnification TEXT, counting_chamber TEXT, 
  notes TEXT
); CREATE TEMP TABLE curated_uncertainty (
  original_id INTEGER PRIMARY KEY NOT NULL, 
  uncertain INTEGER NOT NULL CHECK (
    uncertain IN (0, 1)
  )
); CREATE TEMP TABLE curated_notes (
  original_id INTEGER PRIMARY KEY NOT NULL, 
  abundance_notes TEXT, coordinates_notes TEXT, 
  event_datetime_notes TEXT, depth_notes TEXT, 
  taxonomy_notes TEXT, event_notes TEXT, 
  sample_notes TEXT, taxon_record_notes TEXT
); 
-- Uncomment only the imports for supplied inputs. Empty curation tables mean no overrides.
-- .import --csv --skip 1 --schema temp source_metadata.csv curated_source
-- .import --csv --skip 1 --schema temp sampling_protocol_curated.csv curated_sampling
-- .import --csv --skip 1 --schema temp taxonomic_uncertainty_curated.csv curated_uncertainty
-- .import --csv --skip 1 --schema temp record_notes_curated.csv curated_notes
-- Import the dictionary only when metadata_dictionary is empty.
-- .import --csv --skip 1 data_dictionary.csv metadata_dictionary
UPDATE 
  curated_source 
SET 
  subsource = NULLIF(subsource, ''), 
  metadata = NULLIF(metadata, ''), 
  url = NULLIF(url, ''), 
  notes = NULLIF(notes, ''); CREATE UNIQUE INDEX temp.uq_curated_source ON curated_source (
  source, 
  COALESCE(subsource, '')
); 
UPDATE 
  curated_sampling 
SET 
  count_method = NULLIF(count_method, ''), 
  enumeration_method = NULLIF(enumeration_method, ''), 
  sampling_method = NULLIF(sampling_method, ''), 
  sampling_method_notes = NULLIF(sampling_method_notes, ''), 
  concentration_method = NULLIF(concentration_method, ''), 
  concentration_method_notes = NULLIF(concentration_method_notes, ''), 
  volume_sampled = NULLIF(volume_sampled, ''), 
  microscopy_method = NULLIF(microscopy_method, ''), 
  microscopy_notes = NULLIF(microscopy_notes, ''), 
  preservation = NULLIF(preservation, ''), 
  volume_concentrated = NULLIF(volume_concentrated, ''), 
  volume_counted = NULLIF(volume_counted, ''), 
  magnification = NULLIF(magnification, ''), 
  counting_chamber = NULLIF(counting_chamber, ''), 
  notes = NULLIF(notes, ''); 
UPDATE 
  curated_notes 
SET 
  abundance_notes = NULLIF(abundance_notes, ''), 
  coordinates_notes = NULLIF(coordinates_notes, ''), 
  event_datetime_notes = NULLIF(event_datetime_notes, ''), 
  depth_notes = NULLIF(depth_notes, ''), 
  taxonomy_notes = NULLIF(taxonomy_notes, ''), 
  event_notes = NULLIF(event_notes, ''), 
  sample_notes = NULLIF(sample_notes, ''), 
  taxon_record_notes = NULLIF(taxon_record_notes, ''); 
UPDATE 
  metadata_dictionary 
SET 
  possible_values = NULLIF(possible_values, ''), 
  missing_value_interpretation = NULLIF(
    missing_value_interpretation, ''
  ), 
  notes = NULLIF(notes, ''); CREATE TEMP TABLE migration_assertion (
  check_name TEXT PRIMARY KEY, 
  failures INTEGER NOT NULL CHECK (failures = 0)
); INSERT INTO migration_assertion 
SELECT 
  'curated_sampling_source_ids', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      curated_sampling AS c 
    WHERE 
      NOT EXISTS (
        SELECT 
          1 
        FROM 
          abundance_full AS f 
        WHERE 
          f.id = c.original_id
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'curated_uncertainty_source_ids', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      curated_uncertainty AS c 
    WHERE 
      NOT EXISTS (
        SELECT 
          1 
        FROM 
          abundance_full AS f 
        WHERE 
          f.id = c.original_id
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'curated_notes_source_ids', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      curated_notes AS c 
    WHERE 
      NOT EXISTS (
        SELECT 
          1 
        FROM 
          abundance_full AS f 
        WHERE 
          f.id = c.original_id
      )
  ); 
-----------------------------------------
-- 3. Filter observations and curation --
-----------------------------------------
CREATE TABLE 
  abundance_full_pre AS 
SELECT 
  f.id AS id, 
  f.latitude AS latitude, 
  f.longitude AS longitude, 
  NULLIF(f.coord_precision, '') AS coord_precision, 
  f.year AS year, 
  f.month AS month, 
  f.day AS day, 
  NULLIF(f.datetime, '') AS datetime, 
  NULLIF(f.date_precision, '') AS date_precision, 
  f.day_is_imputed AS day_is_imputed, 
  f.depth AS depth, 
  NULLIF(f.depth_label, '') AS depth_label, 
  f.depth_label_min AS depth_label_min, 
  f.depth_label_max AS depth_label_max, 
  f.t_aphiaid AS t_aphiaid, 
  NULLIF(f.t_kingdom, '') AS t_kingdom, 
  NULLIF(f.t_subkingdom, '') AS t_subkingdom, 
  NULLIF(f.t_infrakingdom, '') AS t_infrakingdom, 
  NULLIF(f.t_phylum, '') AS t_phylum, 
  NULLIF(f.t_subphylum, '') AS t_subphylum, 
  NULLIF(f.t_infraphylum, '') AS t_infraphylum, 
  NULLIF(f.t_parvphylum, '') AS t_parvphylum, 
  NULLIF(f.t_gigaclass, '') AS t_gigaclass, 
  NULLIF(f.t_superclass, '') AS t_superclass, 
  NULLIF(f.t_class, '') AS t_class, 
  NULLIF(f.t_subclass, '') AS t_subclass, 
  NULLIF(f.t_infraclass, '') AS t_infraclass, 
  NULLIF(f.t_subterclass, '') AS t_subterclass, 
  NULLIF(f.t_superorder, '') AS t_superorder, 
  NULLIF(f.t_order, '') AS t_order, 
  NULLIF(f.t_suborder, '') AS t_suborder, 
  NULLIF(f.t_infraorder, '') AS t_infraorder, 
  NULLIF(f.t_parvorder, '') AS t_parvorder, 
  NULLIF(f.t_superfamily, '') AS t_superfamily, 
  NULLIF(f.t_family, '') AS t_family, 
  NULLIF(f.t_subfamily, '') AS t_subfamily, 
  NULLIF(f.t_tribe, '') AS t_tribe, 
  NULLIF(f.t_genus, '') AS t_genus, 
  NULLIF(f.t_subgenus, '') AS t_subgenus, 
  NULLIF(f.t_section, '') AS t_section, 
  NULLIF(f.t_subsection, '') AS t_subsection, 
  NULLIF(f.t_species, '') AS t_species, 
  NULLIF(f.t_subspecies, '') AS t_subspecies, 
  NULLIF(f.t_variety, '') AS t_variety, 
  NULLIF(f.t_forma, '') AS t_forma, 
  NULLIF(f.t_rank, '') AS t_rank, 
  NULLIF(f.t_taxonomic_status, '') AS t_taxonomic_status, 
  NULLIF(
    CASE 
      WHEN cp.original_id IS NOT NULL THEN cp.count_method ELSE f.s_count_method
    END, 
    ''
  ) AS s_count_method, 
  NULLIF(
    CASE 
      WHEN cp.original_id IS NOT NULL THEN cp.enumeration_method ELSE f.s_enumeration_method
    END, 
    ''
  ) AS s_enumeration_method, 
  NULLIF(
    CASE 
      WHEN cp.original_id IS NOT NULL THEN cp.sampling_method ELSE f.s_sampling_method
    END, 
    ''
  ) AS s_sampling_method, 
  NULLIF(
    CASE 
      WHEN cp.original_id IS NOT NULL THEN cp.sampling_method_notes ELSE f.s_sampling_method_notes
    END, 
    ''
  ) AS s_sampling_method_notes, 
  NULLIF(
    CASE 
      WHEN cp.original_id IS NOT NULL THEN cp.concentration_method ELSE f.s_concentration_method
    END, 
    ''
  ) AS s_concentration_method, 
  NULLIF(
    CASE 
      WHEN cp.original_id IS NOT NULL THEN cp.concentration_method_notes ELSE f.s_concentration_method_notes
    END, 
    ''
  ) AS s_concentration_method_notes, 
  NULLIF(
    CASE 
      WHEN cp.original_id IS NOT NULL THEN cp.volume_sampled ELSE f.s_volume_sampled
    END, 
    ''
  ) AS s_volume_sampled, 
  NULLIF(
    CASE 
      WHEN cp.original_id IS NOT NULL THEN cp.microscopy_method ELSE f.s_microscopy_method
    END, 
    ''
  ) AS s_microscopy_method, 
  NULLIF(
    CASE 
      WHEN cp.original_id IS NOT NULL THEN cp.microscopy_notes ELSE f.s_microscopy_notes
    END, 
    ''
  ) AS s_microscopy_notes, 
  NULLIF(
    CASE 
      WHEN cp.original_id IS NOT NULL THEN cp.preservation ELSE f.s_preservation
    END, 
    ''
  ) AS s_preservation, 
  NULLIF(
    CASE 
      WHEN cp.original_id IS NOT NULL THEN cp.volume_concentrated ELSE f.s_volume_concentrated
    END, 
    ''
  ) AS s_volume_concentrated, 
  NULLIF(
    CASE 
      WHEN cp.original_id IS NOT NULL THEN cp.volume_counted ELSE f.s_volume_counted
    END, 
    ''
  ) AS s_volume_counted, 
  NULLIF(
    CASE 
      WHEN cp.original_id IS NOT NULL THEN cp.magnification ELSE f.s_magnification
    END, 
    ''
  ) AS s_magnification, 
  NULLIF(
    CASE 
      WHEN cp.original_id IS NOT NULL THEN cp.counting_chamber ELSE f.s_counting_chamber
    END, 
    ''
  ) AS s_counting_chamber, 
  NULLIF(f.source, '') AS source, 
  NULLIF(f.subsource, '') AS subsource, 
  NULLIF(f.t_taxon, '') AS t_taxon, 
  NULLIF(f.t_taxon_clean, '') AS t_taxon_clean, 
  NULLIF(f.t_matched_name, '') AS t_matched_name, 
  f.t_matched_aphiaid AS t_matched_aphiaid, 
  NULLIF(f.t_accepted_name, '') AS t_accepted_name, 
  f.t_accepted_aphiaid AS t_accepted_aphiaid, 
  NULLIF(f.t_resolution_method, '') AS t_resolution_method, 
  CASE 
    WHEN cu.original_id IS NOT NULL THEN cu.uncertain ELSE f.n_uncertain
  END AS n_uncertain, 
  f.n_ltr AS n_ltr, 
  NULLIF(f.n_aggregates, '') AS n_aggregates, 
  NULLIF(f.n_size_info, '') AS n_size_info, 
  NULLIF(f.n_morphotypes, '') AS n_morphotypes, 
  NULLIF(f.n_life_stages, '') AS n_life_stages, 
  NULLIF(f.n_infraspecific_taxa, '') AS n_infraspecific_taxa, 
  NULLIF(f.n_trophic, '') AS n_trophic, 
  NULLIF(f.n_local_authority, '') AS n_local_authority, 
  f.abundance AS abundance, 
  cn.abundance_notes AS abundance_notes, 
  cn.coordinates_notes AS coordinates_notes, 
  cn.event_datetime_notes AS event_datetime_notes, 
  cn.depth_notes AS depth_notes, 
  cn.taxonomy_notes AS taxonomy_notes, 
  cn.event_notes AS event_notes, 
  cn.sample_notes AS sample_notes, 
  cn.taxon_record_notes AS taxon_record_notes, 
  cp.notes AS sampling_protocol_notes 
FROM 
  abundance_full AS f 
  LEFT JOIN curated_sampling AS cp ON cp.original_id = f.id 
  LEFT JOIN curated_uncertainty AS cu ON cu.original_id = f.id 
  LEFT JOIN curated_notes AS cn ON cn.original_id = f.id 
WHERE 
  f.n_flag_for_removal IS NOT 1; CREATE UNIQUE INDEX ix_abundance_full_pre_original_id ON abundance_full_pre(id); INSERT INTO migration_assertion 
SELECT 
  'source_ids_not_null', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      abundance_full_pre 
    WHERE 
      id IS NULL
  ); INSERT INTO migration_assertion 
SELECT 
  'retained_source_rows', 
  (
    SELECT 
      ABS(
        (
          SELECT 
            COUNT(*) 
          FROM 
            abundance_full 
          WHERE 
            n_flag_for_removal IS NOT 1
        ) - (
          SELECT 
            COUNT(*) 
          FROM 
            abundance_full_pre
        )
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'time_format', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      abundance_full_pre 
    WHERE 
      datetime IS NULL 
      OR length(datetime) <> 8 
      OR datetime NOT GLOB '[0-2][0-9]:[0-5][0-9]:[0-5][0-9]' 
      OR CAST(
        substr(datetime, 1, 2) AS INTEGER
      ) > 23
  ); 
CREATE INDEX 
  ix_abundance_full_pre_source ON abundance_full_pre(source, subsource); INSERT INTO migration_assertion 
SELECT 
  'curated_source_coverage', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      (
        SELECT 
          DISTINCT source, 
          COALESCE(subsource, '') AS subsource 
        FROM 
          abundance_full_pre
      ) AS p 
    WHERE 
      EXISTS (
        SELECT 
          1 
        FROM 
          curated_source
      ) 
      AND NOT EXISTS (
        SELECT 
          1 
        FROM 
          curated_source AS c 
        WHERE 
          c.source = p.source 
          AND COALESCE(c.subsource, '') = p.subsource
      )
  ); 
-------------------------------------
-- 4. Deterministic relational IDs --
-------------------------------------
ALTER TABLE 
  abundance_full_pre 
ADD 
  COLUMN coordinates_id INTEGER; 
ALTER TABLE 
  abundance_full_pre 
ADD 
  COLUMN event_datetime_id INTEGER; 
ALTER TABLE 
  abundance_full_pre 
ADD 
  COLUMN depth_id INTEGER; 
ALTER TABLE 
  abundance_full_pre 
ADD 
  COLUMN taxonomy_id INTEGER; 
ALTER TABLE 
  abundance_full_pre 
ADD 
  COLUMN sampling_protocol_id INTEGER; 
ALTER TABLE 
  abundance_full_pre 
ADD 
  COLUMN source_id INTEGER; 
ALTER TABLE 
  abundance_full_pre 
ADD 
  COLUMN event_id INTEGER; 
ALTER TABLE 
  abundance_full_pre 
ADD 
  COLUMN sample_id INTEGER; 
ALTER TABLE 
  abundance_full_pre 
ADD 
  COLUMN taxon_record_id INTEGER; 
ALTER TABLE 
  abundance_full_pre 
ADD 
  COLUMN abundance_id INTEGER; CREATE TEMP TABLE id_map AS 
SELECT 
  id AS original_id, 
  DENSE_RANK() OVER (
    ORDER BY 
      latitude, 
      longitude, 
      COALESCE(coord_precision, '')
  ) AS relational_id 
FROM 
  abundance_full_pre; CREATE UNIQUE INDEX temp.ix_id_map ON id_map(original_id); 
UPDATE 
  abundance_full_pre AS p 
SET 
  coordinates_id = (
    SELECT 
      relational_id 
    FROM 
      id_map AS m 
    WHERE 
      m.original_id = p.id
  ); 
DROP 
  TABLE id_map; CREATE TEMP TABLE id_map AS 
SELECT 
  id AS original_id, 
  DENSE_RANK() OVER (
    ORDER BY 
      year, 
      month, 
      day, 
      COALESCE(datetime, ''), 
      COALESCE(date_precision, ''), 
      day_is_imputed
  ) AS relational_id 
FROM 
  abundance_full_pre; CREATE UNIQUE INDEX temp.ix_id_map ON id_map(original_id); 
UPDATE 
  abundance_full_pre AS p 
SET 
  event_datetime_id = (
    SELECT 
      relational_id 
    FROM 
      id_map AS m 
    WHERE 
      m.original_id = p.id
  ); 
DROP 
  TABLE id_map; CREATE TEMP TABLE id_map AS 
SELECT 
  id AS original_id, 
  DENSE_RANK() OVER (
    ORDER BY 
      depth, 
      COALESCE(depth_label, ''), 
      COALESCE(depth_label_min, -1.0), 
      COALESCE(depth_label_max, -1.0)
  ) AS relational_id 
FROM 
  abundance_full_pre; CREATE UNIQUE INDEX temp.ix_id_map ON id_map(original_id); 
UPDATE 
  abundance_full_pre AS p 
SET 
  depth_id = (
    SELECT 
      relational_id 
    FROM 
      id_map AS m 
    WHERE 
      m.original_id = p.id
  ); 
DROP 
  TABLE id_map; CREATE TEMP TABLE id_map AS 
SELECT 
  id AS original_id, 
  DENSE_RANK() OVER (
    ORDER BY 
      t_aphiaid
  ) AS relational_id 
FROM 
  abundance_full_pre; CREATE UNIQUE INDEX temp.ix_id_map ON id_map(original_id); 
UPDATE 
  abundance_full_pre AS p 
SET 
  taxonomy_id = (
    SELECT 
      relational_id 
    FROM 
      id_map AS m 
    WHERE 
      m.original_id = p.id
  ); 
DROP 
  TABLE id_map; CREATE TEMP TABLE id_map AS 
SELECT 
  id AS original_id, 
  DENSE_RANK() OVER (
    ORDER BY 
      COALESCE(s_count_method, ''), 
      COALESCE(s_enumeration_method, ''), 
      COALESCE(s_sampling_method, ''), 
      COALESCE(s_sampling_method_notes, ''), 
      COALESCE(s_concentration_method, ''), 
      COALESCE(
        s_concentration_method_notes, ''
      ), 
      COALESCE(s_volume_sampled, ''), 
      COALESCE(s_microscopy_method, ''), 
      COALESCE(s_microscopy_notes, ''), 
      COALESCE(s_preservation, ''), 
      COALESCE(s_volume_concentrated, ''), 
      COALESCE(s_volume_counted, ''), 
      COALESCE(s_magnification, ''), 
      COALESCE(s_counting_chamber, '')
  ) AS relational_id 
FROM 
  abundance_full_pre; CREATE UNIQUE INDEX temp.ix_id_map ON id_map(original_id); 
UPDATE 
  abundance_full_pre AS p 
SET 
  sampling_protocol_id = (
    SELECT 
      relational_id 
    FROM 
      id_map AS m 
    WHERE 
      m.original_id = p.id
  ); 
DROP 
  TABLE id_map; CREATE TEMP TABLE id_map AS 
SELECT 
  id AS original_id, 
  DENSE_RANK() OVER (
    ORDER BY 
      source, 
      COALESCE(subsource, '')
  ) AS relational_id 
FROM 
  abundance_full_pre; CREATE UNIQUE INDEX temp.ix_id_map ON id_map(original_id); 
UPDATE 
  abundance_full_pre AS p 
SET 
  source_id = (
    SELECT 
      relational_id 
    FROM 
      id_map AS m 
    WHERE 
      m.original_id = p.id
  ); 
DROP 
  TABLE id_map; CREATE TEMP TABLE id_map AS 
SELECT 
  id AS original_id, 
  DENSE_RANK() OVER (
    ORDER BY 
      event_datetime_id, 
      coordinates_id, 
      source_id
  ) AS relational_id 
FROM 
  abundance_full_pre; CREATE UNIQUE INDEX temp.ix_id_map ON id_map(original_id); 
UPDATE 
  abundance_full_pre AS p 
SET 
  event_id = (
    SELECT 
      relational_id 
    FROM 
      id_map AS m 
    WHERE 
      m.original_id = p.id
  ); 
DROP 
  TABLE id_map; CREATE TEMP TABLE id_map AS 
SELECT 
  id AS original_id, 
  DENSE_RANK() OVER (
    ORDER BY 
      event_id, 
      depth_id
  ) AS relational_id 
FROM 
  abundance_full_pre; CREATE UNIQUE INDEX temp.ix_id_map ON id_map(original_id); 
UPDATE 
  abundance_full_pre AS p 
SET 
  sample_id = (
    SELECT 
      relational_id 
    FROM 
      id_map AS m 
    WHERE 
      m.original_id = p.id
  ); 
DROP 
  TABLE id_map; CREATE TEMP TABLE id_map AS 
SELECT 
  id AS original_id, 
  DENSE_RANK() OVER (
    ORDER BY 
      taxonomy_id, 
      COALESCE(t_taxon, ''), 
      COALESCE(t_taxon_clean, ''), 
      COALESCE(t_matched_name, ''), 
      COALESCE(t_matched_aphiaid, -1), 
      COALESCE(t_accepted_name, ''), 
      COALESCE(t_accepted_aphiaid, -1), 
      COALESCE(t_resolution_method, ''), 
      COALESCE(n_uncertain, -1), 
      COALESCE(n_ltr, -1), 
      COALESCE(n_aggregates, ''), 
      COALESCE(n_size_info, ''), 
      COALESCE(n_morphotypes, ''), 
      COALESCE(n_life_stages, ''), 
      COALESCE(n_infraspecific_taxa, ''), 
      COALESCE(n_trophic, ''), 
      COALESCE(n_local_authority, '')
  ) AS relational_id 
FROM 
  abundance_full_pre; CREATE UNIQUE INDEX temp.ix_id_map ON id_map(original_id); 
UPDATE 
  abundance_full_pre AS p 
SET 
  taxon_record_id = (
    SELECT 
      relational_id 
    FROM 
      id_map AS m 
    WHERE 
      m.original_id = p.id
  ); 
DROP 
  TABLE id_map; CREATE TEMP TABLE id_map AS 
SELECT 
  id AS original_id, 
  DENSE_RANK() OVER (
    ORDER BY 
      sample_id, 
      taxon_record_id, 
      sampling_protocol_id, 
      abundance
  ) AS relational_id 
FROM 
  abundance_full_pre; CREATE UNIQUE INDEX temp.ix_id_map ON id_map(original_id); 
UPDATE 
  abundance_full_pre AS p 
SET 
  abundance_id = (
    SELECT 
      relational_id 
    FROM 
      id_map AS m 
    WHERE 
      m.original_id = p.id
  ); 
DROP 
  TABLE id_map; 
---------------------------------------------
-- 5. Load dimensions and dependent tables --
---------------------------------------------
INSERT INTO migration_assertion 
SELECT 
  'coordinates_single_representation', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      (
        SELECT 
          coordinates_id 
        FROM 
          (
            SELECT 
              DISTINCT coordinates_id, 
              latitude, 
              longitude, 
              coord_precision 
            FROM 
              abundance_full_pre
          ) 
        GROUP BY 
          coordinates_id 
        HAVING 
          COUNT(*) > 1
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'coordinates_consistent_notes', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      (
        SELECT 
          coordinates_id 
        FROM 
          abundance_full_pre 
        GROUP BY 
          coordinates_id 
        HAVING 
          COUNT(DISTINCT coordinates_notes) > 1
      )
  ); INSERT INTO coordinates (
  id, latitude, longitude, coord_precision, 
  notes
) 
SELECT 
  coordinates_id, 
  latitude, 
  longitude, 
  coord_precision, 
  MAX(coordinates_notes) 
FROM 
  abundance_full_pre 
GROUP BY 
  coordinates_id; INSERT INTO migration_assertion 
SELECT 
  'event_datetime_single_representation', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      (
        SELECT 
          event_datetime_id 
        FROM 
          (
            SELECT 
              DISTINCT event_datetime_id, 
              year, 
              month, 
              day, 
              datetime, 
              date_precision, 
              day_is_imputed 
            FROM 
              abundance_full_pre
          ) 
        GROUP BY 
          event_datetime_id 
        HAVING 
          COUNT(*) > 1
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'event_datetime_consistent_notes', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      (
        SELECT 
          event_datetime_id 
        FROM 
          abundance_full_pre 
        GROUP BY 
          event_datetime_id 
        HAVING 
          COUNT(DISTINCT event_datetime_notes) > 1
      )
  ); INSERT INTO event_datetime (
  id, year, month, day, time_value, date_precision, 
  day_is_imputed, notes
) 
SELECT 
  event_datetime_id, 
  year, 
  month, 
  day, 
  datetime, 
  date_precision, 
  day_is_imputed, 
  MAX(event_datetime_notes) 
FROM 
  abundance_full_pre 
GROUP BY 
  event_datetime_id; INSERT INTO migration_assertion 
SELECT 
  'depth_single_representation', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      (
        SELECT 
          depth_id 
        FROM 
          (
            SELECT 
              DISTINCT depth_id, 
              depth, 
              depth_label, 
              depth_label_min, 
              depth_label_max 
            FROM 
              abundance_full_pre
          ) 
        GROUP BY 
          depth_id 
        HAVING 
          COUNT(*) > 1
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'depth_consistent_notes', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      (
        SELECT 
          depth_id 
        FROM 
          abundance_full_pre 
        GROUP BY 
          depth_id 
        HAVING 
          COUNT(DISTINCT depth_notes) > 1
      )
  ); INSERT INTO depth (
  id, depth, depth_label, depth_label_min, 
  depth_label_max, notes
) 
SELECT 
  depth_id, 
  depth, 
  depth_label, 
  depth_label_min, 
  depth_label_max, 
  MAX(depth_notes) 
FROM 
  abundance_full_pre 
GROUP BY 
  depth_id; INSERT INTO migration_assertion 
SELECT 
  'taxonomy_single_representation', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      (
        SELECT 
          taxonomy_id 
        FROM 
          (
            SELECT 
              DISTINCT taxonomy_id, 
              t_aphiaid, 
              t_kingdom, 
              t_subkingdom, 
              t_infrakingdom, 
              t_phylum, 
              t_subphylum, 
              t_infraphylum, 
              t_parvphylum, 
              t_gigaclass, 
              t_superclass, 
              t_class, 
              t_subclass, 
              t_infraclass, 
              t_subterclass, 
              t_superorder, 
              t_order, 
              t_suborder, 
              t_infraorder, 
              t_parvorder, 
              t_superfamily, 
              t_family, 
              t_subfamily, 
              t_tribe, 
              t_genus, 
              t_subgenus, 
              t_section, 
              t_subsection, 
              t_species, 
              t_subspecies, 
              t_variety, 
              t_forma, 
              t_rank, 
              t_taxonomic_status 
            FROM 
              abundance_full_pre
          ) 
        GROUP BY 
          taxonomy_id 
        HAVING 
          COUNT(*) > 1
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'taxonomy_consistent_notes', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      (
        SELECT 
          taxonomy_id 
        FROM 
          abundance_full_pre 
        GROUP BY 
          taxonomy_id 
        HAVING 
          COUNT(DISTINCT taxonomy_notes) > 1
      )
  ); INSERT INTO taxonomy (
  id, aphiaid, kingdom, subkingdom, infrakingdom, 
  phylum, subphylum, infraphylum, parvphylum, 
  gigaclass, superclass, class, subclass, 
  infraclass, subterclass, superorder, 
  taxon_order, suborder, infraorder, 
  parvorder, superfamily, family, subfamily, 
  tribe, genus, subgenus, section, subsection, 
  species, subspecies, variety, forma, 
  rank, taxonomic_status, notes
) 
SELECT 
  taxonomy_id, 
  t_aphiaid, 
  t_kingdom, 
  t_subkingdom, 
  t_infrakingdom, 
  t_phylum, 
  t_subphylum, 
  t_infraphylum, 
  t_parvphylum, 
  t_gigaclass, 
  t_superclass, 
  t_class, 
  t_subclass, 
  t_infraclass, 
  t_subterclass, 
  t_superorder, 
  t_order, 
  t_suborder, 
  t_infraorder, 
  t_parvorder, 
  t_superfamily, 
  t_family, 
  t_subfamily, 
  t_tribe, 
  t_genus, 
  t_subgenus, 
  t_section, 
  t_subsection, 
  t_species, 
  t_subspecies, 
  t_variety, 
  t_forma, 
  t_rank, 
  t_taxonomic_status, 
  MAX(taxonomy_notes) 
FROM 
  abundance_full_pre 
GROUP BY 
  taxonomy_id; INSERT INTO migration_assertion 
SELECT 
  'sampling_protocol_single_representation', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      (
        SELECT 
          sampling_protocol_id 
        FROM 
          (
            SELECT 
              DISTINCT sampling_protocol_id, 
              s_count_method, 
              s_enumeration_method, 
              s_sampling_method, 
              s_sampling_method_notes, 
              s_concentration_method, 
              s_concentration_method_notes, 
              s_volume_sampled, 
              s_microscopy_method, 
              s_microscopy_notes, 
              s_preservation, 
              s_volume_concentrated, 
              s_volume_counted, 
              s_magnification, 
              s_counting_chamber 
            FROM 
              abundance_full_pre
          ) 
        GROUP BY 
          sampling_protocol_id 
        HAVING 
          COUNT(*) > 1
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'sampling_protocol_consistent_notes', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      (
        SELECT 
          sampling_protocol_id 
        FROM 
          abundance_full_pre 
        GROUP BY 
          sampling_protocol_id 
        HAVING 
          COUNT(
            DISTINCT sampling_protocol_notes
          ) > 1
      )
  ); INSERT INTO sampling_protocol (
  id, count_method, enumeration_method, 
  sampling_method, sampling_method_notes, 
  concentration_method, concentration_method_notes, 
  volume_sampled, microscopy_method, 
  microscopy_notes, preservation, 
  volume_concentrated, volume_counted, 
  magnification, counting_chamber, 
  notes
) 
SELECT 
  sampling_protocol_id, 
  s_count_method, 
  s_enumeration_method, 
  s_sampling_method, 
  s_sampling_method_notes, 
  s_concentration_method, 
  s_concentration_method_notes, 
  s_volume_sampled, 
  s_microscopy_method, 
  s_microscopy_notes, 
  s_preservation, 
  s_volume_concentrated, 
  s_volume_counted, 
  s_magnification, 
  s_counting_chamber, 
  MAX(sampling_protocol_notes) 
FROM 
  abundance_full_pre 
GROUP BY 
  sampling_protocol_id; INSERT INTO migration_assertion 
SELECT 
  'source_single_representation', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      (
        SELECT 
          source_id 
        FROM 
          (
            SELECT 
              DISTINCT source_id, 
              source, 
              subsource 
            FROM 
              abundance_full_pre
          ) 
        GROUP BY 
          source_id 
        HAVING 
          COUNT(*) > 1
      )
  ); INSERT INTO source (
  id, source, subsource, metadata, url, 
  notes
) 
SELECT 
  p.source_id, 
  p.source, 
  p.subsource, 
  c.metadata, 
  c.url, 
  c.notes 
FROM 
  (
    SELECT 
      DISTINCT source_id, 
      source, 
      subsource 
    FROM 
      abundance_full_pre
  ) AS p 
  LEFT JOIN curated_source AS c ON c.source = p.source 
  AND COALESCE(c.subsource, '') = COALESCE(p.subsource, ''); INSERT INTO migration_assertion 
SELECT 
  'event_single_representation', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      (
        SELECT 
          event_id 
        FROM 
          (
            SELECT 
              DISTINCT event_id, 
              event_datetime_id, 
              coordinates_id, 
              source_id 
            FROM 
              abundance_full_pre
          ) 
        GROUP BY 
          event_id 
        HAVING 
          COUNT(*) > 1
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'event_consistent_notes', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      (
        SELECT 
          event_id 
        FROM 
          abundance_full_pre 
        GROUP BY 
          event_id 
        HAVING 
          COUNT(DISTINCT event_notes) > 1
      )
  ); INSERT INTO event (
  id, event_datetime_id, coordinates_id, 
  source_id, notes
) 
SELECT 
  event_id, 
  event_datetime_id, 
  coordinates_id, 
  source_id, 
  MAX(event_notes) 
FROM 
  abundance_full_pre 
GROUP BY 
  event_id; INSERT INTO migration_assertion 
SELECT 
  'sample_single_representation', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      (
        SELECT 
          sample_id 
        FROM 
          (
            SELECT 
              DISTINCT sample_id, 
              event_id, 
              depth_id 
            FROM 
              abundance_full_pre
          ) 
        GROUP BY 
          sample_id 
        HAVING 
          COUNT(*) > 1
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'sample_consistent_notes', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      (
        SELECT 
          sample_id 
        FROM 
          abundance_full_pre 
        GROUP BY 
          sample_id 
        HAVING 
          COUNT(DISTINCT sample_notes) > 1
      )
  ); INSERT INTO sample (id, event_id, depth_id, notes) 
SELECT 
  sample_id, 
  event_id, 
  depth_id, 
  MAX(sample_notes) 
FROM 
  abundance_full_pre 
GROUP BY 
  sample_id; INSERT INTO migration_assertion 
SELECT 
  'taxon_record_single_representation', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      (
        SELECT 
          taxon_record_id 
        FROM 
          (
            SELECT 
              DISTINCT taxon_record_id, 
              taxonomy_id, 
              t_taxon, 
              t_taxon_clean, 
              t_matched_name, 
              t_matched_aphiaid, 
              t_accepted_name, 
              t_accepted_aphiaid, 
              t_resolution_method, 
              n_uncertain, 
              n_ltr, 
              n_aggregates, 
              n_size_info, 
              n_morphotypes, 
              n_life_stages, 
              n_infraspecific_taxa, 
              n_trophic, 
              n_local_authority 
            FROM 
              abundance_full_pre
          ) 
        GROUP BY 
          taxon_record_id 
        HAVING 
          COUNT(*) > 1
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'taxon_record_consistent_notes', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      (
        SELECT 
          taxon_record_id 
        FROM 
          abundance_full_pre 
        GROUP BY 
          taxon_record_id 
        HAVING 
          COUNT(DISTINCT taxon_record_notes) > 1
      )
  ); INSERT INTO taxon_record (
  id, taxonomy_id, taxon_raw, taxon_clean, 
  taxon_matched, aphiaid_matched, 
  taxon_accepted, aphiaid_accepted, 
  resolution_method, uncertain, ltr, 
  aggregates, size_info, morphotype, 
  life_stage, infraspecific_taxa, 
  trophic, local_authority, notes
) 
SELECT 
  taxon_record_id, 
  taxonomy_id, 
  t_taxon, 
  t_taxon_clean, 
  t_matched_name, 
  t_matched_aphiaid, 
  t_accepted_name, 
  t_accepted_aphiaid, 
  t_resolution_method, 
  n_uncertain, 
  n_ltr, 
  n_aggregates, 
  n_size_info, 
  n_morphotypes, 
  n_life_stages, 
  n_infraspecific_taxa, 
  n_trophic, 
  n_local_authority, 
  MAX(taxon_record_notes) 
FROM 
  abundance_full_pre 
GROUP BY 
  taxon_record_id; INSERT INTO migration_assertion 
SELECT 
  'abundance_single_representation', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      (
        SELECT 
          abundance_id 
        FROM 
          (
            SELECT 
              DISTINCT abundance_id, 
              sample_id, 
              taxon_record_id, 
              sampling_protocol_id, 
              abundance 
            FROM 
              abundance_full_pre
          ) 
        GROUP BY 
          abundance_id 
        HAVING 
          COUNT(*) > 1
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'abundance_consistent_notes', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      (
        SELECT 
          abundance_id 
        FROM 
          abundance_full_pre 
        GROUP BY 
          abundance_id 
        HAVING 
          COUNT(DISTINCT abundance_notes) > 1
      )
  ); INSERT INTO abundance (
  id, sample_id, taxon_record_id, sampling_protocol_id, 
  abundance, notes
) 
SELECT 
  abundance_id, 
  sample_id, 
  taxon_record_id, 
  sampling_protocol_id, 
  abundance, 
  MAX(abundance_notes) 
FROM 
  abundance_full_pre 
GROUP BY 
  abundance_id; 
-----------------------------------
-- 6. Denormalized view creation --
-----------------------------------
CREATE VIEW 
  v_abundance_flat AS 
SELECT 
  ab.id AS abundance_id, 
  src.id AS source_id, 
  src.source, 
  src.subsource, 
  src.metadata AS source_metadata, 
  src.url AS source_url, 
  ev.id AS event_id, 
  edt.id AS event_datetime_id, 
  edt.year, 
  edt.month, 
  edt.day, 
  edt.time_value, 
  edt.date_precision, 
  edt.day_is_imputed, 
  c.id AS coordinates_id, 
  c.latitude, 
  c.longitude, 
  c.coord_precision, 
  sm.id AS sample_id, 
  d.id AS depth_id, 
  d.depth, 
  d.depth_label, 
  d.depth_label_min, 
  d.depth_label_max, 
  ab.abundance, 
  ab.notes AS abundance_notes, 
  sp.id AS sampling_protocol_id, 
  sp.count_method, 
  sp.enumeration_method, 
  sp.sampling_method, 
  sp.sampling_method_notes, 
  sp.concentration_method, 
  sp.concentration_method_notes, 
  sp.volume_sampled, 
  sp.microscopy_method, 
  sp.microscopy_notes, 
  sp.preservation, 
  sp.volume_concentrated, 
  sp.volume_counted, 
  sp.magnification, 
  sp.counting_chamber, 
  tr.id AS taxon_record_id, 
  tr.taxon_raw, 
  tr.taxon_clean, 
  tr.taxon_matched, 
  tr.aphiaid_matched, 
  tr.taxon_accepted, 
  tr.aphiaid_accepted, 
  tr.resolution_method, 
  tr.uncertain, 
  tr.ltr, 
  tr.aggregates, 
  tr.size_info, 
  tr.morphotype, 
  tr.life_stage, 
  tr.infraspecific_taxa, 
  tr.trophic, 
  tr.local_authority, 
  tx.id AS taxonomy_id, 
  tx.aphiaid, 
  tx.kingdom, 
  tx.subkingdom, 
  tx.infrakingdom, 
  tx.phylum, 
  tx.subphylum, 
  tx.infraphylum, 
  tx.parvphylum, 
  tx.gigaclass, 
  tx.superclass, 
  tx.class, 
  tx.subclass, 
  tx.infraclass, 
  tx.subterclass, 
  tx.superorder, 
  tx.taxon_order, 
  tx.suborder, 
  tx.infraorder, 
  tx.parvorder, 
  tx.superfamily, 
  tx.family, 
  tx.subfamily, 
  tx.tribe, 
  tx.genus, 
  tx.subgenus, 
  tx.section, 
  tx.subsection, 
  tx.species, 
  tx.subspecies, 
  tx.variety, 
  tx.forma, 
  tx.rank, 
  tx.taxonomic_status 
FROM 
  abundance ab 
  JOIN sample sm ON sm.id = ab.sample_id 
  JOIN event ev ON ev.id = sm.event_id 
  JOIN depth d ON d.id = sm.depth_id 
  JOIN event_datetime edt ON edt.id = ev.event_datetime_id 
  JOIN coordinates c ON c.id = ev.coordinates_id 
  JOIN source src ON src.id = ev.source_id 
  JOIN taxon_record tr ON tr.id = ab.taxon_record_id 
  JOIN taxonomy tx ON tx.id = tr.taxonomy_id 
  JOIN sampling_protocol sp ON sp.id = ab.sampling_protocol_id; 
-------------------------------------------------
-- 7. Validate counts and reconstructed values --
-------------------------------------------------
INSERT INTO migration_assertion 
SELECT 
  'coordinates_row_count', 
  (
    SELECT 
      ABS(
        (
          SELECT 
            COUNT(*) 
          FROM 
            coordinates
        ) - (
          SELECT 
            COUNT(DISTINCT coordinates_id) 
          FROM 
            abundance_full_pre
        )
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'coordinates_values', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      abundance_full_pre AS p 
      LEFT JOIN coordinates AS r ON r.id = p.coordinates_id 
    WHERE 
      r.id IS NULL 
      OR r.latitude IS NOT p.latitude 
      OR r.longitude IS NOT p.longitude 
      OR r.coord_precision IS NOT p.coord_precision 
      OR (
        p.coordinates_notes IS NOT NULL 
        AND r.notes IS NOT p.coordinates_notes
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'event_datetime_row_count', 
  (
    SELECT 
      ABS(
        (
          SELECT 
            COUNT(*) 
          FROM 
            event_datetime
        ) - (
          SELECT 
            COUNT(DISTINCT event_datetime_id) 
          FROM 
            abundance_full_pre
        )
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'event_datetime_values', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      abundance_full_pre AS p 
      LEFT JOIN event_datetime AS r ON r.id = p.event_datetime_id 
    WHERE 
      r.id IS NULL 
      OR r.year IS NOT p.year 
      OR r.month IS NOT p.month 
      OR r.day IS NOT p.day 
      OR r.time_value IS NOT p.datetime 
      OR r.date_precision IS NOT p.date_precision 
      OR r.day_is_imputed IS NOT p.day_is_imputed 
      OR (
        p.event_datetime_notes IS NOT NULL 
        AND r.notes IS NOT p.event_datetime_notes
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'depth_row_count', 
  (
    SELECT 
      ABS(
        (
          SELECT 
            COUNT(*) 
          FROM 
            depth
        ) - (
          SELECT 
            COUNT(DISTINCT depth_id) 
          FROM 
            abundance_full_pre
        )
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'depth_values', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      abundance_full_pre AS p 
      LEFT JOIN depth AS r ON r.id = p.depth_id 
    WHERE 
      r.id IS NULL 
      OR r.depth IS NOT p.depth 
      OR r.depth_label IS NOT p.depth_label 
      OR r.depth_label_min IS NOT p.depth_label_min 
      OR r.depth_label_max IS NOT p.depth_label_max 
      OR (
        p.depth_notes IS NOT NULL 
        AND r.notes IS NOT p.depth_notes
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'taxonomy_row_count', 
  (
    SELECT 
      ABS(
        (
          SELECT 
            COUNT(*) 
          FROM 
            taxonomy
        ) - (
          SELECT 
            COUNT(DISTINCT taxonomy_id) 
          FROM 
            abundance_full_pre
        )
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'taxonomy_values', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      abundance_full_pre AS p 
      LEFT JOIN taxonomy AS r ON r.id = p.taxonomy_id 
    WHERE 
      r.id IS NULL 
      OR r.aphiaid IS NOT p.t_aphiaid 
      OR r.kingdom IS NOT p.t_kingdom 
      OR r.subkingdom IS NOT p.t_subkingdom 
      OR r.infrakingdom IS NOT p.t_infrakingdom 
      OR r.phylum IS NOT p.t_phylum 
      OR r.subphylum IS NOT p.t_subphylum 
      OR r.infraphylum IS NOT p.t_infraphylum 
      OR r.parvphylum IS NOT p.t_parvphylum 
      OR r.gigaclass IS NOT p.t_gigaclass 
      OR r.superclass IS NOT p.t_superclass 
      OR r.class IS NOT p.t_class 
      OR r.subclass IS NOT p.t_subclass 
      OR r.infraclass IS NOT p.t_infraclass 
      OR r.subterclass IS NOT p.t_subterclass 
      OR r.superorder IS NOT p.t_superorder 
      OR r.taxon_order IS NOT p.t_order 
      OR r.suborder IS NOT p.t_suborder 
      OR r.infraorder IS NOT p.t_infraorder 
      OR r.parvorder IS NOT p.t_parvorder 
      OR r.superfamily IS NOT p.t_superfamily 
      OR r.family IS NOT p.t_family 
      OR r.subfamily IS NOT p.t_subfamily 
      OR r.tribe IS NOT p.t_tribe 
      OR r.genus IS NOT p.t_genus 
      OR r.subgenus IS NOT p.t_subgenus 
      OR r.section IS NOT p.t_section 
      OR r.subsection IS NOT p.t_subsection 
      OR r.species IS NOT p.t_species 
      OR r.subspecies IS NOT p.t_subspecies 
      OR r.variety IS NOT p.t_variety 
      OR r.forma IS NOT p.t_forma 
      OR r.rank IS NOT p.t_rank 
      OR r.taxonomic_status IS NOT p.t_taxonomic_status 
      OR (
        p.taxonomy_notes IS NOT NULL 
        AND r.notes IS NOT p.taxonomy_notes
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'sampling_protocol_row_count', 
  (
    SELECT 
      ABS(
        (
          SELECT 
            COUNT(*) 
          FROM 
            sampling_protocol
        ) - (
          SELECT 
            COUNT(DISTINCT sampling_protocol_id) 
          FROM 
            abundance_full_pre
        )
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'sampling_protocol_values', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      abundance_full_pre AS p 
      LEFT JOIN sampling_protocol AS r ON r.id = p.sampling_protocol_id 
    WHERE 
      r.id IS NULL 
      OR r.count_method IS NOT p.s_count_method 
      OR r.enumeration_method IS NOT p.s_enumeration_method 
      OR r.sampling_method IS NOT p.s_sampling_method 
      OR r.sampling_method_notes IS NOT p.s_sampling_method_notes 
      OR r.concentration_method IS NOT p.s_concentration_method 
      OR r.concentration_method_notes IS NOT p.s_concentration_method_notes 
      OR r.volume_sampled IS NOT p.s_volume_sampled 
      OR r.microscopy_method IS NOT p.s_microscopy_method 
      OR r.microscopy_notes IS NOT p.s_microscopy_notes 
      OR r.preservation IS NOT p.s_preservation 
      OR r.volume_concentrated IS NOT p.s_volume_concentrated 
      OR r.volume_counted IS NOT p.s_volume_counted 
      OR r.magnification IS NOT p.s_magnification 
      OR r.counting_chamber IS NOT p.s_counting_chamber 
      OR (
        p.sampling_protocol_notes IS NOT NULL 
        AND r.notes IS NOT p.sampling_protocol_notes
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'source_row_count', 
  (
    SELECT 
      ABS(
        (
          SELECT 
            COUNT(*) 
          FROM 
            source
        ) - (
          SELECT 
            COUNT(DISTINCT source_id) 
          FROM 
            abundance_full_pre
        )
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'source_values', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      abundance_full_pre AS p 
      LEFT JOIN source AS r ON r.id = p.source_id 
    WHERE 
      r.id IS NULL 
      OR r.source IS NOT p.source 
      OR r.subsource IS NOT p.subsource
  ); INSERT INTO migration_assertion 
SELECT 
  'event_row_count', 
  (
    SELECT 
      ABS(
        (
          SELECT 
            COUNT(*) 
          FROM 
            event
        ) - (
          SELECT 
            COUNT(DISTINCT event_id) 
          FROM 
            abundance_full_pre
        )
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'event_values', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      abundance_full_pre AS p 
      LEFT JOIN event AS r ON r.id = p.event_id 
    WHERE 
      r.id IS NULL 
      OR r.event_datetime_id IS NOT p.event_datetime_id 
      OR r.coordinates_id IS NOT p.coordinates_id 
      OR r.source_id IS NOT p.source_id 
      OR (
        p.event_notes IS NOT NULL 
        AND r.notes IS NOT p.event_notes
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'sample_row_count', 
  (
    SELECT 
      ABS(
        (
          SELECT 
            COUNT(*) 
          FROM 
            sample
        ) - (
          SELECT 
            COUNT(DISTINCT sample_id) 
          FROM 
            abundance_full_pre
        )
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'sample_values', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      abundance_full_pre AS p 
      LEFT JOIN sample AS r ON r.id = p.sample_id 
    WHERE 
      r.id IS NULL 
      OR r.event_id IS NOT p.event_id 
      OR r.depth_id IS NOT p.depth_id 
      OR (
        p.sample_notes IS NOT NULL 
        AND r.notes IS NOT p.sample_notes
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'taxon_record_row_count', 
  (
    SELECT 
      ABS(
        (
          SELECT 
            COUNT(*) 
          FROM 
            taxon_record
        ) - (
          SELECT 
            COUNT(DISTINCT taxon_record_id) 
          FROM 
            abundance_full_pre
        )
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'taxon_record_values', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      abundance_full_pre AS p 
      LEFT JOIN taxon_record AS r ON r.id = p.taxon_record_id 
    WHERE 
      r.id IS NULL 
      OR r.taxonomy_id IS NOT p.taxonomy_id 
      OR r.taxon_raw IS NOT p.t_taxon 
      OR r.taxon_clean IS NOT p.t_taxon_clean 
      OR r.taxon_matched IS NOT p.t_matched_name 
      OR r.aphiaid_matched IS NOT p.t_matched_aphiaid 
      OR r.taxon_accepted IS NOT p.t_accepted_name 
      OR r.aphiaid_accepted IS NOT p.t_accepted_aphiaid 
      OR r.resolution_method IS NOT p.t_resolution_method 
      OR r.uncertain IS NOT p.n_uncertain 
      OR r.ltr IS NOT p.n_ltr 
      OR r.aggregates IS NOT p.n_aggregates 
      OR r.size_info IS NOT p.n_size_info 
      OR r.morphotype IS NOT p.n_morphotypes 
      OR r.life_stage IS NOT p.n_life_stages 
      OR r.infraspecific_taxa IS NOT p.n_infraspecific_taxa 
      OR r.trophic IS NOT p.n_trophic 
      OR r.local_authority IS NOT p.n_local_authority 
      OR (
        p.taxon_record_notes IS NOT NULL 
        AND r.notes IS NOT p.taxon_record_notes
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'abundance_row_count', 
  (
    SELECT 
      ABS(
        (
          SELECT 
            COUNT(*) 
          FROM 
            abundance
        ) - (
          SELECT 
            COUNT(DISTINCT abundance_id) 
          FROM 
            abundance_full_pre
        )
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'abundance_values', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      abundance_full_pre AS p 
      LEFT JOIN abundance AS r ON r.id = p.abundance_id 
    WHERE 
      r.id IS NULL 
      OR r.sample_id IS NOT p.sample_id 
      OR r.taxon_record_id IS NOT p.taxon_record_id 
      OR r.sampling_protocol_id IS NOT p.sampling_protocol_id 
      OR r.abundance IS NOT p.abundance 
      OR (
        p.abundance_notes IS NOT NULL 
        AND r.notes IS NOT p.abundance_notes
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'source_metadata_values', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      source AS r 
      JOIN curated_source AS c ON c.source = r.source 
      AND COALESCE(c.subsource, '') = COALESCE(r.subsource, '') 
    WHERE 
      r.metadata IS NOT c.metadata 
      OR r.url IS NOT c.url 
      OR r.notes IS NOT c.notes
  ); INSERT INTO migration_assertion 
SELECT 
  'flat_row_count', 
  (
    SELECT 
      ABS(
        (
          SELECT 
            COUNT(*) 
          FROM 
            abundance
        ) - (
          SELECT 
            COUNT(*) 
          FROM 
            v_abundance_flat
        )
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'foreign_keys_before_cleanup', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      pragma_foreign_key_check
  ); 
SELECT 
  'source_rows' AS metric, 
  COUNT(*) AS value 
FROM 
  abundance_full 
UNION ALL 
SELECT 
  'retained_source_rows', 
  COUNT(*) 
FROM 
  abundance_full_pre 
UNION ALL 
SELECT 
  'distinct_abundance_rows', 
  COUNT(*) 
FROM 
  abundance 
UNION ALL 
SELECT 
  'collapsed_exact_abundance_duplicates', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      abundance_full_pre
  ) - (
    SELECT 
      COUNT(*) 
    FROM 
      abundance
  ); 
------------------------------
-- 8. Remove useless tables --
------------------------------
DROP 
  TABLE abundance_full_pre; 
DROP 
  TABLE abundance_full; 
---------------------------------------------------------
-- 9. Verify the release object and database integrity --
---------------------------------------------------------
INSERT INTO migration_assertion 
SELECT 
  'unexpected_tables', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      sqlite_schema 
    WHERE 
      type = 'table' 
      AND name NOT LIKE 'sqlite_%' 
      AND name NOT IN (
        'coordinates', 'event_datetime', 
        'depth', 'taxonomy', 'sampling_protocol', 
        'source', 'event', 'sample', 'taxon_record', 
        'abundance', 'metadata_dictionary'
      )
  ); INSERT INTO migration_assertion 
SELECT 
  'unexpected_views', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      sqlite_schema 
    WHERE 
      type = 'view' 
      AND name <> 'v_abundance_flat'
  ); INSERT INTO migration_assertion 
SELECT 
  'foreign_keys_after_cleanup', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      pragma_foreign_key_check
  ); INSERT INTO migration_assertion 
SELECT 
  'database_integrity', 
  (
    SELECT 
      COUNT(*) 
    FROM 
      pragma_integrity_check 
    WHERE 
      integrity_check <> 'ok'
  ); 
SELECT 
  check_name, 
  failures 
FROM 
  migration_assertion 
ORDER BY 
  check_name; 
COMMIT; 

PRAGMA foreign_key_check; 
PRAGMA integrity_check;
