# Phyto-MAP data-processing scripts

<!-- badges: start -->
[![Version: 1.0.0](https://img.shields.io/badge/version-1.0.0-blue.svg)](https://github.com/RayonVina/phytomap-code/releases/tag/v1.0.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23219719.svg)](https://doi.org/10.5281/zenodo.23219719)
[![Last commit](https://img.shields.io/github/last-commit/RayonVina/phytomap-code?style=flat-square&color=b4befe)](https://github.com/RayonVina/phytomap-code/commits/main/)
<!-- badges: end -->

## Purpose and scope

This collection contains edited and standardised versions of scripts used during the development of Phyto-MAP, a curated database of marine phytoplankton abundance records. It supports inspection of the retained programmatic components of source-data processing, taxonomic curation and relational database construction.

The collection is not a self-contained pipeline for regenerating the published database. Original input files, the development database and some intermediate resources are not included. Some curation and correction steps were performed directly in the working database without separate scripts being retained.

The archived scripts therefore document the retained programmatic components, not a complete executable history of every change made to the data.

The database release currently being prepared is Phyto-MAP v2.1.0. The database and its release-specific preparation materials are managed separately from this processing-code collection. Code and dataset version numbers are independent.

## Author and attribution

The author responsible for this edited and standardised script collection is Fernando Rayón Viña:

- GitHub: [RayonVina](https://github.com/RayonVina)
- ORCID: [0000-0002-1622-2180](https://orcid.org/0000-0002-1622-2180)

Credits for earlier script versions and contributions recorded in individual file headers are retained. Collection-level authorship does not replace those existing attributions or any applicable third-party notices.

## Files and processing stages

The collection is organised as follows. During manuscript preparation these files are held in `code/`; in the independent code repository this collection will form the repository root.

```text
Phyto-MAP processing-code collection/
  README.md
  LICENSE
  01 Data extraction/
    01 <source>.R
    ...
  02 master_taxonomy.R
  03 phytomap_relational.sql
```

The stage-01 directory groups the source-specific scripts for readability. It does not contain the original source data or recreate their original working directories.

| Stage | Files | Purpose |
|---|---|---|
| 01 | `01 Data extraction/01 <source>.R` | Source-specific import, curation and preparation of records for the working database. |
| 02 | `02 master_taxonomy.R` | Cleaning taxonomic designations, resolving names and exporting the working taxonomic mapping. |
| 03 | `03 phytomap_relational.sql` | Creating and populating the relational structure from the assembled working data, with checks encoded in the SQL script. |

The stage-01 scripts are independent of one another, although they can write to a shared working database. The main workflow proceeds through stage 01, stage 02 and stage 03, with intermediate assembly and direct database operations between these stages. This ordering does not mean that running the archived files consecutively will reconstruct the complete database. Some stage-01 scripts also contain taxonomic processing.

### Stage 01: source-specific processing

Depending on the source, scripts reshape abundance tables, join sampling metadata, convert abundance units, prepare temporal and spatial information, handle source-specific missing-value conventions and assign source or subsource information. Their outputs are imported into `abundance_raw`, a staging table in the development database.

Transformations differ between sources. An operation in one script must not be interpreted as a rule applied to all contributing datasets. Abundance is standardised to cells L^-1. Retained zeros denote source-specific non-detection under the relevant sampling and analytical protocol, not universal ecological absence. Interpreting source-specific missing values as zeros requires support from the source documentation.

### Stage 02: taxonomic processing

The taxonomic script uses `phytaxr`, working dictionaries and external taxonomic services. It includes automated processing and interactive review of unresolved candidate matches. Existing dictionary contents and curator decisions are inputs to the process, not information that can be recovered from the script alone.

`unique_taxon_clean` is the R object used for the processed taxonomic records. `master_taxonomy` is the working database table to which the resulting mapping is exported. Historical working dictionary names, including `tax_sinom` and `tax_dict`, refer to development resources rather than published dataset tables.

### Intermediate assembly and direct curation

The `master_taxonomy` table was incorporated into the development database, and its information was added to `abundance_raw` through a join. Temporal values ingested as timestamps were subsequently separated into the corresponding date and time columns. Not all operations involved in these steps were retained as scripts; the exact join command and a complete execution history are not supplied here.

Data were subsequently assembled in `abundance_full`, which served as the input to relational database construction. The archived scripts do not provide a complete executable implementation of every operation connecting `abundance_raw`, `master_taxonomy` and `abundance_full`. Missing commands or logs have not been retrospectively invented.

`abundance_raw` is a non-final staging table, not necessarily a SQL TEMP table. Its contents can already include source-specific transformations; its name does not imply that every value is an untouched original observation. Neither `abundance_raw` nor `abundance_full` exists in the final published relational model.

### Stage 03: relational construction

The SQL script creates and populates the relational structure from an existing `abundance_full` table. It includes optional imports for author-curated inputs and checks for counts, reconstructed values, relationships and database integrity. It removes intermediate tables as part of its workflow.

The stage-03 input differs from the initial ingestion output. In particular, the script expects separate date fields and an eight-character time string in the working field named `datetime`. This is not the full timestamp representation generated by some stage-01 scripts. Intermediate preparation is therefore necessary before this SQL script can be run.

Checks encoded in the script are distinct from evidence that they were successfully executed for a particular data release. Release-specific validation evidence belongs with the corresponding dataset documentation.

## Original organisation of inputs

During source processing, scripts were associated with source-specific working directories containing the input data, metadata and supporting documentation for that source. The scripts are collected together here for inspection, but those original directory contents are not included.

The following schematic illustrates this organisation; it is not an exact inventory or a layout shared by every source:

```text
original-working-project/
  <source-directory>/
    <source-processing-script>.R
    <input-data-files>
    <associated-metadata>
    <supporting-documentation>
  <another-source-directory>/
    ...
```

Path handling varies. Some scripts locate inputs relative to the script using `this.path`; others locate a named source directory using `here`; some retain explicit development paths. Individual scripts may read local inputs, access online resources or combine both approaches.

Readers wishing to execute a script must inspect its input-loading block and supply the expected files, directory context and working database. Collecting scripts into one directory does not reproduce their original execution context. Locally prepared input files are not guaranteed to be identical to files currently downloadable from a source provider.

Source materials remain subject to their own access conditions, licences and attribution requirements. Their exclusion from this collection does not imply that every input is restricted or unavailable. Source provenance and bibliographic information should be consulted in the documentation accompanying the dataset.

## Software and environment

The following environment was reported during preparation of this collection:

| Component | Reported value |
|---|---|
| R | 4.6.1 (2026-06-24) |
| Platform | `x86_64-pc-linux-gnu` |
| Operating system | Ubuntu 22.04.3 LTS |
| Locale | `C.UTF-8`, except `LC_NUMERIC=C` |
| System time zone | `America/Halifax` |
| Numerical libraries | OpenBLAS; LAPACK 3.10.0 |
| phytaxr | Installed package version 0.3.2 |

This is a reported environment snapshot, not proof of a single environment used throughout development. The system time zone is not evidence that source sampling times were recorded or harmonised in that time zone.

The supplied session report lists `renv` 1.3.0 among loaded namespaces. It does not establish that a complete dependency lockfile was used or preserved.

Dependencies are declared within individual scripts. Those visible in the inspected code include `tidyverse`, `readxl`, `here`, `this.path`, `lubridate`, `janitor`, `DBI`, `RSQLite`, `furrr`, `remotes` and `phytaxr`; this list is not exhaustive. A complete historical package-version inventory is unavailable, and versions that cannot be recovered are not assigned retrospectively.

`phytaxr` is developed at [RayonVina/phytaxr](https://github.com/RayonVina/phytaxr). The then-current version was used during development rather than one fixed revision throughout the project. The reported installed version is 0.3.2; `RemoteRef` and `RemoteSha` were not available in its package metadata. Consequently, neither a historical Git commit nor a uniform revision for all processing steps can be established from that report. The retained taxonomic script contains an installation block targeting `main`; running it can replace the local package with a different revision.

The relational SQL script contains SQLite command-line directives, including `.bail`, and is intended for the SQLite command-line shell. The historical SQLite version has not been established.

## Execution limitations and safety

These scripts modify a working database. Do not run them against the published database or the only copy of development data.

Several stage-01 scripts and the stage-02 connections use `PHYTOMAP_DB` to identify the working SQLite database. Other path logic may remain in individual scripts and must be checked before execution. Expected staging and dictionary tables must already exist where required.

Scripts may clear the R workspace, change the working directory, install software, access external services, append observations, replace working tables or remove intermediate tables. Repeated execution is not guaranteed to be idempotent. Independent stage-01 scripts are not presented as a concurrent-write workflow.

The collection has been reviewed and edited by its author. That review is not a claim that a complete end-to-end regeneration test was performed. Input availability, manual decisions and unscripted operations limit executable reconstruction of the historical workflow.

## Licence and citation

The code collection is distributed under the MIT License; see `LICENSE`. Existing third-party notices and permissions must be preserved where
applicable. This licence does not apply to the separately distributed database or excluded source data.

The latest archived release is available on Zenodo:
https://doi.org/10.5281/zenodo.23219719

This concept DOI represents all archived versions and always resolves to the latest archived release. Citation metadata are provided in `CITATION.cff`.

The development repository is available at:
https://github.com/RayonVina/phytomap-code

The Phyto-MAP dataset is distributed separately and should be cited independently.
