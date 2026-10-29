# Input files

Processed data used as input by the scripts in `scripts/`. The raw sequencing reads are in the NCBI Sequence Read Archive under BioProject [PRJNA1191164](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1191164); the processing steps that produced these tables are described in the main README and in the Methods of the manuscript.

Sample identifiers (`Sample_id`) are the same in every file and match the BioSample records of PRJNA1191164.

| File | Content | Used by |
|---|---|---|
| `metadata.csv` | One row per sequenced sample. Columns: `Sample_id`; `Status` (clinical status: healthy, respiratory, unspecific, scour, respiratory_scours); `Age_days` (age in days, where recorded as a single value; otherwise empty); `Age_group` (<25 or >=25 days; empty where it could not be assigned); `Farm_Number` (anonymised farm code); `Total_paired_reads` (paired reads sequenced); `has_virus`, `has_bacteria`, `has_host` (TRUE/FALSE: sample has viral, bacterial and host expression data). | all |
| `taxa_long.csv` | Microbial abundance in long format, one row per sample and taxon. Columns: `Sample_id`; `taxon` (viral species or bacterial genus, GTDB r220, without the `g__` prefix); `taxon_label` (abbreviation used in the figures); `type` (`virus` or `bacteria`); `value` (viruses: normalised abundance; bacteria: read counts). Sample–taxon combinations with no reads may be absent; the scripts treat them as zero. | 01, 02, 03a, 04a, 05b |
| `host_counts.csv.gz` | Host gene-level read counts (nf-core/rnaseq, STAR–RSEM, ARS-UCD1.2). Column `gene_id` (RefSeq gene identifier, `gene-` prefix) followed by one column per sample. | 01, 03a, 04a |
| `all_merged_meta_gtdb_genus.txt` | MetaPhlAn 4 read counts for all bacterial genera (GTDB r220), tab-separated. Column `clade_name` (`g__` prefix) followed by one column per sample. | 01 |
| `input_annotations.rds` | Cache of the Ensembl BioMart query (release 115, ARS-UCD1.3): human orthologues and Entrez Gene IDs for the annotated host genes. An R list with elements `basic`, `entrez`, `genes`, `method`, `release`, `host` and `date`. With this file, `05a_gene_annotation.R` reproduces the published annotation without querying Ensembl. | 05a |
