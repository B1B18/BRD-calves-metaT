############################################################
# 01_permanova.R
#
# PERMANOVA of nasal profiles against clinical status, age group
# and farm of origin, for three data layers:
#
#   Section A  Viral species     normalised abundance; log1p;
#                                Euclidean distance
#   Section B  Bacterial genera  ALL genera in the GTDB genus table
#                                (00_data/all_merged_meta_gtdb_genus.txt);
#                                reads per 6e7 paired reads; log1p;
#                                Euclidean distance
#   Section C  Host transcriptome genes retained by filterByExpr;
#                                TMM logCPM; Euclidean distance
#
#   All three layers use Euclidean distance on a log scale without
#   per-feature standardisation, so each feature contributes in proportion
#   to its variance: abundant, variable taxa (and genes) carry more weight
#   than rare ones, and R2 is the share of total log-scale variance
#   explained. Per-feature z-scoring (equal weight per feature) is kept as
#   an option (zscore_features) for sensitivity analysis.
#
# Samples
#   One common set for all three layers: sequenced nasal swabs ("_N_" in the
#   sample ID) with a
#   recorded clinical status, farm and age group (<25 / >=25 days),
#   present in the viral, bacterial and host matrices. All clinical
#   statuses are included; no case-control restriction. Calves with no
#   reads for the retained features are kept.
#
#   Age group: Age_group in metadata.csv (age in days where it is a single
#   value; otherwise the recorded category if it falls clearly on one side
#   of 25 days; derived upstream). Calves without an age group are excluded.
#
# Features (selected independently of clinical status, across all
# sequenced calves)
#   Viruses : target species detected (normalised abundance > 5, as in
#             02 and 04) in >= 5 calves
#   Bacteria: all genera in the GTDB genus table detected (> 3 reads per
#             6e7 paired reads) in >= 1 sequenced calf. Values at or below
#             the detection threshold are set to 0 before the distance
#             (bact_zero_below_detect), matching taxa_long.csv, where no
#             non-zero value fell at or below the threshold.
#
# Model (all layers)
#   distance ~ Status + Age_category + Farm_Number
#   adonis2, 9,999 permutations, marginal terms (by = "margin"):
#   each term is tested after accounting for the other two. vegan permutes
#   the residuals of the reduced model (the model without the term tested),
#   so farm differences are removed before permuting.
#
# Sensitivity analysis (farm clustering)
#   Same model, with permutations restricted to calves within the same farm
#   (permute::how(blocks = Farm_Number)). Reported for Status and Age only:
#   farm cannot be tested when calves are permuted only within farms.
#   Informative farms (with more than one status / age group) are counted.
#
# Diagnostics
#   Design-matrix rank; homogeneity of multivariate dispersion
#   (betadisper + permutest) for clinical status and age group.
#
# Inputs
#   input_files/metadata.csv        Status, Age_group, Farm_Number, Total_paired_reads
#   input_files/taxa_long.csv       viral normalised abundance
#   input_files/host_counts.csv.gz  RSEM gene counts
#   00_data/all_merged_meta_gtdb_genus.txt  bacterial read counts, genera x samples
#
# Outputs (01_permanova/)
#   log.txt                            full console log + session information
#   Table2_permanova_manuscript.csv    Table 2 as formatted for the manuscript
#   Table2_permanova_summary.csv       Table 2 with Df and exact p-values
#   permanova_<layer>.csv              full adonis2 table per layer
#   permanova_within_farm_permutations.csv  sensitivity: within-farm permutations
#   dispersion_summary.csv             betadisper tests, all layers
#   permanova_common_samples.csv       calves analysed
#   permanova_sample_availability.csv  calves per layer before intersection
#   virus_feature_filter.csv           prevalence of every target virus
#   bacteria_feature_filter.csv        prevalence of every genus
############################################################


############################################################
# 0 — SETUP
############################################################

# Run from the repository root: all paths below are relative to it.

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(edgeR)
})

in_dir    <- "input_files"
bact_file <- file.path("00_data", "all_merged_meta_gtdb_genus.txt")
out_dir   <- "01_permanova"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- parameters ----
permanova_seed <- 42
n_perm         <- 9999     # set to 999 for quick checks
n_cores        <- max(1, parallel::detectCores() - 1)

virus_targets <- c(
  "Bovine adenovirus 3", "Bovine coronavirus", "Bovine nebovirus",
  "Bovine picornavirus", "Bovine rhinitis A", "Bovine rhinitis B",
  "Bovine torovirus", "BVDV 1", "Bovine respiratory and syncytial virus",
  "Bovine parainfluenza 3", "Bovine nidovirus", "Influenza D virus"
)
virus_detect_above <- 5        # normalised abundance (detected, as in 02 and 04)
virus_min_calves   <- 5

bact_reads_scale       <- 6e7  # reads per 6e7 paired reads
bact_detect_above      <- 3
bact_min_calves        <- 1    # all genera detected in at least one sequenced calf
bact_zero_below_detect <- TRUE # set values <= bact_detect_above to 0 before the distance

zscore_features <- FALSE       # TRUE = z-score each microbial feature after log1p (sensitivity only)

# ---- log: console output and messages/warnings go to log.txt ----
# If the script stops with an error, run  sink(); sink(type = "message")
log_con <- file(file.path(out_dir, "log.txt"), open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

options(width = 250, readr.show_progress = FALSE)
cat("01_permanova.R\n")
cat("Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Virus targets:", paste(virus_targets, collapse = ", "), "\n")
cat("Permutations:", n_perm, "| cores:", n_cores, "| seed:", permanova_seed, "\n")
cat("Virus filter: detected (>", virus_detect_above, ") in >=", virus_min_calves, "calves\n")
cat("Bacteria input:", bact_file, "\n")
cat("Bacteria filter: >", bact_detect_above, "reads per", format(bact_reads_scale, scientific = TRUE),
    "paired reads in >=", bact_min_calves, "calves | sub-threshold values zeroed:",
    bact_zero_below_detect, "\n")
cat("Microbial features z-scored before distance:", zscore_features, "\n")


############################################################
# 1 — HELPERS
############################################################

# stop with the message shown on the console AND written to the log
fail <- function(...) {
  msg <- paste0(...)
  cat("\nERROR:", msg, "\n")
  sink(type = "message"); sink()
  stop(msg, call. = FALSE)
}

# taxa x samples matrix from taxa_long.csv
taxa_matrix <- function(taxa, layer) {
  d <- taxa_long %>% filter(type == layer, taxon %in% taxa)
  missing <- setdiff(taxa, d$taxon)
  if (length(missing) > 0) fail("Taxa absent from taxa_long.csv (", layer, "): ", paste(missing, collapse = ", "))
  d %>% dplyr::select(taxon, Sample_id, value) %>%
    pivot_wider(names_from = Sample_id, values_from = value, values_fill = 0) %>%
    column_to_rownames("taxon") %>% as.matrix() %>% .[taxa, , drop = FALSE]
}

prepare_meta <- function(meta) {
  meta %>%
    mutate(Status       = droplevels(factor(Status)),
           Age_category = factor(Age_category, levels = c("<25", ">=25")),
           Farm_Number  = droplevels(factor(Farm_Number)))
}

# status-independent prevalence filter across the sequenced cohort
select_features <- function(mat, detect_above, min_calves, cohort_ids, layer, audit_file) {
  m <- mat[rownames(mat) %in% cohort_ids, , drop = FALSE]
  n_det <- colSums(m > detect_above)
  audit <- tibble(
    feature           = colnames(m),
    n_calves_detected = n_det,
    pct_calves        = round(100 * n_det / nrow(m), 1),
    retained          = n_det >= min_calves
  ) %>%
    arrange(desc(n_calves_detected))
  write_csv(audit, file.path(out_dir, audit_file))
  cat("\n", layer, " features: ", nrow(audit), " evaluated in ", nrow(m),
      " sequenced calves; ", sum(audit$retained), " retained\n", sep = "")
  print(audit %>% filter(retained), n = Inf)
  if (!any(audit$retained)) fail("No ", layer, " features pass the prevalence filter.")
  audit$feature[audit$retained]
}

# log1p (optionally z-scored per feature), Euclidean distance
microbial_distance <- function(mat, meta, zscore = zscore_features) {
  stopifnot(identical(rownames(mat), meta$Sample_id))
  feature_sd <- apply(mat, 2, sd)
  keep <- is.finite(feature_sd) & feature_sd > 0
  if (any(!keep)) cat("Removing zero-variance features (", sum(!keep), "): ",
                      paste(colnames(mat)[!keep], collapse = ", "), "\n", sep = "")
  mat <- mat[, keep, drop = FALSE]
  if (ncol(mat) == 0) fail("No variable microbial features remain.")
  mat_log <- log1p(mat)
  if (zscore) mat_log <- scale(mat_log)
  if (any(!is.finite(mat_log))) fail("Non-finite values after transformation.")

  # share of total variance carried by each feature (what drives the distance)
  v <- apply(mat_log, 2, var)
  share <- sort(v / sum(v), decreasing = TRUE)
  cat("Features used for distance:", ncol(mat_log),
      "| z-scored:", zscore,
      "| top 5 features carry", round(100 * sum(head(share, 5)), 1), "% of total variance\n")
  print(round(100 * head(share, 10), 1))
  dist(mat_log, method = "euclidean")
}

# homogeneity of multivariate dispersion for one grouping
run_dispersion <- function(dist_mat, meta, term, layer) {
  group <- droplevels(factor(meta[[term]]))
  if (nlevels(group) < 2) {
    warning("Dispersion test skipped for ", layer, " / ", term, ": fewer than two groups.")
    return(NULL)
  }
  bd <- betadisper(dist_mat, group)
  set.seed(permanova_seed)
  pt <- permutest(bd, permutations = n_perm, parallel = n_cores)
  cat("\nDispersion —", layer, "/", term, "\n")
  print(pt)
  tab <- pt$tab
  tibble(layer = layer, grouping = term,
         Df = tab["Groups", "Df"], F = tab["Groups", "F"],
         p = tab["Groups", "Pr(>F)"])
}

# marginal PERMANOVA + diagnostics
run_permanova <- function(dist_mat, meta, layer, out_file) {
  
  meta <- prepare_meta(meta)
  
  if (attr(dist_mat, "Size") != nrow(meta)) fail(layer, ": distance size does not match metadata.")
  if (!is.null(labels(dist_mat)) && !identical(labels(dist_mat), meta$Sample_id)) {
    fail(layer, ": distance-matrix sample order does not match metadata.")
  }
  for (v in c("Status", "Age_category", "Farm_Number")) {
    if (nlevels(droplevels(meta[[v]])) < 2) fail(layer, ": ", v, " has fewer than two levels.")
  }
  
  design <- model.matrix(~ Status + Age_category + Farm_Number, data = meta)
  design_rank <- qr(design)$rank
  
  cat("\n====================================\n")
  cat(" PERMANOVA —", layer, "\n")
  cat("====================================\n")
  cat("Samples:", nrow(meta), "| Farms:", nlevels(meta$Farm_Number),
      "| model-matrix rank:", design_rank, "of", ncol(design), "\n")
  print(table(Status = meta$Status, Age = meta$Age_category))
  if (design_rank < ncol(design)) {
    warning(layer, ": model matrix is rank deficient (", design_rank, " of ", ncol(design),
            "). Check confounding among status, age and farm.")
  }
  
  set.seed(permanova_seed)
  res <- adonis2(dist_mat ~ Status + Age_category + Farm_Number,
                 data = meta, permutations = n_perm, by = "margin",
                 parallel = n_cores)
  print(res)
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("term") %>%
    mutate(layer = layer, n_samples = nrow(meta),
           n_farms = nlevels(meta$Farm_Number), .before = 1)
  write_csv(res_df, file.path(out_dir, out_file))
  
  # ---- sensitivity: permutations restricted within farms ----
  informative <- meta %>% group_by(Farm_Number) %>%
    summarise(status = n_distinct(Status) > 1, age = n_distinct(Age_category) > 1, .groups = "drop")
  set.seed(permanova_seed)
  res_w <- adonis2(dist_mat ~ Status + Age_category + Farm_Number,
                   data = meta, permutations = how(blocks = meta$Farm_Number, nperm = n_perm),
                   by = "margin", parallel = n_cores)
  cat("\nSensitivity — permutations restricted within farms (Status and Age only)\n")
  cat("Farms with > 1 clinical status:", sum(informative$status), "of", nrow(informative),
      "| with both age groups:", sum(informative$age), "\n")
  print(res_w)
  within_df <- as.data.frame(res_w) %>%
    rownames_to_column("term") %>%
    filter(term %in% c("Status", "Age_category")) %>%
    transmute(layer = layer, term, Df, R2, F, p_within_farm = `Pr(>F)`,
              p_unrestricted = res_df$`Pr(>F)`[match(term, res_df$term)],
              informative_farms = if_else(term == "Status", sum(informative$status), sum(informative$age)),
              n_farms = nrow(informative))

  dispersion <- bind_rows(
    run_dispersion(dist_mat, meta, "Status", layer),
    run_dispersion(dist_mat, meta, "Age_category", layer)
  )
  
  list(permanova = res_df, dispersion = dispersion, within_farm = within_df)
}


############################################################
# 2 — SAMPLES AND METADATA
############################################################

metadata_seq <- read_csv(file.path(in_dir, "metadata.csv"), col_types = cols(
  Sample_id = col_character(), Status = col_character(), Age_days = col_double(),
  Age_group = col_character(), Farm_Number = col_character(), Total_paired_reads = col_double(),
  has_virus = col_logical(), has_bacteria = col_logical(), has_host = col_logical()))
if (anyDuplicated(metadata_seq$Sample_id)) fail("Duplicate Sample_id values in metadata.csv.")

taxa_long <- read_csv(file.path(in_dir, "taxa_long.csv"), col_types = cols(
  Sample_id = col_character(), taxon = col_character(), taxon_label = col_character(),
  type = col_character(), value = col_double()))

# cohort for the feature filter: all sequenced samples with metadata
sequenced_ids <- metadata_seq$Sample_id
cat("\nSequenced samples with metadata:", length(sequenced_ids), "\n")

metadata_all <- metadata_seq %>% filter(str_detect(Sample_id, "_N_"))
cat("Nasal-labelled samples (_N_ in the sample ID):", nrow(metadata_all), "\n")

cat("\nClinical status:\n")
metadata_all %>% count(Status) %>% print(n = Inf)
cat("\nAge group:\n")
metadata_all %>% count(Age_group) %>% print(n = Inf)

# ---- eligible calves ----
meta_perm <- metadata_all %>%
  filter(!is.na(Age_group), !is.na(Status), !is.na(Farm_Number)) %>%
  mutate(Age_category = Age_group)

bad_age <- setdiff(unique(meta_perm$Age_category), c("<25", ">=25"))
if (length(bad_age)) fail("Unexpected Age_group values: ", paste(bad_age, collapse = ", "))

cat("\nWith status, farm and age group:", nrow(meta_perm), "\n")
cat("Excluded — age group unknown:", sum(is.na(metadata_all$Age_group)),
    "| status missing:", sum(is.na(metadata_all$Status)),
    "| farm missing:", sum(is.na(metadata_all$Farm_Number)), "\n")


############################################################
# 3 — DATA LAYERS
############################################################

# ---- 3A viral species (samples x viruses) ----
virus_mat <- t(taxa_matrix(virus_targets, "virus"))

# ---- 3B bacterial genera (samples x genera), reads per 6e7 paired reads ----
bact_raw <- read_tsv(bact_file, col_types = cols(clade_name = col_character(), .default = col_double()))

non_genus <- bact_raw$clade_name[!str_starts(bact_raw$clade_name, "g__")]
if (length(non_genus) > 0) cat("\nNon-genus rows dropped from", bact_file, ":",
                               paste(non_genus, collapse = ", "), "\n")

bact_raw <- bact_raw %>%
  filter(str_starts(clade_name, "g__")) %>%
  mutate(clade_name = str_remove(clade_name, "^g__"))
if (anyDuplicated(bact_raw$clade_name)) fail("Duplicate genus names in ", bact_file)

bacteria_mat <- bact_raw %>% column_to_rownames("clade_name") %>% as.matrix() %>% t()
if (anyNA(bacteria_mat)) fail("NA values in ", bact_file)
if (any(bacteria_mat < 0)) fail("Negative counts in ", bact_file)
cat("\nGTDB genus table:", ncol(bacteria_mat), "genera x", nrow(bacteria_mat), "samples\n")

extra <- setdiff(rownames(bacteria_mat), metadata_seq$Sample_id)
if (length(extra) > 0) {
  cat("Samples in the genus table but not in metadata.csv (dropped):", paste(extra, collapse = ", "), "\n")
  bacteria_mat <- bacteria_mat[!rownames(bacteria_mat) %in% extra, , drop = FALSE]
}

reads <- metadata_seq$Total_paired_reads[match(rownames(bacteria_mat), metadata_seq$Sample_id)]
if (any(!is.finite(reads) | reads <= 0)) fail("Missing or non-positive Total_paired_reads for bacterial samples.")
bacteria_mat <- sweep(bacteria_mat, 1, reads, "/") * bact_reads_scale

if (bact_zero_below_detect) {
  n_sub <- sum(bacteria_mat > 0 & bacteria_mat <= bact_detect_above)
  cat("Bacterial values > 0 and <=", bact_detect_above, "reads per",
      format(bact_reads_scale, scientific = TRUE), "set to 0:", n_sub,
      "of", sum(bacteria_mat > 0), "non-zero values\n")
  bacteria_mat[bacteria_mat <= bact_detect_above] <- 0
}

# ---- 3C host gene counts (genes x samples) ----
host_mat <- read_csv(file.path(in_dir, "host_counts.csv.gz"),
                     col_types = cols(gene_id = col_character(), .default = col_double())) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()
if (anyNA(host_mat)) fail("NA values in host_counts.csv.gz.")


############################################################
# 4 — FEATURE SELECTION (independent of clinical status)
############################################################

virus_features <- select_features(virus_mat, virus_detect_above, virus_min_calves,
                                  sequenced_ids, "Virus", "virus_feature_filter.csv")

bact_features  <- select_features(bacteria_mat, bact_detect_above, bact_min_calves,
                                  sequenced_ids, "Bacteria", "bacteria_feature_filter.csv")


############################################################
# 5 — ONE COMMON SAMPLE SET
############################################################

sample_availability <- tibble(
  layer           = c("Virus", "Bacteria", "Host expression"),
  eligible_calves = nrow(meta_perm),
  present         = c(sum(meta_perm$Sample_id %in% rownames(virus_mat)),
                      sum(meta_perm$Sample_id %in% rownames(bacteria_mat)),
                      sum(meta_perm$Sample_id %in% colnames(host_mat)))
)
cat("\nEligible calves present in each layer:\n")
print(sample_availability)
write_csv(sample_availability, file.path(out_dir, "permanova_sample_availability.csv"))

for (lyr in list(list("Virus", rownames(virus_mat)),
                 list("Bacteria", rownames(bacteria_mat)),
                 list("Host expression", colnames(host_mat)))) {
  miss <- setdiff(meta_perm$Sample_id, lyr[[2]])
  if (length(miss) > 0) cat("Missing from", lyr[[1]], ":", paste(miss, collapse = ", "), "\n")
}

common_ids <- meta_perm$Sample_id[
  meta_perm$Sample_id %in% rownames(virus_mat) &
    meta_perm$Sample_id %in% rownames(bacteria_mat) &
    meta_perm$Sample_id %in% colnames(host_mat)
]
if (length(common_ids) < 10) fail("Fewer than 10 calves in the common sample set.")

meta_common <- meta_perm %>%
  filter(Sample_id %in% common_ids) %>%
  arrange(match(Sample_id, common_ids))
stopifnot(identical(meta_common$Sample_id, common_ids))

cat("\n====================================\n")
cat(" COMMON SAMPLE SET\n")
cat("====================================\n")
cat("Calves:", nrow(meta_common), "| Farms:", n_distinct(meta_common$Farm_Number), "\n")
print(table(Status = meta_common$Status, Age = meta_common$Age_category))

write_csv(meta_common %>% dplyr::select(Sample_id, Status, Age_category, Farm_Number),
          file.path(out_dir, "permanova_common_samples.csv"))

virus_mat    <- virus_mat[common_ids, virus_features, drop = FALSE]
bacteria_mat <- bacteria_mat[common_ids, bact_features, drop = FALSE]
host_mat     <- host_mat[, common_ids, drop = FALSE]


############################################################
# SECTION A — VIRAL SPECIES
############################################################

res_virus <- run_permanova(microbial_distance(virus_mat, meta_common), meta_common,
                           layer = "Virus", out_file = "permanova_virus.csv")


############################################################
# SECTION B — BACTERIAL GENERA
############################################################

res_bact <- run_permanova(microbial_distance(bacteria_mat, meta_common), meta_common,
                          layer = "Bacteria", out_file = "permanova_bacteria_genus.csv")


############################################################
# SECTION C — HOST TRANSCRIPTOME
############################################################

y_all <- DGEList(counts = host_mat)

# filterByExpr uses group sizes to set the minimum number of samples;
# it does not select genes for association with a group
keep_genes <- filterByExpr(y_all, group = factor(meta_common$Status))
cat("\nHost genes:", nrow(y_all), "| kept after filterByExpr:", sum(keep_genes), "\n")

y_all <- calcNormFactors(y_all[keep_genes, , keep.lib.sizes = FALSE], method = "TMM")
host_logcpm <- cpm(y_all, log = TRUE, prior.count = 2)
if (any(!is.finite(host_logcpm))) fail("Non-finite values in host logCPM matrix.")

res_host <- run_permanova(dist(t(host_logcpm), method = "euclidean"), meta_common,
                          layer = "Host expression", out_file = "permanova_host_expression.csv")


############################################################
# 6 — TABLE 2 AND DISPERSION SUMMARY
############################################################

perm_all <- bind_rows(res_bact$permanova, res_virus$permanova, res_host$permanova) %>%
  filter(term %in% c("Age_category", "Farm_Number", "Status")) %>%
  mutate(term = dplyr::recode(term, Age_category = "Age", Farm_Number = "Farm"))

# exact values
table2_summary <- perm_all %>%
  dplyr::select(layer, n_samples, n_farms, term, Df, R2, p = `Pr(>F)`) %>%
  mutate(R2 = round(R2, 4)) %>%
  pivot_wider(names_from = term, values_from = c(Df, R2, p), names_glue = "{term}_{.value}")
write_csv(table2_summary, file.path(out_dir, "Table2_permanova_summary.csv"))

# manuscript format
fmt_p <- function(p) case_when(is.na(p) ~ NA_character_,
                               p < 0.001 ~ "<0.001",
                               TRUE ~ formatC(p, format = "f", digits = 3))

table2_manuscript <- perm_all %>%
  mutate(Category = dplyr::recode(layer, Bacteria = "Bacteria genus"),
         R2 = formatC(R2, format = "f", digits = 3),
         p  = fmt_p(`Pr(>F)`)) %>%
  dplyr::select(Category, n_samples, n_farms, term, R2, p) %>%
  pivot_wider(names_from = term, values_from = c(R2, p), names_glue = "{term} {.value}") %>%
  dplyr::select(Category, n_samples, n_farms,
                `Age R2`, `Age p`, `Farm R2`, `Farm p`, `Status R2`, `Status p`) %>%
  mutate(Category = factor(Category, levels = c("Bacteria genus", "Virus", "Host expression"))) %>%
  arrange(Category)
write_csv(table2_manuscript, file.path(out_dir, "Table2_permanova_manuscript.csv"))

dispersion_summary <- bind_rows(res_bact$dispersion, res_virus$dispersion, res_host$dispersion) %>%
  mutate(F = round(F, 3))
write_csv(dispersion_summary, file.path(out_dir, "dispersion_summary.csv"))

within_farm <- bind_rows(res_bact$within_farm, res_virus$within_farm, res_host$within_farm) %>% as_tibble() %>%
  mutate(across(c(R2, F), ~ round(.x, 4)))
write_csv(within_farm, file.path(out_dir, "permanova_within_farm_permutations.csv"))

cat("\n====================================\n")
cat(" TABLE 2 — PERMANOVA (marginal terms)\n")
cat("====================================\n")
print(table2_manuscript, width = Inf)

cat("\n====================================\n")
cat(" SENSITIVITY — PERMUTATIONS WITHIN FARMS (Status, Age)\n")
cat("====================================\n")
print(within_farm, width = Inf)

cat("\n====================================\n")
cat(" DISPERSION (betadisper)\n")
cat("====================================\n")
print(dispersion_summary, width = Inf)


############################################################
# 7 — SESSION INFORMATION AND CLOSE LOG
############################################################

cat("\n====================================\n")
cat(" SESSION INFORMATION\n")
cat("====================================\n")
print(sessionInfo())

cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Outputs in:", out_dir, "\n")

sink(type = "message")
sink()
close(log_con)


############################################################
# 8 — DISPERSION DETAIL (printed to the console only)
# Mean distance to centroid and pairwise dispersion tests by clinical status.
############################################################

disp_detail <- function(d, meta, layer) {
  bd <- betadisper(d, meta$Status)
  cat("\n", layer, "— mean distance to centroid by status\n")
  print(round(tapply(bd$distances, meta$Status, mean), 2))
  set.seed(permanova_seed)
  print(permutest(bd, pairwise = TRUE, permutations = n_perm)$pairwise$permuted)
}
mc <- prepare_meta(meta_common)
disp_detail(microbial_distance(bacteria_mat, mc), mc, "Bacteria")
disp_detail(microbial_distance(virus_mat, mc),    mc, "Virus")

# excluding calves with no reads for the retained taxa
keep <- rowSums(bacteria_mat) > 0
disp_detail(microbial_distance(bacteria_mat[keep, ], mc[keep, ]), mc[keep, ], "Bacteria (reads > 0)")