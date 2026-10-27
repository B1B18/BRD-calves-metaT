############################################################
# 03a_host_DE.R
#
# Host differential expression (limma-voom) in nasal swab
# metatranscriptomes by clinical status and taxon detection.
#
# Pipeline position
#   03a_host_DE.R                 <- this script
#   04a (networks)                reads the full DE results (overlap with the networks)
#   05a_gene_annotation.R         reads 03a outputs
#   05b_figures.R                 reads 03a and 05a outputs
#
# Calves
#   Sequenced calves with clinical status healthy or respiratory, a
#   single-value age in days and a farm number. Calves with other statuses
#   (unspecific, scour, respiratory_scours) are not included.
#
# Taxon detection (per taxon, across all sequenced samples with viral and
# bacterial data, independent of clinical status)
#   abundance == 0                    -> negative
#   0 < abundance < Q25 of non-zero   -> excluded from that taxon's models
#   abundance >= Q25 of non-zero      -> positive
#   (taxa with < q25_min_nonzero non-zero samples are not classified)
#   Viruses : normalised abundance
#   Bacteria: reads per bact_reads_scale paired reads (the classes do not
#             depend on the scaling constant)
#   A taxon is modelled only if >= min_positive_taxon calves in the DE set
#   are positive.
#
# Gene filtering (separately for each model, so genes expressed only in a
#   subset of calves, e.g. taxon-positive calves, are retained)
#   kept if CPM >= min_count_filter / median library size (in millions) in
#   >= min_prop_filter of the calves in the smallest group compared by that
#   model, and >= min_total_count reads in total (filterByExpr logic, with
#   the minimum sample size set to 70% of the smallest compared group);
#   TMM normalisation after filtering, within each model.
#
# Models (limma-voom; farm as blocking factor via duplicateCorrelation, voom
#   and duplicateCorrelation run twice as in the limma user guide; eBayes;
#   BH-adjusted FDR across all genes tested in that model, within each test)
#
#   respiratory_vs_healthy                     ~ 0 + Status + Age_days
#       respiratory vs healthy calves, taxon detection ignored
#       (filter groups: healthy, respiratory)
#   positive_vs_negative                       ~ Status + detection + Age_days
#       taxon-positive vs taxon-negative calves, adjusted for clinical
#       status (detection coefficient; healthy and respiratory calves
#       weighted by the information they carry); the detection effect is
#       assumed to be the same in healthy and respiratory calves
#       (filter groups: taxon-positive, taxon-negative)
#
#   Reporting: DE = FDR < fdr_cutoff and |log2FC| > lfc_cutoff;
#   high fold change = FDR < fdr_cutoff and |log2FC| > lfc_cutoff_high.
#
# Inputs (input_files/)
#   metadata.csv        Status, Age_days, Farm_Number, Total_paired_reads
#   taxa_long.csv       viral normalised abundance, bacterial read counts
#   host_counts.csv.gz  RSEM gene counts
#
# Outputs (03a_host_DE/)
#   log.txt                                          parameters, checks, numbers
#                                                    for the Results, session info
#   03a_objects.rds                                  hand-off to 05a, 05b
#   host_DE_all_results.csv   all genes, all tests
#   host_DE_significant_genes.csv  FDR < fdr_cutoff, |log2FC| > lfc_cutoff
#   host_DE_significant_genes_high_FC.csv            FDR < fdr_cutoff, |log2FC| > lfc_cutoff_high
#   DE_summary_by_taxon_comparison.csv               DE gene counts per test
#   pathogen_DE_sample_counts_q25.csv                detection classes
#   host_expression_sample_summary.csv               every test, run or skipped
#   QC_MDS_status.png, QC_MDS_farm.png, QC_voom_mean_variance.png
############################################################


############################################################
# 0 — PARAMETERS (all user-adjustable settings; see README)
############################################################

# Run from the repository root: all paths below are relative to it.

in_dir  <- "input_files"
out_dir <- "03a_host_DE"

# taxa analysed (taxon names in taxa_long.csv; genera with "g__")
virus_targets <- c("Bovine adenovirus 3", "Bovine coronavirus", "Bovine nebovirus", "Bovine picornavirus",
                   "Bovine rhinitis A", "Bovine rhinitis B", "Bovine torovirus", "BVDV 1",
                   "Influenza D virus")
bacteria_targets <- c("g__Prevotella", "g__Alloprevotella", "g__Moraxella", "g__Moraxella_A", "g__Mannheimia",
                      "g__Mycoplasmopsis_A", "g__Mesomycoplasma", "g__Pasteurella",
                      "g__Streptococcus", "g__Burkholderia", "g__Psychrobacter")

# detection classes
bact_reads_scale <- 6e7    # bacterial reads per this many paired reads (does not change classes)
q25_min_nonzero  <- 5      # taxa with fewer non-zero samples are not classified

# calves and model eligibility
status_levels_de   <- c("healthy", "respiratory")   # statuses in the DE models (reference first)
min_positive_taxon <- 10   # positive calves needed to model a taxon

# gene filtering, per model (filterByExpr defaults except the minimum sample size)
min_count_filter <- 10     # reads; converted to CPM using the median library size of the model's calves
min_prop_filter  <- 0.7    # gene expressed in >= this proportion of the smallest compared group
min_total_count  <- 15     # total reads across the model's calves

# DE thresholds (read by 05a and 05b from 03a_objects.rds)
fdr_cutoff      <- 0.05
lfc_cutoff      <- 1       # DE: counts of up- and down-regulated genes
lfc_cutoff_high <- 2       # high fold change: genes described individually

# QC
mds_top_genes <- 500       # genes used for the MDS plots


############################################################
# 1 — SETUP, LOG AND HELPERS
############################################################

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(tidyverse)
})

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_con <- file(file.path(out_dir, "log.txt"), open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
options(width = 250, max.print = 1e6, readr.show_progress = FALSE)

fail <- function(...) {
  msg <- paste0(...)
  cat("\nERROR:", msg, "\n")
  sink(type = "message"); sink()
  stop(msg, call. = FALSE)
}

section <- function(title) {
  options(width = 250, max.print = 1e6)   # RStudio resets the width when the pane is resized
  cat("\n\n============================================================\n")
  cat(" ", title, "\n")
  cat("============================================================\n")
}

log_parameters <- function(names) {
  section("PARAMETERS")
  for (nm in names) {
    val <- get(nm)
    cat(sprintf("%-22s = %s\n", nm,
                if (length(val) == 0) "(none)" else if (is.numeric(val)) paste(format(val, scientific = FALSE, trim = TRUE), collapse = ", ")
                else paste(val, collapse = ", ")))
  }
}

cat("03a_host_DE.R\nRun:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

params_03a <- c("virus_targets", "bacteria_targets", "bact_reads_scale", "q25_min_nonzero",
                "status_levels_de", "min_positive_taxon",
                "min_count_filter", "min_prop_filter", "min_total_count", "fdr_cutoff",
                "lfc_cutoff", "lfc_cutoff_high", "mds_top_genes")
log_parameters(params_03a)

# taxa x samples matrix from taxa_long.csv (bacterial genera given with "g__")
taxa_matrix <- function(taxa, layer) {
  d <- taxa_long %>% filter(type == layer, taxon %in% str_remove(taxa, "^g__"))
  missing <- setdiff(str_remove(taxa, "^g__"), d$taxon)
  if (length(missing) > 0) fail("Taxa absent from taxa_long.csv (", layer, "): ", paste(missing, collapse = ", "))
  m <- d %>% dplyr::select(taxon, Sample_id, value) %>%
    pivot_wider(names_from = Sample_id, values_from = value, values_fill = 0) %>%
    column_to_rownames("taxon") %>% as.matrix()
  m <- m[str_remove(taxa, "^g__"), , drop = FALSE]
  rownames(m) <- taxa
  m
}

# names used in the DE output files (pathogen_clean): viral abbreviation from
# taxa_long.csv (taxon_label); genera keep "g__"
clean_pathogen_names <- function(x) {
  lab <- virus_labels[str_trim(x)]
  unname(if_else(is.na(lab), str_replace_all(str_trim(x), " ", "_"), lab))
}

# names used in the figures (viral abbreviations, bacterial genera)
taxon_label <- function(x) str_remove(clean_pathogen_names(x), "^g__")

classify_pathogen_q25 <- function(x, min_nonzero = q25_min_nonzero) {
  nonzero <- x[x > 0 & !is.na(x)]
  if (length(nonzero) < min_nonzero) return(rep(NA_character_, length(x)))
  q25 <- as.numeric(quantile(nonzero, 0.25, na.rm = TRUE, names = FALSE))
  case_when(is.na(x)        ~ NA_character_,
            x == 0          ~ "negative",
            x > 0 & x < q25 ~ "exclude",
            x >= q25        ~ "positive")
}


############################################################
# 2 — SAMPLES AND METADATA
############################################################

section("SAMPLES AND METADATA")

meta_seq <- read_csv(file.path(in_dir, "metadata.csv"), col_types = cols(
  Sample_id = col_character(), Status = col_character(), Age_days = col_double(),
  Age_group = col_character(), Farm_Number = col_character(), Total_paired_reads = col_double(),
  has_virus = col_logical(), has_bacteria = col_logical(), has_host = col_logical()))
if (anyDuplicated(meta_seq$Sample_id)) fail("Duplicate Sample_id values in metadata.csv.")
taxa_long <- read_csv(file.path(in_dir, "taxa_long.csv"), col_types = cols(
  Sample_id = col_character(), taxon = col_character(), taxon_label = col_character(),
  type = col_character(), value = col_double()))
n_dup <- taxa_long %>% count(type, taxon, Sample_id) %>% filter(n > 1) %>% nrow()
if (n_dup > 0) fail("Duplicate type/taxon/Sample_id rows in taxa_long.csv: ", n_dup)
virus_labels <- taxa_long %>% filter(type == "virus") %>% distinct(taxon, taxon_label) %>%
  { setNames(.$taxon_label, .$taxon) }

cat("Sequenced samples with metadata:", nrow(meta_seq), "\n")
cat("\nClinical status:\n")
meta_seq %>% count(Status) %>% print(n = Inf)
cat("\nSamples without a single-value age in days:", sum(is.na(meta_seq$Age_days)), "\n")

n_no_farm <- meta_seq %>%
  filter(!is.na(Age_days), Status %in% status_levels_de,
         is.na(Farm_Number) | trimws(as.character(Farm_Number)) == "") %>% nrow()
cat("\nCalves excluded from DE for missing Farm_Number:", n_no_farm, "\n")

metadata_age <- meta_seq %>%
  filter(!is.na(Age_days), Status %in% status_levels_de,
         !is.na(Farm_Number), trimws(as.character(Farm_Number)) != "") %>%
  mutate(Age_category = if_else(Age_days < 25, "<25", ">=25"))

cat("\nCalves eligible for DE (", paste(status_levels_de, collapse = "/"),
    ", age in days, farm recorded):", nrow(metadata_age),
    "| farms:", n_distinct(metadata_age$Farm_Number), "\n")
print(metadata_age %>% count(Status))


############################################################
# 3 — TAXON DETECTION CLASSES (Q25 rule)
############################################################

section("TAXON DETECTION CLASSES")

virus_mat <- taxa_matrix(virus_targets, "virus")
bact_mat  <- taxa_matrix(bacteria_targets, "bacteria")   # read counts

# ---- one sample set: samples flagged in metadata as having viral and bacterial data ----
# taxa_long.csv may omit zero rows, so a calf with data but no target taxon detected
# has no column in virus_mat / bact_mat; it is added with zeros (taxon-negative)
pathogen_samples <- meta_seq$Sample_id[coalesce(meta_seq$has_virus, FALSE) & coalesce(meta_seq$has_bacteria, FALSE)]
if (length(pathogen_samples) == 0) fail("No samples flagged has_virus and has_bacteria in metadata.csv.")

add_missing <- function(m) {
  miss <- setdiff(pathogen_samples, colnames(m))
  if (length(miss) > 0) m <- cbind(m, matrix(0, nrow(m), length(miss), dimnames = list(rownames(m), miss)))
  m[, pathogen_samples, drop = FALSE]
}
n_no_virus_row <- length(setdiff(pathogen_samples, colnames(virus_mat)))
n_no_bact_row  <- length(setdiff(pathogen_samples, colnames(bact_mat)))
n_unflagged    <- length(setdiff(union(colnames(virus_mat), colnames(bact_mat)), pathogen_samples))
virus_mat <- add_missing(virus_mat)
bact_mat  <- add_missing(bact_mat)

cat("Samples with viral and bacterial data (metadata flags):", length(pathogen_samples),
    "| without viral data:", sum(!coalesce(meta_seq$has_virus, FALSE)),
    "| without bacterial data:", sum(!coalesce(meta_seq$has_bacteria, FALSE)), "\n")
cat("  no row for any target virus in taxa_long.csv (added as zero):", n_no_virus_row,
    "| no row for any target bacterium (added as zero):", n_no_bact_row, "\n")
cat("  samples in taxa_long.csv for these taxa but not flagged (not used):", n_unflagged, "\n")


reads <- meta_seq$Total_paired_reads[match(pathogen_samples, meta_seq$Sample_id)]
if (any(!is.finite(reads) | reads <= 0)) fail("Missing or non-positive Total_paired_reads.")
bact_norm <- sweep(bact_mat[, pathogen_samples, drop = FALSE], 2, reads, "/") * bact_reads_scale

pathogen_mat <- rbind(virus_mat[, pathogen_samples, drop = FALSE], bact_norm)
stopifnot(identical(colnames(pathogen_mat), pathogen_samples))

pathogen_classified <- pathogen_mat %>%
  as.data.frame() %>%
  rownames_to_column("pathogen") %>%
  pivot_longer(-pathogen, names_to = "Sample_id", values_to = "abundance") %>%
  left_join(meta_seq %>% dplyr::select(Sample_id, Status, Age_days, Farm_Number), by = "Sample_id") %>%
  group_by(pathogen) %>%
  mutate(class = classify_pathogen_q25(abundance),
         pathogen_clean = clean_pathogen_names(pathogen),
         taxon = taxon_label(pathogen)) %>%
  ungroup()

summary_table <- pathogen_classified %>%
  filter(!is.na(class)) %>%
  group_by(pathogen, pathogen_clean, Status, class) %>%
  summarise(n_samples = n(), .groups = "drop") %>%
  pivot_wider(names_from = class, values_from = n_samples, values_fill = 0) %>%
  arrange(pathogen_clean, Status)
write_csv(summary_table, file.path(out_dir, "pathogen_DE_sample_counts_q25.csv"))

cat("\nDetection classes per taxon and status (all sequenced samples):\n")
print(as.data.frame(summary_table), row.names = FALSE)


############################################################
# 4 — HOST COUNTS: filtering, TMM normalisation
############################################################

section("HOST COUNTS")

host_mat_full <- read_csv(file.path(in_dir, "host_counts.csv.gz"),
                          col_types = cols(gene_id = col_character(), .default = col_double())) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()
if (anyNA(host_mat_full)) fail("NA values in host_counts.csv.gz.")

common_samples <- Reduce(intersect, list(colnames(host_mat_full), pathogen_samples, metadata_age$Sample_id))
if (length(common_samples) < 10) fail("Too few common samples after alignment.")

host_mat <- host_mat_full[, common_samples, drop = FALSE]
meta_host <- metadata_age %>%
  filter(Sample_id %in% common_samples) %>%
  arrange(match(Sample_id, common_samples))
stopifnot(identical(meta_host$Sample_id, colnames(host_mat)))

cat("Calves in the DE models:", ncol(host_mat), "| farms:", n_distinct(meta_host$Farm_Number), "\n")
print(table(Status = meta_host$Status, Age = meta_host$Age_category))

y_raw <- DGEList(counts = host_mat)   # all genes, DE calves; filtered per model

# filter genes for one model and TMM-normalise
#   samples: calves in the model; group_sizes: sizes of the groups the model compares
prepare_y <- function(samples, group_sizes) {
  y <- y_raw[, samples]
  n_min   <- ceiling(min_prop_filter * min(group_sizes))
  cpm_cut <- min_count_filter / median(y$samples$lib.size) * 1e6
  keep <- rowSums(cpm(y) >= cpm_cut) >= n_min & rowSums(y$counts) >= min_total_count
  if (sum(keep) < 100) stop("too few genes retained after filtering (", sum(keep), ")")
  list(y = calcNormFactors(y[keep, , keep.lib.sizes = FALSE]),
       n_genes = sum(keep), n_min = n_min, cpm_cut = cpm_cut)
}

# disease-model gene set: used for QC and handed to 05a/05b
prep_all <- prepare_y(meta_host$Sample_id, table(meta_host$Status))
y_all <- prep_all$y
cat("\nHost genes:", nrow(y_raw), "| kept for the disease model:", prep_all$n_genes,
    "(CPM >=", round(prep_all$cpm_cut, 2), "in >=", prep_all$n_min, "calves)\n")

# expression for the figures (05b): TMM log2 CPM across ALL sequenced calves with
# host, taxon and metadata records (all statuses, including unknown age)
display_samples <- Reduce(intersect, list(colnames(host_mat_full), pathogen_samples, meta_seq$Sample_id))
y_disp <- DGEList(counts = host_mat_full[, display_samples, drop = FALSE])
y_disp <- calcNormFactors(y_disp[rowSums(y_disp$counts) > 0, , keep.lib.sizes = FALSE])
logcpm_disp <- cpm(y_disp, log = TRUE, prior.count = 2)
cat("Calves available for the figures (all statuses):", ncol(logcpm_disp), "\n")


############################################################
# 5 — QC: library sizes, TMM, MDS, mean-variance; farm balance
############################################################

section("QC: LIBRARY SIZES, NORMALISATION, MDS")
qc <- meta_host %>%
  mutate(lib_size = y_all$samples$lib.size[match(Sample_id, colnames(y_all))],
         norm_factor = y_all$samples$norm.factors[match(Sample_id, colnames(y_all))])
cat("Library size (reads assigned to filtered genes) and TMM factor by status:\n")
qc %>% group_by(Status) %>%
  summarise(calves = n(), lib_median = round(median(lib_size)), lib_min = round(min(lib_size)),
            lib_max = round(max(lib_size)), tmm_median = round(median(norm_factor), 2),
            tmm_min = round(min(norm_factor), 2), tmm_max = round(max(norm_factor), 2), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

logcpm_qc <- cpm(y_all, log = TRUE, prior.count = 2)
mds <- plotMDS(logcpm_qc, top = mds_top_genes, plot = FALSE)
mds_df <- tibble(Sample_id = colnames(logcpm_qc), dim1 = mds$x, dim2 = mds$y) %>%
  left_join(meta_host %>% dplyr::select(Sample_id, Status, Age_days, Farm_Number), by = "Sample_id")
p_mds_status <- ggplot(mds_df, aes(dim1, dim2, colour = Status)) +
  geom_point(alpha = 0.8) +
  scale_colour_manual(values = c(healthy = "#2CA02C", respiratory = "#1F77B4")) +
  labs(x = paste0("Dim 1 (", round(100 * mds$var.explained[1], 1), "%)"),
       y = paste0("Dim 2 (", round(100 * mds$var.explained[2], 1), "%)"),
       title = paste0("MDS of host expression (top ", mds_top_genes, " genes), by clinical status")) +
  theme_bw()
p_mds_farm <- ggplot(mds_df, aes(dim1, dim2, colour = factor(Farm_Number))) +
  geom_point(alpha = 0.8) +
  labs(x = "Dim 1", y = "Dim 2", colour = "Farm", title = "MDS of host expression, by farm") +
  theme_bw() + theme(legend.position = "none")
ggsave(file.path(out_dir, "QC_MDS_status.png"), p_mds_status, width = 7, height = 6, dpi = 300)
ggsave(file.path(out_dir, "QC_MDS_farm.png"), p_mds_farm, width = 7, height = 6, dpi = 300)
png(file.path(out_dir, "QC_voom_mean_variance.png"), width = 7, height = 6, units = "in", res = 300)
invisible(voom(y_all, model.matrix(~ 0 + Status + Age_days, data = meta_host), plot = TRUE))
invisible(dev.off())
cat("Saved: QC_MDS_status.png, QC_MDS_farm.png, QC_voom_mean_variance.png\n")

section("FARM BALANCE")
cat("Farms:", n_distinct(meta_host$Farm_Number), "| calves per farm: median",
    median(table(meta_host$Farm_Number)), ", range", paste(range(table(meta_host$Farm_Number)), collapse = "-"), "\n")
cat("Farms by the clinical statuses present:\n")
meta_host %>% group_by(Farm_Number) %>%
  summarise(statuses = paste(sort(unique(Status)), collapse = " + "), .groups = "drop") %>%
  count(statuses, name = "farms") %>% arrange(desc(farms)) %>% as.data.frame() %>% print(row.names = FALSE)


############################################################
# 6 — DIFFERENTIAL EXPRESSION
############################################################

section("DIFFERENTIAL EXPRESSION")

# Fit one limma-voom model (farm as block) and extract a set of contrasts.
# contrast_defs: named list; each element has weights (named by design column),
# g1/g2 (reference / test group labels) and n1/n2 (calves).
fit_contrasts <- function(y, meta, design, contrast_defs, pathogen_id, pathogen_clean_name) {
  stopifnot(identical(colnames(y), meta$Sample_id))
  farm   <- factor(meta$Farm_Number)
  # two passes (limma user guide): the farm correlation changes the voom weights
  v      <- voom(y, design, plot = FALSE)
  corfit <- duplicateCorrelation(v, design, block = farm)
  v      <- voom(y, design, plot = FALSE, block = farm, correlation = corfit$consensus)
  corfit <- duplicateCorrelation(v, design, block = farm)
  fit    <- lmFit(v, design, block = farm, correlation = corfit$consensus)
  cm <- sapply(contrast_defs, function(cd) {
    w <- setNames(rep(0, ncol(design)), colnames(design))
    if (!all(names(cd$weights) %in% colnames(design)))
      stop("contrast columns not in design: ", paste(setdiff(names(cd$weights), colnames(design)), collapse = ", "))
    w[names(cd$weights)] <- cd$weights
    w
  })
  cm <- matrix(cm, nrow = ncol(design), dimnames = list(colnames(design), names(contrast_defs)))
  fit2 <- eBayes(contrasts.fit(fit, cm))
  res <- bind_rows(lapply(names(contrast_defs), function(cn) {
    cd <- contrast_defs[[cn]]
    topTable(fit2, coef = cn, number = Inf, sort.by = "P") %>%
      rownames_to_column("gene") %>%
      dplyr::rename(PValue = P.Value, FDR = adj.P.Val) %>%
      mutate(pathogen = pathogen_id, pathogen_clean = pathogen_clean_name, comparison = cn,
             group_1_name = cd$g1, group_2_name = cd$g2,
             n_group_1 = cd$n1, n_group_2 = cd$n2, farm_correlation = corfit$consensus)
  }))
  list(res = res, farm_correlation = corfit$consensus)
}

log_row <- function(pathogen_id, comparison, g1 = NA, n1 = NA, g2 = NA, n2 = NA,
                    status = "success", reason = NA, farm_correlation = NA_real_, genes_tested = NA_integer_) {
  tibble(pathogen = pathogen_id,
         pathogen_clean = if (pathogen_id == "none") "all" else clean_pathogen_names(pathogen_id),
         comparison = comparison, group_1 = g1, n_group_1 = as.integer(n1),
         group_2 = g2, n_group_2 = as.integer(n2), status = status, reason = reason,
         farm_correlation = farm_correlation, genes_tested = as.integer(genes_tested))
}

# filter genes for the model, fit, store and log; a failed fit is logged and skipped
#   group_sizes: sizes of the groups the model compares (sets the gene filter)
run_fit <- function(key, meta, design, defs, p, group_sizes) {
  if (qr(design)$rank < ncol(design)) {
    for (cn in names(defs)) comparison_log[[length(comparison_log) + 1]] <<-
      log_row(p, cn, defs[[cn]]$g1, defs[[cn]]$n1, defs[[cn]]$g2, defs[[cn]]$n2,
              status = "skipped", reason = "design not estimable")
    return(invisible(NULL))
  }
  prep <- tryCatch(prepare_y(meta$Sample_id, group_sizes),
                   error = function(e) { message(p, " / ", key, ": ", conditionMessage(e)); NULL })
  out <- if (is.null(prep)) NULL else
    tryCatch(fit_contrasts(prep$y, meta, design, defs, p, if (p == "none") "all" else clean_pathogen_names(p)),
             error = function(e) { message(p, " / ", key, ": model failed — ", conditionMessage(e)); NULL })
  if (is.null(out)) {
    for (cn in names(defs)) comparison_log[[length(comparison_log) + 1]] <<-
      log_row(p, cn, defs[[cn]]$g1, defs[[cn]]$n1, defs[[cn]]$g2, defs[[cn]]$n2,
              status = "skipped", reason = "model failed")
    return(invisible(NULL))
  }
  results_list[[paste(p, key, sep = " | ")]] <<- out$res
  for (cn in names(defs)) comparison_log[[length(comparison_log) + 1]] <<-
    log_row(p, cn, defs[[cn]]$g1, defs[[cn]]$n1, defs[[cn]]$g2, defs[[cn]]$n2,
            farm_correlation = out$farm_correlation, genes_tested = prep$n_genes)
  out$n_genes <- prep$n_genes; out$n_min <- prep$n_min
  invisible(out)
}

results_list   <- list()
comparison_log <- list()

# ---- disease model: respiratory vs healthy, taxon detection ignored ----
dis <- meta_host %>% mutate(group = factor(Status, levels = status_levels_de))
n_dis <- table(dis$group)
out <- run_fit("disease", dis, model.matrix(~ 0 + group + Age_days, data = dis),
               list(respiratory_vs_healthy = list(weights = c(grouprespiratory = 1, grouphealthy = -1),
                                                  g1 = "healthy", g2 = "respiratory",
                                                  n1 = n_dis[["healthy"]], n2 = n_dis[["respiratory"]])),
               "none", group_sizes = n_dis)
if (is.null(out))
  fail("Disease model (respiratory_vs_healthy) not fitted: design not estimable, too few genes ",
       "after filtering, or model error (see the message above, if any).")
cat("Disease model fitted — genes tested:", out$n_genes, "| farm correlation:", round(out$farm_correlation, 3), "\n")

# ---- per-taxon models ----
cell_levels <- as.vector(t(outer(status_levels_de, c("negative", "positive"), paste, sep = "_")))
cell_counts <- list()

for (p in unique(pathogen_classified$pathogen)) {
  dat <- pathogen_classified %>%
    filter(pathogen == p, Sample_id %in% common_samples, class %in% c("negative", "positive")) %>%
    distinct(Sample_id, .keep_all = TRUE) %>%
    dplyr::select(Sample_id, class) %>%
    inner_join(meta_host %>% dplyr::select(Sample_id, Status, Age_days, Farm_Number), by = "Sample_id") %>%
    arrange(match(Sample_id, colnames(y_all))) %>%
    mutate(Status    = factor(Status, levels = status_levels_de),
           detection = factor(class, levels = c("negative", "positive")),
           cell      = paste(Status, class, sep = "_"))

  n_cell <- table(factor(dat$cell, levels = cell_levels))
  cell_counts[[p]] <- as_tibble(as.list(n_cell)) %>% mutate(pathogen = p, .before = 1)
  n_pos <- sum(dat$detection == "positive")
  n_neg <- sum(dat$detection == "negative")

  if (n_pos < min_positive_taxon) {
    comparison_log[[length(comparison_log) + 1]] <-
      log_row(p, "positive_vs_negative", "negative", n_neg, "positive", n_pos,
              status = "skipped", reason = paste0("fewer than ", min_positive_taxon, " positive calves"))
    next
  }

  # positive vs negative, adjusted for clinical status
  out_a <- run_fit("additive", dat,
                   model.matrix(~ Status + detection + Age_days, data = dat),
                   list(positive_vs_negative = list(weights = c(detectionpositive = 1),
                                                    g1 = "negative (adjusted for status)", n1 = n_neg,
                                                    g2 = "positive (adjusted for status)", n2 = n_pos)),
                   p, group_sizes = c(n_neg, n_pos))

  cat("Fitted:", taxon_label(p), "| positive:", n_pos,
      "| genes tested:", if (is.null(out_a)) "NA" else paste0(out_a$n_genes, " (>= ", out_a$n_min, " calves)"),
      "| farm correlation:", if (is.null(out_a)) "NA" else round(out_a$farm_correlation, 3), "\n")
}

cat("\nCalves per status x taxon class (excluded low-abundance calves not counted):\n")
bind_rows(cell_counts) %>%
  mutate(taxon = taxon_label(pathogen), positive_total = healthy_positive + respiratory_positive, .before = 1) %>%
  dplyr::select(-pathogen) %>% arrange(desc(positive_total)) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nFarm clustering of taxon detection (calves in the models):\n")
pathogen_classified %>%
  filter(Sample_id %in% common_samples, class %in% c("negative", "positive")) %>%
  group_by(taxon, Farm_Number) %>%
  summarise(pos = sum(class == "positive"), neg = sum(class == "negative"), .groups = "drop") %>%
  group_by(taxon) %>%
  summarise(positive_calves = sum(pos), farms_with_positives = sum(pos > 0),
            farms_with_pos_and_neg = sum(pos > 0 & neg > 0),
            pct_positives_on_top3_farms = round(100 * sum(sort(pos, decreasing = TRUE)[1:min(3, n())]) /
                                                  max(sum(pos), 1)),
            .groups = "drop") %>%
  arrange(desc(positive_calves)) %>%
  as.data.frame() %>% print(row.names = FALSE)

if (length(results_list) == 0) fail("No DE comparisons were estimable.")

DE_all <- bind_rows(results_list) %>%
  as_tibble() %>%
  mutate(taxon = if_else(pathogen == "none", "all (taxon ignored)", taxon_label(pathogen)))

DE_sig      <- DE_all %>% filter(FDR < fdr_cutoff, abs(logFC) > lfc_cutoff)
DE_sig_high <- DE_all %>% filter(FDR < fdr_cutoff, abs(logFC) > lfc_cutoff_high)

comparison_summary <- bind_rows(comparison_log) %>% as_tibble() %>% mutate(tested = status == "success")


############################################################
# 7 — DE OUTPUT FILES
############################################################

write_csv(DE_all %>% dplyr::select(-taxon),
          file.path(out_dir, "host_DE_all_results.csv"))
write_csv(DE_sig %>% dplyr::select(-taxon),
          file.path(out_dir, "host_DE_significant_genes.csv"))
write_csv(DE_sig_high %>% dplyr::select(-taxon),
          file.path(out_dir, "host_DE_significant_genes_high_FC.csv"))
write_csv(comparison_summary, file.path(out_dir, "host_expression_sample_summary.csv"))


############################################################
# 8 — DE RESULTS SUMMARIES
############################################################

section("TESTS RUN AND SKIPPED")
comparison_summary %>% count(comparison, status, reason) %>% as.data.frame() %>% print(row.names = FALSE)
cat("\nGroup sizes per test (successful runs):\n")
comparison_summary %>%
  filter(status == "success") %>%
  dplyr::select(pathogen_clean, comparison, group_1, n_group_1, group_2, n_group_2, genes_tested, farm_correlation) %>%
  mutate(farm_correlation = round(farm_correlation, 3)) %>%
  arrange(comparison, pathogen_clean) %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\nSkipped tests:\n")
comparison_summary %>%
  filter(status != "success") %>%
  dplyr::select(pathogen_clean, comparison, group_1, n_group_1, group_2, n_group_2, reason) %>%
  arrange(comparison, pathogen_clean) %>%
  as.data.frame() %>% print(row.names = FALSE)

section("DE GENE COUNTS PER TAXON AND COMPARISON")
cat("Genes tested per test (filtered per model): range",
    paste(range(count(DE_all, pathogen, comparison)$n), collapse = " to "), "\n")
cat("DE at FDR <", fdr_cutoff, "with |log2FC| >", lfc_cutoff, "(DE_*) and >", lfc_cutoff_high, "(high_*)\n",
    "FDR is controlled within each test (taxon x comparison), across all filtered genes\n\n")
de_summary <- DE_all %>%
  group_by(pathogen_clean, comparison) %>%
  summarise(
    genes_tested = n(),
    FDR_only  = sum(FDR < fdr_cutoff),
    DE        = sum(FDR < fdr_cutoff & abs(logFC) > lfc_cutoff),
    DE_up     = sum(FDR < fdr_cutoff & logFC > lfc_cutoff),
    DE_down   = sum(FDR < fdr_cutoff & logFC < -lfc_cutoff),
    high      = sum(FDR < fdr_cutoff & abs(logFC) > lfc_cutoff_high),
    high_up   = sum(FDR < fdr_cutoff & logFC > lfc_cutoff_high),
    high_down = sum(FDR < fdr_cutoff & logFC < -lfc_cutoff_high),
    .groups = "drop"
  ) %>%
  arrange(comparison, desc(DE))
write_csv(de_summary, file.path(out_dir, "DE_summary_by_taxon_comparison.csv"))
print(as.data.frame(de_summary), row.names = FALSE)


############################################################
# 9 — HAND-OFF TO 05a, 05b
############################################################

saveRDS(list(
  params = mget(params_03a),
  y_all = y_all,                          # filtered, TMM-normalised DGEList (DE calves)
  meta_host = meta_host,                  # DE calves
  meta_seq = meta_seq,                    # all sequenced calves
  pathogen_classified = pathogen_classified,
  logcpm_disp = logcpm_disp,              # TMM log2 CPM, all sequenced calves (figures)
  DE_all = DE_all,
  comparison_summary = comparison_summary,
  created = Sys.time()
), file.path(out_dir, "03a_objects.rds"))
cat("\nSaved: 03a_objects.rds (hand-off to 05a, 05b)\n")

section("SESSION INFORMATION")
print(sessionInfo())
cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

sink(type = "message")
sink()
close(log_con)
