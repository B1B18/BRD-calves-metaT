############################################################
# 04a_networks_DE_genes.R
#
# Host–taxon correlation networks, one per clinical status
# (healthy, respiratory). Computes the networks only; the network figures
# are drawn in 05b_figures.R.
#
#   Calves   : all calves flagged in metadata.csv as having host, viral and
#              bacterial data (has_host, has_virus, has_bacteria), split by
#              clinical status. Eligibility comes from the metadata flags, not
#              from taxa_long.csv: a calf with viral (bacterial) data but no row
#              for any target virus (bacterium) has zero abundance for all of
#              them and is kept (counts written to the log).
#   Features : target viral species (normalised abundance), selected
#              bacterial genera (log10 reads per 1e8 paired reads + 1),
#              host genes (log2 TMM-normalised CPM + 1). Correlations are
#              Spearman (rank-based), so these transformations do not change rho.
#   Taxa     : kept if detected in >= MIN_PATHOGEN_SAMPLES calves of the status group
#              (viruses: normalised abundance > VIRUS_PREV_MIN; bacteria: any read)
#   Genes    : all host genes, kept in a status group if TMM-CPM > 1 in
#              >= MIN_GENE_CALVES calves of that group (independent of the DE
#              filtering, so genes expressed only in a few calves, e.g. those
#              carrying a taxon, are retained). Genes are NOT pre-selected by DE.
#   Edges    : Spearman correlation of each taxon with every other taxon
#              and host gene; BH adjustment over all pairs tested within a
#              network; kept if |rho| >= COR_THRESH and adjusted p <= P_THRESH.
#              Correlations are unadjusted (marginal): not adjusted for farm,
#              age or other covariates.
#   Stability: kept if rho has the same sign and |rho| >= COR_THRESH in
#              >= STABILITY_THRESHOLD of N_BOOTSTRAP_ITER bootstrap resamples
#              of the calves (doRNG, BOOTSTRAP_SEED)
#   Strict   : within the network, edges also stable at COR_THRESH_STRICT are flagged
#              (stable_strict: |rho| >= COR_THRESH_STRICT and same sign with
#              |rho| >= COR_THRESH_STRICT in >= STABILITY_THRESHOLD of the SAME
#              resamples; stability_strict = that fraction). The flag does not change
#              the network: an edge stable at COR_THRESH_STRICT is also stable at COR_THRESH.
#   DE overlap (two separate, descriptive summaries; not enrichment tests)
#     gene level : correlated genes (>= 1 stable taxon-gene edge) that were
#                  tested in 03a and are DE at FDR < DE_FDR in any comparison,
#                  no fold-change threshold; shown next to the same proportion
#                  for all genes tested in that network (background)
#     edge level : each stable taxon-gene edge matched to the 03a
#                  positive_vs_negative test for the SAME taxon and gene; counted
#                  only where that gene was tested in that taxon's model; DE at
#                  FDR < DE_FDR, and direction (sign of log2FC = sign of rho)
#
# Inputs
#   input_files/metadata.csv, input_files/taxa_long.csv, input_files/host_counts.csv.gz
#   03a_host_DE/host_DE_all_results.csv    all DE results
#
# Outputs (04a_networks_DE_genes/)
#   log.txt                                   console log, numbers for the
#                                             Results, session information
#   network_<status>_stable_bootstrap_edges.csv
#   network_<status>_stable_bootstrap_nodes.csv  node type (virus / bacteria / host_gene)
#   network_edges_all_statuses.csv            all stable edges, all networks
#   network_parameters.csv                    thresholds used (read by 05b for the figure title)
#   network_genes_DE_overlap.csv              gene level: correlated genes also DE (FDR only),
#                                             with the network background
#   network_edges_DE_same_taxon.csv           edge level: each taxon-gene edge with the
#                                             03a result for the same taxon and gene
#   network_edges_DE_same_taxon_summary.csv   edge level: counts per network and taxon
############################################################


############################################################
# 0 — SETUP
############################################################

# Run from the repository root: all paths below are relative to it.

suppressPackageStartupMessages({
  library(edgeR)
  library(psych)
  library(foreach)
  library(doParallel)
  library(doRNG)
  library(tidyverse)
})

out_dir <- "04a_networks_DE_genes"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
in_dir  <- "input_files"
de_dir  <- "03a_host_DE"            # DE results (03a_host_DE.R)

# ---- parameters ----
STATUS_MODES         <- c("healthy", "respiratory")
METHOD               <- "spearman"
COR_THRESH           <- 0.4       # network edges
COR_THRESH_STRICT    <- 0.5       # stricter level, flagged within the network (bold in Figure 7)
P_THRESH             <- 0.05      # BH-adjusted p-value
BOOTSTRAP_SEED       <- 555
N_BOOTSTRAP_ITER     <- 1000
STABILITY_THRESHOLD  <- 0.8
MIN_PATHOGEN_SAMPLES <- 5
MIN_GENE_CALVES      <- 5        # gene kept if CPM > 1 in >= this many calves of the status
VIRUS_PREV_MIN       <- 5        # virus detected if normalised abundance > 5 (as in 01, 02)
BACT_PREV_MIN        <- log10(0.0001 + 1)  # bacteria on log10 scale: detected = any read
DE_FDR               <- 0.05     # DE overlap: FDR only, no fold-change threshold

selected_bacteria <- c(
  "Mesomycoplasma", "Moraxella", "Moraxella_A", "Mannheimia",
  "Mycoplasmopsis_A", "Pasteurella", "Streptococcus",
  "Psychrobacter", "Burkholderia", "Prevotella", "Alloprevotella"
)

virus_targets <- c(
  "Bovine adenovirus 3", "Bovine coronavirus", "Bovine nebovirus", "Bovine picornavirus",
  "Bovine rhinitis A", "Bovine rhinitis B", "Bovine torovirus", "BVDV 1",
  "Influenza D virus"
)

de_all_file <- file.path(de_dir, "host_DE_all_results.csv")

# ---- log ----
options(width = 250, readr.show_progress = FALSE)
log_con <- file(file.path(out_dir, "log.txt"), open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat("04a_networks_DE_genes.R\n")
cat("Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("\nPARAMETERS\n")
for (nm in c("STATUS_MODES", "METHOD", "COR_THRESH", "COR_THRESH_STRICT", "P_THRESH", "BOOTSTRAP_SEED",
             "N_BOOTSTRAP_ITER", "STABILITY_THRESHOLD", "MIN_PATHOGEN_SAMPLES", "MIN_GENE_CALVES",
             "VIRUS_PREV_MIN", "BACT_PREV_MIN", "DE_FDR", "selected_bacteria", "virus_targets"))
  cat(sprintf("  %-22s = %s\n", nm, paste(get(nm), collapse = ", ")))


############################################################
# 1 — HELPERS
############################################################

fail <- function(...) {
  msg <- paste0(...)
  cat("\nERROR:", msg, "\n")
  sink(type = "message"); sink()
  stop(msg, call. = FALSE)
}

section <- function(title) {
  options(width = 250, max.print = 1e6)
  cat("\n\n============================================================\n")
  cat(" ", title, "\n")
  cat("============================================================\n")
}

# taxa x samples matrix from taxa_long.csv (only samples with >= 1 row for these taxa)
taxa_matrix <- function(taxa, layer) {
  d <- taxa_long %>% filter(type == layer, taxon %in% taxa)
  missing <- setdiff(taxa, d$taxon)
  if (length(missing) > 0) fail("Taxa absent from taxa_long.csv (", layer, "): ", paste(missing, collapse = ", "))
  d %>% dplyr::select(taxon, Sample_id, value) %>%
    pivot_wider(names_from = Sample_id, values_from = value, values_fill = 0) %>%
    column_to_rownames("taxon") %>% as.matrix() %>% .[taxa, , drop = FALSE]
}

# add eligible samples with no rows for any of the taxa (zero abundance) and
# return the matrix with columns in the order of `samples`
fill_samples <- function(m, samples) {
  miss <- setdiff(samples, colnames(m))
  if (length(miss) > 0)
    m <- cbind(m, matrix(0, nrow(m), length(miss), dimnames = list(rownames(m), miss)))
  m[, samples, drop = FALSE]
}

# same taxon names as in 03a (viral abbreviations from taxa_long.csv, genus without g__)
taxon_label <- function(x) {
  lab <- virus_labels[str_trim(x)]
  unname(if_else(is.na(lab), str_remove(str_trim(x), "^g__"), lab))
}


############################################################
# 2 — METADATA, TAXA AND HOST EXPRESSION
############################################################

section("DATA")

metadata <- read_csv(file.path(in_dir, "metadata.csv"), col_types = cols(
  Sample_id = col_character(), Status = col_character(), Age_days = col_double(),
  Age_group = col_character(), Farm_Number = col_character(), Total_paired_reads = col_double(),
  has_virus = col_logical(), has_bacteria = col_logical(), has_host = col_logical()))
if (anyDuplicated(metadata$Sample_id)) fail("Duplicate Sample_id values in metadata.csv.")
taxa_long <- read_csv(file.path(in_dir, "taxa_long.csv"), col_types = cols(
  Sample_id = col_character(), taxon = col_character(), taxon_label = col_character(),
  type = col_character(), value = col_double()))
n_dup <- taxa_long %>% count(type, taxon, Sample_id) %>% filter(n > 1) %>% nrow()
if (n_dup > 0) fail("Duplicate type/taxon/Sample_id rows in taxa_long.csv: ", n_dup)
virus_labels <- taxa_long %>% filter(type == "virus") %>% distinct(taxon, taxon_label) %>%
  { setNames(.$taxon_label, .$taxon) }

virus_mat <- taxa_matrix(virus_targets, "virus")          # viruses x samples, normalised abundance
bact_mat  <- taxa_matrix(selected_bacteria, "bacteria")   # genera x samples, read counts

host_mat <- read_csv(file.path(in_dir, "host_counts.csv.gz"),
                     col_types = cols(gene_id = col_character(), .default = col_double())) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()
if (anyNA(host_mat)) fail("NA values in host_counts.csv.gz.")

# ---- eligible calves: from the metadata flags ----
eligible <- metadata %>%
  filter(coalesce(has_host, FALSE), coalesce(has_virus, FALSE), coalesce(has_bacteria, FALSE))
common_samples <- eligible$Sample_id
if (length(common_samples) < 10) fail("Too few eligible calves: ", length(common_samples))

no_host <- setdiff(common_samples, colnames(host_mat))
if (length(no_host) > 0)
  fail("Calves flagged has_host but absent from host_counts.csv.gz: ", paste(no_host, collapse = ", "))

absent_virus <- setdiff(common_samples, colnames(virus_mat))
absent_bact  <- setdiff(common_samples, colnames(bact_mat))
cat("Calves with host, viral and bacterial data (metadata flags):", length(common_samples), "\n")
cat("  of which no row for any target virus in taxa_long.csv (zero abundance):", length(absent_virus), "\n")
cat("  of which no row for any selected bacterium in taxa_long.csv (zero reads):", length(absent_bact), "\n")
cat("Calves in taxa_long.csv but not flagged as having viral / bacterial data:",
    length(setdiff(colnames(virus_mat), metadata$Sample_id[coalesce(metadata$has_virus, FALSE)])), "/",
    length(setdiff(colnames(bact_mat),  metadata$Sample_id[coalesce(metadata$has_bacteria, FALSE)])), "\n")
print(eligible %>% count(Status))

reads <- metadata$Total_paired_reads[match(common_samples, metadata$Sample_id)]
if (any(!is.finite(reads) | reads <= 0)) fail("Missing or non-positive Total_paired_reads.")

virus_df <- t(fill_samples(virus_mat, common_samples))
colnames(virus_df) <- taxon_label(colnames(virus_df))
bact_df  <- t(log10(sweep(fill_samples(bact_mat, common_samples), 2, reads, "/") * 1e8 + 1))
host_mat <- host_mat[, common_samples, drop = FALSE]
stopifnot(identical(rownames(virus_df), common_samples),
          identical(rownames(bact_df), common_samples))

# ---- host genes: all genes with any reads; filtered per status in run_network ----
# TMM on all genes and calves
gene_universe <- rownames(host_mat)[rowSums(host_mat) > 0]
cpm_mat  <- cpm(calcNormFactors(DGEList(counts = host_mat[gene_universe, , drop = FALSE]), method = "TMM"))
host_log <- log2(cpm_mat + 1)
cat("\nHost genes with any reads:", length(gene_universe), "of", nrow(host_mat), "\n")

features_all <- cbind(virus_df, bact_df, t(host_log))
status_all   <- metadata$Status[match(common_samples, metadata$Sample_id)]

if (COR_THRESH_STRICT < COR_THRESH) fail("COR_THRESH_STRICT must be >= COR_THRESH.")
if (!file.exists(de_all_file)) fail("DE results from 03a not found: ", de_all_file)
DE_all <- read_csv(de_all_file, show_col_types = FALSE)


############################################################
# 3 — NETWORK FOR ONE CLINICAL STATUS
############################################################

run_network <- function(status) {

  section(paste("NETWORK —", status))

  dat <- features_all[which(status_all == status), , drop = FALSE]
  dat[is.na(dat)] <- 0
  if (nrow(dat) < 10) { cat("Too few calves:", nrow(dat), "\n"); return(NULL) }

  # ---- taxa detected in >= MIN_PATHOGEN_SAMPLES calves ----
  virus_prev <- colSums(dat[, colnames(virus_df), drop = FALSE] > VIRUS_PREV_MIN)
  bact_prev  <- colSums(dat[, colnames(bact_df),  drop = FALSE] > BACT_PREV_MIN)
  prevalence <- c(virus_prev, bact_prev)
  pathogen_features <- names(prevalence)[prevalence >= MIN_PATHOGEN_SAMPLES]
  virus_features    <- intersect(colnames(virus_df), pathogen_features)
  bacteria_features <- intersect(colnames(bact_df),  pathogen_features)

  cat("Calves:", nrow(dat), "\n")
  cat("Calves with each taxon (kept if >=", MIN_PATHOGEN_SAMPLES, "):\n")
  print(tibble(taxon = names(prevalence), calves = prevalence, kept = names(prevalence) %in% pathogen_features), n = Inf)

  # ---- host genes: CPM > 1 (log2 CPM + 1 > 1) in >= MIN_GENE_CALVES calves ----
  host_keep <- gene_universe[colSums(dat[, gene_universe, drop = FALSE] > 1) >= MIN_GENE_CALVES]
  cat("Host genes tested:", length(host_keep), "of", length(gene_universe), "\n")
  if (length(pathogen_features) == 0) { cat("No taxa retained.\n"); return(NULL) }

  dat_pathogen <- dat[, pathogen_features, drop = FALSE]
  dat_all      <- dat[, c(pathogen_features, host_keep), drop = FALSE]
  stopifnot(identical(rownames(dat_pathogen), rownames(dat_all)))

  # ---- correlations: taxa x (taxa + host genes); BH over all pairs, one row per pair ----
  ct <- psych::corr.test(dat_pathogen, dat_all, method = METHOD, adjust = "none", ci = FALSE)

  edges_full <- as.data.frame(ct$r) %>%
    rownames_to_column("from") %>%
    pivot_longer(-from, names_to = "to", values_to = "corr") %>%
    mutate(p = as.vector(t(ct$p))) %>%
    # self-pairs and the second ordering of each taxon-taxon pair removed BEFORE the
    # BH adjustment, so the correction counts each pair tested exactly once
    filter(from != to) %>%
    mutate(edge_key = paste(pmin(from, to), pmax(from, to), sep = "_")) %>%
    distinct(edge_key, .keep_all = TRUE) %>%
    filter(!is.na(corr)) %>%
    mutate(p_adj = p.adjust(p, method = "BH")) %>%
    arrange(p_adj, desc(abs(corr)))

  edges <- edges_full %>% filter(abs(corr) >= COR_THRESH, p_adj <= P_THRESH)
  cat("Pairs tested:", nrow(edges_full), "| edges with |rho| >=", COR_THRESH, "and BH p <=", P_THRESH, ":", nrow(edges), "\n")
  if (nrow(edges) == 0) return(NULL)

  # ---- bootstrap stability: same sign and |rho| >= COR_THRESH (network) and
  #      >= COR_THRESH_STRICT (flag), counted on the same resamples ----
  n_tax <- ncol(dat_pathogen)
  # parameters copied to local variables: PSOCK workers receive the function's local
  # objects automatically, but not global ones (METHOD, COR_THRESH, COR_THRESH_STRICT)
  cor_method <- METHOD
  cor_thresh <- COR_THRESH
  cor_strict <- COR_THRESH_STRICT
  boot_vars  <- c("dat_pathogen", "dat_all", "cor_method", "cor_thresh", "cor_strict")
  set.seed(BOOTSTRAP_SEED)
  counts <- foreach(i = seq_len(N_BOOTSTRAP_ITER), .combine = "+",
                    .export = boot_vars, .noexport = setdiff(ls(), boot_vars)) %dorng% {
    idx <- sample(nrow(dat_pathogen), replace = TRUE)
    rc  <- suppressWarnings(cor(dat_pathogen[idx, , drop = FALSE], dat_all[idx, , drop = FALSE], method = cor_method))
    rc[is.na(rc)] <- 0                                     # constant in a resample: not passed
    # rows: positive / negative pass at cor_thresh, then positive / negative pass at cor_strict
    rbind((rc >= cor_thresh) * 1, (rc <= -cor_thresh) * 1,
          (rc >= cor_strict) * 1, (rc <= -cor_strict) * 1)
  }
  blk <- function(k) {
    m <- counts[(k - 1) * n_tax + seq_len(n_tax), , drop = FALSE] / N_BOOTSTRAP_ITER
    dimnames(m) <- list(colnames(dat_pathogen), colnames(dat_all))
    m
  }
  pos_stab <- blk(1); neg_stab <- blk(2); pos_strict <- blk(3); neg_strict <- blk(4)

  final_edges <- edges %>%
    mutate(stability        = if_else(corr > 0, pos_stab[cbind(from, to)],   neg_stab[cbind(from, to)]),
           stability_strict = if_else(corr > 0, pos_strict[cbind(from, to)], neg_strict[cbind(from, to)]),
           stable_strict    = abs(corr) >= COR_THRESH_STRICT & stability_strict >= STABILITY_THRESHOLD) %>%
    filter(stability >= STABILITY_THRESHOLD) %>%
    dplyr::select(-edge_key)

  cat("Stable edges (same sign, |rho| >=", COR_THRESH, "in >=", STABILITY_THRESHOLD * 100, "% of resamples):",
      nrow(final_edges), "| of which also stable at |rho| >=", COR_THRESH_STRICT, ":", sum(final_edges$stable_strict), "\n")
  if (nrow(final_edges) == 0) return(NULL)

  nodes <- tibble(name = unique(c(final_edges$from, final_edges$to))) %>%
    mutate(type = case_when(name %in% virus_features    ~ "virus",
                            name %in% bacteria_features ~ "bacteria",
                            name %in% host_keep         ~ "host_gene",
                            TRUE                        ~ "other"),
           gene_clean = str_remove(name, "^gene-"))

  prefix <- file.path(out_dir, paste0("network_", status))
  write_csv(final_edges, paste0(prefix, "_stable_bootstrap_edges.csv"))
  write_csv(nodes,       paste0(prefix, "_stable_bootstrap_nodes.csv"))
  cat("Saved:", basename(paste0(prefix, "_stable_bootstrap_edges.csv")), "and nodes\n")

  list(edges = final_edges, nodes = nodes, n_calves = nrow(dat), host_genes = length(host_keep),
       host_keep = host_keep, taxa = pathogen_features, pairs_tested = nrow(edges_full))
}


############################################################
# 4 — RUN
############################################################

cl <- makeCluster(max(1L, min(8L, parallel::detectCores() - 1L)))
registerDoParallel(cl)
networks <- setNames(lapply(STATUS_MODES, run_network), STATUS_MODES)
stopCluster(cl)
networks <- networks[!vapply(networks, is.null, logical(1))]
if (length(networks) == 0) fail("No network had stable edges.")


############################################################
# 5 — RESULTS SUMMARIES
############################################################

taxa_all <- c(colnames(virus_df), colnames(bact_df))
edges_annot <- bind_rows(lapply(names(networks), function(st) networks[[st]]$edges %>% mutate(status = st))) %>%
  mutate(
    taxon     = if_else(from %in% taxa_all, from, to),
    partner   = if_else(from %in% taxa_all, to, from),
    edge_type = if_else(from %in% taxa_all & to %in% taxa_all, "taxon-taxon", "taxon-gene"),
    gene      = if_else(edge_type == "taxon-gene", partner, NA_character_)
  )
write_csv(edges_annot, file.path(out_dir, "network_edges_all_statuses.csv"))
write_csv(tibble(parameter = c("METHOD", "COR_THRESH", "COR_THRESH_STRICT", "P_THRESH",
                               "STABILITY_THRESHOLD", "N_BOOTSTRAP_ITER", "BOOTSTRAP_SEED"),
                 value = as.character(c(METHOD, COR_THRESH, COR_THRESH_STRICT, P_THRESH,
                                        STABILITY_THRESHOLD, N_BOOTSTRAP_ITER, BOOTSTRAP_SEED))),
          file.path(out_dir, "network_parameters.csv"))

section("SUMMARY PER NETWORK")
cat("Edges are unadjusted Spearman correlations (not adjusted for farm or age).\n")
for (st in names(networks)) {
  e <- edges_annot %>% filter(status == st)
  cat("\n---", st, "---\n")
  cat("Calves:", networks[[st]]$n_calves, "| host genes tested:", networks[[st]]$host_genes,
      "| pairs tested:", networks[[st]]$pairs_tested,
      "\nTaxa tested:", paste(networks[[st]]$taxa, collapse = ", "), "\n")
  cat("Stable edges:", nrow(e), "| taxon-gene:", sum(e$edge_type == "taxon-gene"),
      "| taxon-taxon:", sum(e$edge_type == "taxon-taxon"),
      "| genes correlated with >= 1 taxon:", n_distinct(e$gene[e$edge_type == "taxon-gene"]), "\n")
  cat("Also stable at |rho| >=", COR_THRESH_STRICT, ":", sum(e$stable_strict),
      "| taxon-gene:", sum(e$stable_strict & e$edge_type == "taxon-gene"),
      "| taxon-taxon:", sum(e$stable_strict & e$edge_type == "taxon-taxon"),
      "| genes:", n_distinct(e$gene[e$stable_strict & e$edge_type == "taxon-gene"]), "\n")
  cat("Edges per taxon (taxon-gene; positive / negative; also stable at |rho| >=", COR_THRESH_STRICT, "):\n")
  e %>% filter(edge_type == "taxon-gene") %>%
    group_by(taxon) %>%
    summarise(edges = n(), positive = sum(corr > 0), negative = sum(corr < 0),
              stable_strict = sum(stable_strict), .groups = "drop") %>%
    arrange(desc(edges)) %>% print(n = Inf)
  if (any(e$edge_type == "taxon-taxon")) {
    cat("Taxon-taxon edges:\n")
    e %>% filter(edge_type == "taxon-taxon") %>%
      transmute(from, to, rho = round(corr, 2), p_adj = signif(p_adj, 2), stability,
                stability_strict, stable_strict) %>% print(n = Inf)
  }
}

genes_by_status <- edges_annot %>% filter(edge_type == "taxon-gene") %>% distinct(status, gene)
cat("\nGenes correlated in both networks:",
    genes_by_status %>% count(gene) %>% filter(n > 1) %>% nrow(), "\n")


############################################################
# 6a — GENE LEVEL: CORRELATED GENES ALSO DE IN 03a (FDR only)
############################################################

section(paste0("GENE LEVEL: CORRELATED GENES ALSO DE IN 03a (FDR < ", DE_FDR, ", no fold-change threshold)"))

de_fdr <- DE_all %>% filter(FDR < DE_FDR) %>% distinct(gene, comparison)
de_any_genes <- unique(de_fdr$gene)
tested_genes <- unique(DE_all$gene)
comparisons  <- sort(unique(DE_all$comparison))

# genes: correlated genes; background: all host genes tested in that network
overlap_row <- function(genes, background, label) {
  bg_tested <- background[background %in% tested_genes]
  tibble(network = label,
         correlated_genes = length(genes),
         tested_in_DE = sum(genes %in% tested_genes),
         DE_any_comparison = sum(genes %in% de_any_genes),
         !!!setNames(lapply(comparisons, function(cn) sum(genes %in% de_fdr$gene[de_fdr$comparison == cn])),
                     paste0("DE_", comparisons))) %>%
    mutate(pct_DE_of_tested = round(100 * DE_any_comparison / pmax(tested_in_DE, 1), 1),
           .after = DE_any_comparison) %>%
    mutate(background_genes = length(background),
           background_tested_in_DE = length(bg_tested),
           background_pct_DE = round(100 * sum(bg_tested %in% de_any_genes) / max(length(bg_tested), 1), 1))
}

overlap <- bind_rows(
  lapply(names(networks), function(st)
    overlap_row(genes_by_status$gene[genes_by_status$status == st], networks[[st]]$host_keep, st)),
  overlap_row(unique(genes_by_status$gene),
              unique(unlist(lapply(networks, `[[`, "host_keep"))), "any network")
)
cat("correlated_genes: genes with >= 1 stable taxon-gene edge\n",
    "tested_in_DE: of these, tested in at least one 03a model (genes filtered per model in 03a)\n",
    "DE_any_comparison / DE_<comparison>: DE in any 03a test (FDR <", DE_FDR, ", no fold-change threshold)\n",
    "pct_DE_of_tested: DE_any_comparison as % of tested_in_DE\n",
    "background_*: the same for all host genes tested in that network (descriptive comparison only;\n",
    "  correlated genes tend to be more highly expressed and variable, which also raises DE power)\n\n")
print(as.data.frame(overlap), row.names = FALSE)
write_csv(overlap, file.path(out_dir, "network_genes_DE_overlap.csv"))


############################################################
# 6b — EDGE LEVEL: SAME TAXON AND GENE IN 03a positive_vs_negative
############################################################

section(paste0("EDGE LEVEL: TAXON-GENE EDGES vs 03a positive_vs_negative FOR THE SAME TAXON (FDR < ", DE_FDR, ")"))

de_same <- DE_all %>%
  filter(comparison == "positive_vs_negative") %>%
  transmute(taxon = str_remove(pathogen_clean, "^g__"), gene, logFC_DE = logFC, FDR_DE = FDR)
if (anyDuplicated(de_same[, c("taxon", "gene")])) fail("Duplicate taxon/gene rows in the 03a positive_vs_negative results.")
taxa_modelled <- unique(de_same$taxon)

net_taxa <- unique(edges_annot$taxon[edges_annot$edge_type == "taxon-gene"])
cat("Taxa with taxon-gene edges:", paste(sort(net_taxa), collapse = ", "), "\n")
cat("  of which modelled in 03a:", paste(sort(intersect(net_taxa, taxa_modelled)), collapse = ", "), "\n")
cat("  not modelled in 03a (e.g. too few positive calves):",
    paste(sort(setdiff(net_taxa, taxa_modelled)), collapse = ", "), "\n")
cat("Note: detection definitions differ (networks: virus > ", VIRUS_PREV_MIN,
    " / any bacterial read; 03a: Q25 rule)\n\n", sep = "")

edge_overlap <- edges_annot %>%
  filter(edge_type == "taxon-gene") %>%
  dplyr::select(status, taxon, gene, corr, p_adj, stability, stable_strict) %>%
  left_join(de_same, by = c("taxon", "gene")) %>%
  mutate(taxon_modelled_in_DE = taxon %in% taxa_modelled,
         gene_tested_in_DE    = !is.na(FDR_DE),                       # tested in THIS taxon's model
         DE_same_taxon        = gene_tested_in_DE & FDR_DE < DE_FDR,
         DE_same_direction    = DE_same_taxon & sign(logFC_DE) == sign(corr))
write_csv(edge_overlap, file.path(out_dir, "network_edges_DE_same_taxon.csv"))

summarise_overlap <- function(d) d %>%
  summarise(edges                    = n(),
            edges_taxon_modelled     = sum(taxon_modelled_in_DE),
            edges_gene_tested        = sum(gene_tested_in_DE),
            DE_same_taxon            = sum(DE_same_taxon),
            DE_same_direction        = sum(DE_same_direction),
            pct_DE_of_tested         = round(100 * DE_same_taxon / pmax(edges_gene_tested, 1), 1),
            pct_same_direction_of_DE = round(100 * DE_same_direction / pmax(DE_same_taxon, 1), 1),
            .groups = "drop")

overlap_edges <- bind_rows(
  edge_overlap %>% group_by(status) %>% summarise_overlap() %>% mutate(taxon = "all", .after = status),
  edge_overlap %>% group_by(status, taxon) %>% summarise_overlap()
) %>% arrange(status, taxon != "all", desc(edges))

cat("edges: stable taxon-gene edges\n",
    "edges_taxon_modelled: the taxon has a positive_vs_negative model in 03a\n",
    "edges_gene_tested: the gene was tested in that taxon's model (03a filters genes per model)\n",
    "DE_same_taxon: gene DE for the same taxon (FDR <", DE_FDR, ", no fold-change threshold)\n",
    "DE_same_direction: of these, sign of log2FC (positive vs negative) = sign of rho\n",
    "pct_DE_of_tested: DE_same_taxon as % of edges_gene_tested\n\n")
print(as.data.frame(overlap_edges), row.names = FALSE)
write_csv(overlap_edges, file.path(out_dir, "network_edges_DE_same_taxon_summary.csv"))


############################################################
# 7 — SESSION INFORMATION
############################################################

section("SESSION INFORMATION")
print(sessionInfo())
cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

sink(type = "message")
sink()
close(log_con)
