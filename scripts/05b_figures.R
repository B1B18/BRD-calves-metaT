############################################################
# 05b_figures.R
#
# Figures 6 and 7 and the supplementary table of DE results, from the saved
# outputs of 03a (DE), 04a (correlation networks) and 05a (annotation).
# Nothing is re-fitted.
#
# Pipeline position: run last (after 03a, 04a and 05a_gene_annotation.R).
#
# Heatmaps show the sequenced calves with host, viral and bacterial data (the 03a
# display set, logcpm_disp) and a recorded clinical status (including unspecific,
# scour, respiratory & scour and unknown age); the DE models (03a) used healthy and
# respiratory calves with known age and farm. Figure 6C shows, per taxon, only the
# calves classified as taxon-negative or taxon-positive (03a Q25 rule): calves with
# low abundance (below Q25) or not classified for that taxon are not shown.
# Heatmaps show TMM-normalised log2 CPM, scaled per gene (z-score).
#
#   6A euler_up_down.png            DE genes per taxon in fig6a_comparison
#                                   (FDR < fdr_cutoff, |log2FC| > fig6a_lfc_cutoff);
#                                   fig6a_n_taxa taxa with most genes (default: all), up and down
#   6A upset_up_down.png            the same gene sets as UpSet plots (up | down): genes per exclusive
#                                   intersection of taxa, largest fig6a_upset_max_comb shown
#   6B heatmap_6B_DE_genes.png      genes with FDR < fdr_cutoff and |log2FC| > fig6b_lfc_min in
#                                   fig6b_comparison (all, or the fig6b_n_genes largest |log2FC|);
#                                   rows split up / down, columns split healthy | respiratory | other
#   6C heatmap_combined_panel.png   per taxon, genes with FDR < fdr_cutoff and
#                                   |log2FC| > fig6c_lfc_cutoff in any of fig6c_comparisons
#                                   (optionally capped per taxon, contrast and direction);
#                                   calves classified taxon-negative or taxon-positive only,
#                                   ordered by class, then clinical status, then age
#   6C fig6C_category_bubble.png    summary of the same genes: DE genes per functional category
#                                   and taxon, up | down (fig6c_bubble_comparison); genes per
#                                   taxon checked against the heatmap (differences in the log)
#   7  network_<status>.png         correlation networks from 04a, one file per status (fig7_arrangement):
#                                   every stable edge drawn pastel; edges also stable at the strict 04a
#                                   level in full colour, their genes labelled; node shape = 05a category,
#                                   taxon size = number of stable edges; height fitted to the network
#   Supplementary table             every gene tested in every DE test, with annotation
#
# Inputs
#   03a_host_DE/03a_objects.rds, 05a_gene_annotation/05a_annotation.rds,
#   05a_gene_annotation/gene_annotation_evidence.csv,
#   04a_networks_DE_genes/network_<status>_stable_bootstrap_edges.csv,
#   input_files/taxa_long.csv (viral abbreviations only)
#
# Outputs (05b_figures/)
#   log.txt, the figures above,
#   heatmap_6B_DE_genes.csv (6B genes), fig6C_genes.csv (6C genes, the contrast(s) they
#   were selected from and their DE support), top_up_down_genes_per_pathogen.csv,
#   TableS_DE_results_all.csv.gz (all genes, all tests),
#   TableS_DE_results_significant.csv (FDR < fdr_cutoff and |log2FC| > lfc_cutoff),
#   fig7_edges_<status>.csv (every edge drawn; DE_same_taxon = gene DE for the same taxon in
#   positive_vs_negative at FDR < fdr_cutoff, no fold-change threshold, as in 04a;
#   DE_same_direction = also sign of log2FC = sign of rho),
#   annotation_check_figure_genes.csv (every gene in Figures 6 and 7 with its annotation)
############################################################


############################################################
# 0 — PARAMETERS (see README)
############################################################

# Run from the repository root: all paths below are relative to it.

in_dir    <- "input_files"
de_dir    <- "03a_host_DE"
annot_dir <- "05a_gene_annotation"
net_dir   <- "04a_networks_DE_genes"
out_dir   <- "05b_figures"

# clinical statuses shown in the heatmaps, in display order (03a display set: host, viral and bacterial data)
status_levels <- c("healthy", "unspecific", "respiratory", "scour", "respiratory_scours")

# contrasts, in the order used in the supplementary table
comparison_order <- c("respiratory_vs_healthy", "positive_vs_negative")

# Figure 6A
fig6a_comparison <- "positive_vs_negative"
fig6a_lfc_cutoff <- 1
fig6a_n_taxa     <- Inf    # taxa shown per panel (most genes first); Inf = every taxon with >= 1 gene

# Figure 6B
fig6b_comparison <- "respiratory_vs_healthy"
fig6b_n_genes    <- Inf        # Inf = every gene meeting the criteria; otherwise the largest |log2FC|
fig6b_lfc_min    <- 1          # genes shown: FDR < fdr_cutoff and |log2FC| > this
fig6b_row_order  <- "logFC"    # within the up and down blocks: "logFC" = sorted by log2FC; "cluster"
fig6b_max_age    <- Inf        # Inf = all ages, including unknown age
fig6b_contrast_statuses <- c("healthy", "respiratory")   # column blocks shown side by side;
                               # all other statuses pooled in a third block to the right
fig6b_show_lfc   <- TRUE       # bar of each gene's log2FC to the right of the heatmap
fig6b_category_bar <- TRUE     # top-right inset: genes up / down per category (doubles as the
                               # category key); FALSE = standard Category legend instead

# Figure 6C
fig6c_comparisons       <- "positive_vs_negative"
fig6c_lfc_cutoff        <- 2       # genes with FDR < fdr_cutoff and |log2FC| > this, in any of fig6c_comparisons
fig6c_max_per_direction <- Inf     # Inf = all such genes; otherwise the largest |log2FC| genes
                                   # per taxon, contrast and direction
fig6c_samples           <- "all"   # "all" = all calves in the 03a display set (host, viral and bacterial
                                   # data); "DE" = calves in the DE models; either way, per taxon only the
                                   # calves classified taxon-negative or taxon-positive are shown
fig6c_max_height_px     <- 30000   # whole figure; row height shrinks to fit, row names hidden if too small

# table of candidate genes per taxon
top_table_n <- 5

# Figure 7 (correlation networks, 04a)
fig7_statuses            <- c("healthy", "respiratory")
fig7_min_abs_rho         <- 0.4   # edges drawn: |rho| >= this (taxon-gene and taxon-taxon);
                                  # 0.4 = every stable edge of the 04a networks
fig7_max_genes_per_taxon <- Inf   # taxon-gene edges drawn per taxon (strongest |rho|); Inf = all passing fig7_min_abs_rho
fig7_highlight_rule      <- "strict" # which drawn edges are in full colour (bold); all others pastel:
                                  #   "strict" = edges also stable at the stricter 04a level (column
                                  #             stable_strict: |rho| >= COR_THRESH_STRICT in >= 80% of resamples)
                                  #   "rho"    = edges with observed |rho| >= fig7_highlight_min_abs_rho
                                  #   "top"    = all taxon-taxon edges + per taxon, the fig7_highlight_per_taxon
                                  #             annotated genes with the strongest |rho|
fig7_highlight_min_abs_rho <- 0.5 # rule "rho"
fig7_highlight_per_taxon <- 5     # rule "top" (Inf = every annotated gene edge outside fig7_highlight_exclude)
fig7_highlight_exclude   <- "Other/Uncharacterised"   # rule "top": categories never highlighted
fig7_highlight_DE_only   <- FALSE # all rules: TRUE = full colour only for genes DE (FDR < fdr_cutoff, no
                                  # fold-change threshold, as in 04a) for the same taxon in positive_vs_negative
fig7_label_per_taxon     <- Inf   # gene labels: genes on full-colour edges, at most this many per taxon
                                  # (strongest |rho| first); Inf = label every full-colour gene
fig7_faded_cols          <- c(positive = "#F4C7C3", negative = "#C6DBEF", gene = "#D9D2E9")
fig7_layout              <- "packed_auto" # each connected part laid out with Fruchterman-Reingold, scaled to
                                     #   its number of nodes; the largest part first, the others (largest first):
                                     #   "packed_right" = in a column to its right; "packed_bottom" = in a row below;
                                     #   "packed_auto" = right or below, whichever fills the panel better;
                                     # "fr_packed" = parts packed at random (igraph::merge_coords);
                                     # "fr", "kk" = one layout of the whole network (unconnected parts spread out)
fig7_seed                <- 555      # layout seed (change to rotate / rearrange)
fig7_pack_gap            <- 0.06     # "packed_*": space between parts (fraction of the largest part)
fig7_mirror              <- FALSE    # "packed_*": mirror the largest part left-right
fig7_rotate_largest      <- TRUE     # "packed_*": also try the largest part turned (e.g. a long, thin part laid across)
fig7_target_aspect       <- 0.75     # "packed_*": panel shape (height / width) the arrangement is chosen for; one is
                                     #   changed from "as laid out, parts to the right" only if it fills it > 10% better
fig7_bold_edge_weight    <- 0.1      # "packed_*", "fr_packed": pull of full-colour taxon-gene edges in the layout,
                                     #   relative to 1 for all other edges; < 1 places those genes further out,
                                     #   beyond the pastel genes, so their labels have room (1 = no difference;
                                     #   edge length carries no information either way)
fig7_bold_min_genes      <- 3        # ... only for taxa with at least this many full-colour gene edges (a single
                                     #   labelled gene has room anyway and would only stretch the figure)
fig7_bold_arc            <- 300      # ... and those genes (if linked to that taxon only) spread evenly on a ring
                                     #   around it, over the largest free arc (clear of the other taxa and of shared
                                     #   genes), at most this many degrees; ordered by rho (red to blue); 0 = left
                                     #   where the layout put them
fig7_bold_ring           <- 1.3      # ring radius, relative to the radius of the taxon's pastel genes (90th centile)
fig7_expand              <- 0.06     # space around the network, each side (fraction of its width / height), for labels
fig7_fixed_aspect        <- TRUE     # TRUE = true aspect ratio (circles stay circles); FALSE = stretched to the panel
fig7_taxon_size_by       <- "edges"  # taxon node size (area) proportional to: "edges" = stable taxon-gene edges
                                     #   drawn; "strict" = of these, edges also stable at the strict level;
                                     #   "none" = all taxa the same size
fig7_taxon_size_range    <- c(4, 16) # node size of a taxon with 0 edges / with the most edges (same scale in
                                     #   every figure, so sizes are comparable between networks)
fig7_arrangement         <- "separate" # "separate" = one file per status (network_<status>.png);
                                     # "side" = panels side by side; "stack" = one above the other (one legend)
fig7_legend_categories   <- "own"    # "separate" only: "own" = each figure lists only its own categories;
                                     # "all" = both figures list every category shown in either
                                     # (rho colour and |rho| width scales are always the same in all figures)
fig7_size_in             <- NULL     # NULL = automatic (side 19 x 9, stack 12 x 16; separate: 11 in wide, height
                                     # fitted to the network's shape, never less than the legend);
                                     # or e.g. c(width = 19, height = 9)
fig7_show_title          <- FALSE # figure title listing the display parameters; FALSE = no title
                                  # (the parameters then go in the figure legend; the text is printed in the log)
# CSVs keep every edge drawn, with columns full_colour and labelled


############################################################
# 1 — SETUP, LOG AND HELPERS
############################################################

suppressPackageStartupMessages({
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(eulerr)
  library(gridExtra)
  library(grid)
  library(igraph)
  library(ggraph)
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
  options(width = 250, max.print = 1e6)
  cat("\n\n============================================================\n")
  cat(" ", title, "\n")
  cat("============================================================\n")
}

log_parameters <- function(names) {
  section("PARAMETERS")
  for (nm in names) {
    val <- get(nm)
    cat(sprintf("%-24s = %s\n", nm,
                if (length(val) == 0) "(none)" else if (is.numeric(val)) paste(format(val, scientific = FALSE, trim = TRUE), collapse = ", ")
                else paste(val, collapse = ", ")))
  }
}

cat("05b_figures.R\nRun:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

a_rds <- file.path(de_dir, "03a_objects.rds")
b_rds <- file.path(annot_dir, "05a_annotation.rds")
if (!file.exists(a_rds)) fail("Run 03a_host_DE.R first: ", a_rds, " not found.")
if (!file.exists(b_rds)) fail("Run 05a_gene_annotation.R first: ", b_rds, " not found.")
a <- readRDS(a_rds)
b <- readRDS(b_rds)
fdr_cutoff      <- a$params$fdr_cutoff
lfc_cutoff      <- a$params$lfc_cutoff
lfc_cutoff_high <- if (is.null(a$params$lfc_cutoff_high)) 2 else a$params$lfc_cutoff_high
cat("Inputs:", a_rds, "(created", format(a$created), ") |", b_rds, "(created", format(b$created), ")\n")
cat("From 03a: fdr_cutoff =", fdr_cutoff, "| lfc_cutoff =", lfc_cutoff, "| lfc_cutoff_high =", lfc_cutoff_high,
    "| From 05a: annot_resp_top_n =", b$params$annot_resp_top_n, "\n")
if (b$created < a$created) cat("WARNING: 05a output is older than 03a output — re-run 05a.\n")
# 04a network files must be newer than the 03a output, and 05a (which annotates the network
# genes) newer than the 04a files; run order 03a -> 04a -> 05a -> 05b
net_files <- file.path(net_dir, paste0("network_", fig7_statuses, "_stable_bootstrap_edges.csv"))
net_files <- net_files[file.exists(net_files)]
if (length(net_files) > 0) {
  net_mtime <- file.mtime(net_files)
  cat("04a network files: modified", paste(format(net_mtime), collapse = ", "), "\n")
  if (min(net_mtime) < a$created) cat("WARNING: 04a output is older than 03a output — re-run 04a, then 05a.\n")
  if (b$created < max(net_mtime)) cat("WARNING: 05a output is older than 04a output — re-run 05a.\n")
}

log_parameters(c("status_levels", "comparison_order",
                 "fig6a_comparison", "fig6a_lfc_cutoff", "fig6a_n_taxa",
                 "fig6b_comparison", "fig6b_n_genes", "fig6b_lfc_min", "fig6b_row_order", "fig6b_max_age",
                 "fig6b_contrast_statuses", "fig6b_show_lfc", "fig6b_category_bar",
                 "fig6c_comparisons", "fig6c_lfc_cutoff", "fig6c_max_per_direction", "fig6c_samples",
                 "fig6c_max_height_px", "top_table_n", "fig7_statuses", "fig7_min_abs_rho", "fig7_max_genes_per_taxon",
                 "fig7_highlight_rule", "fig7_highlight_min_abs_rho", "fig7_highlight_per_taxon",
                 "fig7_highlight_exclude", "fig7_highlight_DE_only", "fig7_label_per_taxon",
                 "fig7_faded_cols", "fig7_layout", "fig7_seed", "fig7_pack_gap", "fig7_mirror",
                 "fig7_rotate_largest", "fig7_target_aspect", "fig7_bold_edge_weight", "fig7_bold_min_genes", "fig7_bold_arc", "fig7_bold_ring", "fig7_expand", "fig7_fixed_aspect",
                 "fig7_taxon_size_by", "fig7_taxon_size_range", "fig7_arrangement", "fig7_legend_categories", "fig7_size_in", "fig7_show_title"))

if (!fig6b_row_order %in% c("logFC", "cluster")) fail("fig6b_row_order must be 'logFC' or 'cluster'.")
if (!all(fig6b_contrast_statuses %in% status_levels)) fail("fig6b_contrast_statuses must be in status_levels.")
if (!fig7_highlight_rule %in% c("strict", "rho", "top")) fail("fig7_highlight_rule must be 'strict', 'rho' or 'top'.")
if (!fig7_arrangement %in% c("side", "stack", "separate")) fail("fig7_arrangement must be 'side', 'stack' or 'separate'.")
if (!fig7_legend_categories %in% c("own", "all")) fail("fig7_legend_categories must be 'own' or 'all'.")
if (!fig7_taxon_size_by %in% c("edges", "strict", "none")) fail("fig7_taxon_size_by must be 'edges', 'strict' or 'none'.")
if (!fig6c_samples %in% c("all", "DE")) fail("fig6c_samples must be 'all' or 'DE'.")
if (!is.numeric(fig7_bold_edge_weight) || fig7_bold_edge_weight <= 0) fail("fig7_bold_edge_weight must be > 0.")
if (fig6b_lfc_min < lfc_cutoff && fig6b_n_genes > b$params$annot_resp_top_n)
  fail("fig6b_lfc_min < lfc_cutoff: genes below the DE threshold are annotated by 05a only up to annot_resp_top_n (",
       b$params$annot_resp_top_n, "); increase annot_resp_top_n and re-run 05a.")
missing_cmp <- setdiff(c(fig6a_comparison, fig6b_comparison, fig6c_comparisons), unique(a$DE_all$comparison))
if (length(missing_cmp) > 0) fail("Contrasts not in the 03a results: ", paste(missing_cmp, collapse = ", "))
extra_cmp <- setdiff(unique(a$DE_all$comparison), comparison_order)
if (length(extra_cmp) > 0) fail("03a results contain contrasts not in comparison_order (old 03a output?): ",
                                paste(extra_cmp, collapse = ", "), " — re-run 03a and 05a.")

DE_all       <- a$DE_all
logcpm_disp  <- a$logcpm_disp
meta_seq     <- a$meta_seq
meta_host    <- a$meta_host
tested_genes <- unique(DE_all$gene)     # genes are filtered per model in 03a
pathogen_classified <- a$pathogen_classified
gene_annot   <- b$gene_annot

annotate <- function(df) {
  df %>%
    mutate(gene_clean = str_remove(gene, "^gene-")) %>%
    left_join(gene_annot %>% dplyr::select(gene_clean, gene_label, category), by = "gene_clean") %>%
    mutate(annotated = !is.na(category),
           gene_label = coalesce(gene_label, gene_clean),
           category   = coalesce(category, "Other/Uncharacterised"))
}

DE_taxon <- DE_all %>% filter(FDR < fdr_cutoff, pathogen != "none")


############################################################
# 2 — COLOURS
############################################################

category_palette <- c(
  "Interferon/Antiviral"          = "#EE3377",
  "Adaptive/Leukocyte"            = "#882255",
  "Innate Sensing/Inflammasome"   = "#33BBEE",
  "Cytokine/Inflammation"         = "#CC3311",
  "Innate Humoral Defense"        = "#AA3377",
  "Autophagy/Lysosomal Defense"   = "#44AA99",
  "Ciliary/Mucociliary"           = "#0077BB",
  "Epithelial Barrier"            = "#EE7733",
  "ECM/Matrix/Proteolysis"        = "#009988",
  "Apoptosis/Stress/Remodeling"   = "#999933",
  "Ion/Solute Transport"          = "#BBBBBB",
  "Metabolism"                    = "#332288",
  "Transcription/RNA/Translation" = "#66CCEE",
  "Neuronal/Sensory/Reproductive" = "#000000",
  "Coagulation/Haemostasis"       = "#661100",
  "Signaling/Trafficking/Cytoskeleton/Adhesion"         = "#228833",
  "Innate Immunity (general)"     = "#99DDFF",
  "Immune Response (general)"     = "#EE99AA",
  "Other/Uncharacterised"         = "#F2F2F2"
)
missing_cols <- setdiff(unique(gene_annot$category), names(category_palette))
if (length(missing_cols) > 0) fail("No colour for categories: ", paste(missing_cols, collapse = ", "))

status_cols <- c(healthy = "#2CA02C", unspecific = "#FF7F0E", respiratory = "#1F77B4",
                 scour = "#9467BD", respiratory_scours = "#8C564B")
status_labels <- c(healthy = "healthy", unspecific = "unspecific", respiratory = "respiratory",
                   scour = "diarrhoea", respiratory_scours = "respiratory & diarrhoea")
detection_cols <- c(negative = "#D9D9D9", positive = "#252525")
zscore_col_fun <- colorRamp2(c(-3, 0, 3), c("#313695", "#FFFFFF", "#A50026"))
# age bar in Figure 6C (days; calves without a recorded age are grey)
age_col_fun <- colorRamp2(c(0, 30, 60), c("#F7FCF5", "#74C476", "#00441B"))
age_na_col  <- "grey80"


############################################################
# 3 — FIGURE 6A: EULER DIAGRAMS
############################################################

section(paste0("FIGURE 6A: EULER DIAGRAMS (", fig6a_comparison, "; FDR < ", fdr_cutoff,
               ", |log2FC| > ", fig6a_lfc_cutoff, ")"))

DE_6a <- DE_taxon %>% filter(comparison == fig6a_comparison)

build_gene_sets <- function(df) {
  df %>% distinct(taxon, gene) %>% group_by(taxon) %>%
    summarise(genes = list(gene), .groups = "drop") %>% deframe()
}
region_counts <- function(sets) {
  genes <- unique(unlist(sets))
  if (length(genes) == 0) return(tibble(region = character(), genes = integer()))
  memb  <- vapply(sets, function(s) genes %in% s, logical(length(genes)))
  memb  <- matrix(memb, nrow = length(genes), dimnames = list(NULL, names(sets)))
  tibble(region = apply(memb, 1, function(r) paste(names(sets)[r], collapse = " & "))) %>%
    count(region, name = "genes", sort = TRUE)
}
unique_vs_all <- function(sets) {
  tibble(taxon = names(sets), genes = lengths(sets),
         not_in_any_other_taxon = vapply(names(sets), function(t)
           length(setdiff(sets[[t]], unlist(sets[names(sets) != t]))), integer(1)))
}
top_sets <- function(sets, n) sets[order(-lengths(sets), names(sets))][seq_len(min(n, length(sets)))]

sets_up_all   <- build_gene_sets(DE_6a %>% filter(logFC >  fig6a_lfc_cutoff))
sets_down_all <- build_gene_sets(DE_6a %>% filter(logFC < -fig6a_lfc_cutoff))
sets_up   <- top_sets(sets_up_all, fig6a_n_taxa)
sets_down <- top_sets(sets_down_all, fig6a_n_taxa)

cat("Genes per taxon (all taxa), up and down:\n")
full_join(unique_vs_all(sets_up_all)   %>% dplyr::rename(up = genes, up_not_shared = not_in_any_other_taxon),
          unique_vs_all(sets_down_all) %>% dplyr::rename(down = genes, down_not_shared = not_in_any_other_taxon),
          by = "taxon") %>%
  arrange(desc(coalesce(up, 0L) + coalesce(down, 0L))) %>% as.data.frame() %>% print(row.names = FALSE)
cat("\nTaxa shown, upregulated:  ", paste(names(sets_up), collapse = ", "), "\n")
cat("Taxa shown, downregulated:", paste(names(sets_down), collapse = ", "), "\n")
dropped_up   <- setdiff(names(sets_up_all), names(sets_up))
dropped_down <- setdiff(names(sets_down_all), names(sets_down))
if (length(c(dropped_up, dropped_down)) > 0)
  cat("NOT shown (fig6a_n_taxa =", fig6a_n_taxa, "): up:", paste(dropped_up, collapse = ", "),
      "| down:", paste(dropped_down, collapse = ", "), "\n")
cat("\nUpregulated: genes per exclusive region (numbers in the diagram):\n")
print(as.data.frame(region_counts(sets_up)), row.names = FALSE)
cat("\nDownregulated: genes per exclusive region (numbers in the diagram):\n")
print(as.data.frame(region_counts(sets_down)), row.names = FALSE)

taxa_6a <- union(names(sets_up), names(sets_down))
base_cols  <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
                "#CC79A7", "#999999", "#882255", "#44AA99", "#117733", "#332288")
taxon_cols <- setNames(if (length(taxa_6a) <= length(base_cols)) base_cols[seq_along(taxa_6a)]
                       else grDevices::hcl.colors(length(taxa_6a), "Dark 3"), taxa_6a)
plot_euler <- function(sets, title) {
  if (length(sets) == 0) return(textGrob(paste(title, "\nno genes"), gp = gpar(fontsize = 14)))
  fit <- euler(sets, shape = "ellipse")
  # Euler diagrams with ellipses are approximate for > 3 sets: region areas (not the
  # printed numbers, which are exact counts) may be distorted; report the fit quality
  cat(title, "— Euler fit: stress =", signif(fit$stress, 3), "| diagError =", signif(fit$diagError, 3),
      "(diagError = largest difference, in proportion of total, between a region's area and its gene count)\n")
  p <- plot(fit,
       fills      = list(fill = unname(taxon_cols[names(sets)]), alpha = 0.65),
       edges      = list(col = "grey30", lwd = 1.2),
       labels     = list(font = 2, cex = 0.95, col = "grey10"),
       quantities = list(cex = 0.85, font = 1, col = "grey20"),
       main       = list(label = title, font = 2, cex = 1.1))
  # draw inside a slightly smaller viewport so labels at the edges are not cut off
  grid::gTree(children = grid::gList(p), vp = grid::viewport(width = 0.88, height = 0.88))
}
set.seed(1)
p_up   <- plot_euler(sets_up,   paste0("Upregulated (log2FC > ", fig6a_lfc_cutoff, ")"))
p_down <- plot_euler(sets_down, paste0("Downregulated (log2FC < -", fig6a_lfc_cutoff, ")"))
png(file.path(out_dir, "euler_up_down.png"), width = 14, height = 7, units = "in", res = 300)
grid.arrange(p_up, p_down, ncol = 2,
             top = textGrob("Differentially expressed genes across taxa", gp = gpar(fontsize = 14, fontface = "bold")))
invisible(dev.off())
cat("\nSaved: euler_up_down.png\n")


############################################################
# 3b — FIGURE 6A (UpSet): DE GENES SHARED BETWEEN TAXA, UP AND DOWN
# Same gene sets as the Euler diagrams (sets_up, sets_down: fig6a_comparison,
# FDR < fdr_cutoff, |log2FC| > fig6a_lfc_cutoff). Bars = genes in each exclusive
# intersection ("distinct" mode: the region counts printed above); the dots below
# show which taxa share them; right bars = all DE genes per taxon.
############################################################

# Layout and font sizes. Font sizes are in points AT THE SAVED SIZE (printed in the log). Placed at
# full page width (~7 in, 180 mm) in a multipanel figure, the PNG is reduced by (saved width / 7):
# e.g. a 19-in-wide file is reduced ~2.7x, so 20 pt prints at ~7.5 pt.
fig6a_upset_max_comb    <- 25    # intersections shown per panel, largest first (Inf = all)
fig6a_upset_font        <- 20    # set names, axis titles, panel titles (pt); numbers and axis labels 4 pt smaller
fig6a_upset_col_in      <- 0.20  # width of each intersection column (in): smaller = bars closer together
fig6a_upset_bar_width   <- 0.85  # bar width as a fraction of its column (1 = touching; ComplexHeatmap default 0.6)
fig6a_upset_fixed_in    <- 4.5   # width per panel for set names and set-size bars (in); increase if anything is cut
fig6a_upset_height_in   <- 9     # figure height (in)
fig6a_upset_numbers_rot <- 90    # numbers above the bars: 90 = vertical (fit narrow columns), 0 = horizontal

section(paste0("FIGURE 6A (UpSet): ", fig6a_comparison, "; FDR < ", fdr_cutoff, ", |log2FC| > ", fig6a_lfc_cutoff))

# bacterial genera in italics (03a pathogen IDs start with "g__"); viruses upright
italic_taxa <- DE_all %>% filter(pathogen != "none", str_detect(pathogen, "^g__")) %>% distinct(taxon) %>% pull(taxon)

plot_upset <- function(sets, title, fill) {
  fs <- fig6a_upset_font; fs_small <- fs - 4
  if (length(sets) == 0) return(list(grob = textGrob(paste(title, "\nno genes"), gp = gpar(fontsize = fs)), n_comb = 1))
  m <- make_comb_mat(sets, mode = "distinct")
  set_total <- set_size(m)                  # all DE genes per taxon (right-hand bars)
  n_comb <- length(comb_size(m))
  if (is.finite(fig6a_upset_max_comb) && n_comb > fig6a_upset_max_comb) {
    keep   <- order(comb_size(m), decreasing = TRUE)[seq_len(fig6a_upset_max_comb)]
    hidden <- setdiff(seq_len(n_comb), keep)
    cat(title, "— intersections not shown (smallest):", length(hidden), "with",
        sum(comb_size(m)[hidden]), "genes in total (all listed in the region counts above)\n")
    # taxa with no gene in the shown intersections are left out (their rows would be empty);
    # the matrix is rebuilt from the genes in the shown intersections only, so these keep
    # exactly the same membership and sizes
    kept_codes  <- comb_name(m)[keep]
    in_shown    <- rowSums(sapply(kept_codes, function(cd) as.integer(strsplit(cd, "")[[1]]))) > 0
    sets_shown  <- set_name(m)[in_shown]
    genes_shown <- unlist(lapply(kept_codes, function(cd) extract_comb(m, cd)))
    if (any(!in_shown))
      cat(title, "— taxa whose DE genes are only in intersections not shown (rows left out):",
          paste0(set_name(m)[!in_shown], " (", set_total[!in_shown], " genes)", collapse = ", "), "\n")
    m <- make_comb_mat(lapply(sets[sets_shown], intersect, genes_shown), mode = "distinct")
    # check: same intersections and the same genes as before the rebuild
    stopifnot(setequal(comb_name(m), vapply(kept_codes, function(cd)
      paste(strsplit(cd, "")[[1]][in_shown], collapse = ""), character(1))),
      sum(comb_size(m)) == length(genes_shown))
  }
  set_total <- set_total[set_name(m)]
  cs <- comb_size(m)
  ht <- UpSet(m,
              set_order  = order(set_total, decreasing = TRUE),
              comb_order = order(cs, decreasing = TRUE),
              comb_col = fill, pt_size = unit(fs * 0.22, "mm"), lwd = fs / 8,
              bg_col = c("#F2F2F2", "#FFFFFF"), bg_pt_col = "#D9D9D9",
              row_names_gp = gpar(fontsize = fs, fontface = if_else(set_name(m) %in% italic_taxa, "italic", "plain")),
              column_title = title, column_title_gp = gpar(fontsize = fs + 2, fontface = "bold"),
              top_annotation = upset_top_annotation(m, add_numbers = TRUE, numbers_rot = fig6a_upset_numbers_rot,
                                                    numbers_gp = gpar(fontsize = fs_small), bar_width = fig6a_upset_bar_width,
                                                    gp = gpar(fill = fill, col = NA),
                                                    annotation_name_rot = 90, annotation_name_gp = gpar(fontsize = fs),
                                                    axis_param = list(gp = gpar(fontsize = fs_small)),
                                                    height = unit(fs * 0.3, "cm"), ylim = c(0, max(cs) * 1.3)),
              # right-hand bars: ALL DE genes per taxon (including those in intersections not shown)
              right_annotation = rowAnnotation(
                "Set size" = anno_barplot(set_total, add_numbers = TRUE, numbers_gp = gpar(fontsize = fs_small),
                                          bar_width = fig6a_upset_bar_width, border = FALSE,
                                          gp = gpar(fill = "grey45", col = NA),
                                          axis_param = list(gp = gpar(fontsize = fs_small)),
                                          ylim = c(0, max(set_total) * 1.35), width = unit(fs * 0.2, "cm")),
                annotation_name_side = "bottom", annotation_name_rot = 0, annotation_name_gp = gpar(fontsize = fs)))
  list(grob = grid.grabExpr(draw(ht, padding = unit(c(4, 4, 4, 4), "mm"))), n_comb = length(cs))
}

u_up   <- plot_upset(sets_up,   paste0("Upregulated (log2FC > ", fig6a_lfc_cutoff, ")"),    "#B2182B")
u_down <- plot_upset(sets_down, paste0("Downregulated (log2FC < -", fig6a_lfc_cutoff, ")"), "#2166AC")

# panel widths follow the number of intersections shown (+ room for names and set-size bars)
w_up   <- u_up$n_comb   * fig6a_upset_col_in + fig6a_upset_fixed_in
w_down <- u_down$n_comb * fig6a_upset_col_in + fig6a_upset_fixed_in
png(file.path(out_dir, "upset_up_down.png"), width = w_up + w_down, height = fig6a_upset_height_in, units = "in", res = 300)
grid.arrange(u_up$grob, u_down$grob, ncol = 2, widths = c(w_up, w_down),
             top = textGrob(paste0("Differentially expressed genes shared between taxa (", fig6a_comparison, ")"),
                            gp = gpar(fontsize = fig6a_upset_font + 4, fontface = "bold")))
invisible(dev.off())
cat("Saved: upset_up_down.png (", round(w_up + w_down, 1), "x", fig6a_upset_height_in, "in; at 7 in wide, fonts print at",
    round(fig6a_upset_font * 7 / (w_up + w_down), 1), "pt )\n")



############################################################
# 4 — FIGURE 6B: RESPIRATORY VS HEALTHY, LARGEST |log2FC|
############################################################

section(paste0("FIGURE 6B: ", if (is.finite(fig6b_n_genes)) paste(fig6b_n_genes, "") else "ALL ",
               "DE GENES IN ", fig6b_comparison, " (FDR < ", fdr_cutoff, ", |log2FC| > ", fig6b_lfc_min,
               if (is.finite(fig6b_n_genes)) "; largest |log2FC|" else "", ")"))

cand_6b <- DE_all %>%
  filter(comparison == fig6b_comparison, FDR < fdr_cutoff, abs(logFC) > fig6b_lfc_min)
top_genes_6b <- cand_6b %>%
  annotate() %>%
  filter(gene %in% rownames(logcpm_disp)) %>%
  arrange(desc(abs(logFC)))
if (is.finite(fig6b_n_genes)) top_genes_6b <- top_genes_6b %>% slice_head(n = fig6b_n_genes)
top_genes_6b <- top_genes_6b %>%
  arrange(desc(logFC)) %>%                       # most upregulated first, most downregulated last
  transmute(gene, gene_clean, gene_label, category, annotated, logFC, FDR)
if (nrow(top_genes_6b) == 0) fail("No genes for Figure 6B (genes meeting the criteria: ", nrow(cand_6b), ").")
cat("Genes meeting the criteria:", nrow(cand_6b), "| shown:", nrow(top_genes_6b),
    "| up:", sum(top_genes_6b$logFC > 0), "| down:", sum(top_genes_6b$logFC < 0),
    "| smallest |log2FC| shown:", round(min(abs(top_genes_6b$logFC)), 2), "\n")
if (any(!top_genes_6b$annotated)) fail("Figure 6B genes without annotation (re-run 05a): ",
                                      paste(top_genes_6b$gene_clean[!top_genes_6b$annotated], collapse = ", "))

# columns: the contrasted statuses side by side, then all other statuses pooled; age order within
status_order_6b <- c(fig6b_contrast_statuses, setdiff(status_levels, fig6b_contrast_statuses))
meta_6b <- meta_seq %>%
  filter(Sample_id %in% colnames(logcpm_disp), Status %in% status_levels,
         is.na(Age_days) | Age_days <= fig6b_max_age) %>%
  mutate(Status = droplevels(factor(Status, levels = status_order_6b)),
         block  = if_else(Status %in% fig6b_contrast_statuses, as.character(Status), "other")) %>%
  arrange(Status, is.na(Age_days), Age_days) %>%
  mutate(block = droplevels(factor(block, levels = c(fig6b_contrast_statuses, "other"))))
cat("Calves shown:", nrow(meta_6b), "| with unknown age:", sum(is.na(meta_6b$Age_days)), "\n")
print(table(Status = meta_6b$Status, useNA = "ifany"))

expr_6b <- logcpm_disp[top_genes_6b$gene, meta_6b$Sample_id, drop = FALSE]
keep_6b <- apply(expr_6b, 1, var) > 0
if (any(!keep_6b)) cat("Genes dropped (no variation):", paste(top_genes_6b$gene[!keep_6b], collapse = ", "), "\n")
expr_6b <- expr_6b[keep_6b, , drop = FALSE]
top_genes_6b <- top_genes_6b[keep_6b, ]
row_labels_6b <- make.unique(top_genes_6b$gene_label)
rownames(expr_6b) <- row_labels_6b
z_6b <- t(scale(t(expr_6b)))                      # z-score per gene across all calves shown

# column blocks and titles (n calves); wider gap before the pooled other statuses
block_n_6b      <- table(meta_6b$block)
block_titles_6b <- paste0(if_else(names(block_n_6b) == "other", "Other statuses",
                                  str_to_sentence(status_labels[names(block_n_6b)])),
                          "\n(n = ", as.integer(block_n_6b), ")")
col_gap_6b <- if (length(block_n_6b) < 2) unit(1, "mm") else
  unit(if_else(names(block_n_6b)[-1] == "other", 6, 1.5), "mm")

# row blocks: up- and downregulated in the test group of the contrast
g1_6b <- cand_6b$group_1_name[1]
g2_6b <- cand_6b$group_2_name[1]
dir_6b <- factor(if_else(top_genes_6b$logFC > 0, "up", "down"), levels = c("up", "down"))
dir_n_6b <- table(dir_6b)
dir_titles_6b <- c(up   = paste0("Higher in ", g2_6b, "\n(", dir_n_6b[["up"]], " genes)"),
                   down = paste0("Lower in ", g2_6b, "\n(", dir_n_6b[["down"]], " genes)"))[names(dir_n_6b)[dir_n_6b > 0]]

statuses_6b <- levels(meta_6b$Status)
ha_top_6b <- HeatmapAnnotation(
  Status = meta_6b$Status, Age = meta_6b$Age_days,
  col = list(Status = status_cols[statuses_6b], Age = age_col_fun),
  na_col = age_na_col,
  annotation_label = c("Status", "Age (days)"),
  annotation_name_gp = gpar(fontsize = 10),
  simple_anno_size = unit(4, "mm"), gap = unit(1, "mm"),
  annotation_legend_param = list(
    Status = list(title = "Status", at = statuses_6b, labels = status_labels[statuses_6b]),
    Age    = list(title = "Age (days)", at = c(0, 20, 40, 60))))
categories_6b <- intersect(names(category_palette), unique(top_genes_6b$category))
ha_left_6b <- rowAnnotation(
  Category = factor(top_genes_6b$category, levels = categories_6b),
  col = list(Category = category_palette[categories_6b]),
  show_annotation_name = FALSE, simple_anno_size = unit(4, "mm"),
  annotation_legend_param = list(Category = list(title = "Category")))
ha_right_6b <- if (fig6b_show_lfc) rowAnnotation(
  log2FC = anno_barplot(top_genes_6b$logFC, baseline = 0, bar_width = 0.8,
                        gp = gpar(fill = if_else(top_genes_6b$logFC > 0, "#B2182B", "#2166AC"), col = NA),
                        axis_param = list(gp = gpar(fontsize = 7), labels_rot = 0), width = unit(2.2, "cm")),
  annotation_label = "log2FC", annotation_name_gp = gpar(fontsize = 9), annotation_name_rot = 0) else NULL

ht_6b <- Heatmap(
  z_6b, name = "z-score", col = zscore_col_fun,
  heatmap_legend_param = list(title = "z-score", at = -3:3),
  cluster_columns = FALSE,
  cluster_rows = fig6b_row_order == "cluster", cluster_row_slices = FALSE,
  clustering_distance_rows = "pearson", show_row_dend = FALSE,
  row_split = dir_6b, row_title = unname(dir_titles_6b), row_gap = unit(3, "mm"),
  row_title_gp = gpar(fontsize = 11, fontface = "bold"),
  column_split = meta_6b$block, column_title = block_titles_6b, column_gap = col_gap_6b,
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  show_column_names = FALSE, row_names_gp = gpar(fontsize = 8), row_names_side = "right",
  top_annotation = ha_top_6b, left_annotation = ha_left_6b, right_annotation = ha_right_6b,
  use_raster = FALSE)
extra_legends_6b <- if (any(is.na(meta_6b$Age_days)))
  list(Legend(title = "Age unknown", at = "", legend_gp = gpar(fill = age_na_col))) else list()

fig6b_file  <- "heatmap_6B_DE_genes.png"
fig6b_title <- paste0(str_to_sentence(g2_6b), " vs ", g1_6b, ": ", nrow(z_6b),
                      " DE genes (FDR < ", fdr_cutoff, ", |log2FC| > ", fig6b_lfc_min, ")")
fig6b_h <- max(8, 2.5 + nrow(z_6b) * 0.16)          # inches

# genes up / down per category (the genes in the heatmap)
cat_counts_6b <- top_genes_6b %>%
  mutate(direction = factor(if_else(logFC > 0, "up", "down"), levels = c("up", "down"))) %>%
  count(category, direction) %>%
  complete(category, direction, fill = list(n = 0L)) %>%
  group_by(category) %>% mutate(total = sum(n)) %>% ungroup()
cat("\nFigure 6B genes per category (up / down):\n")
cat_counts_6b %>% pivot_wider(names_from = direction, values_from = n) %>% arrange(desc(total)) %>%
  as.data.frame() %>% print(row.names = FALSE)

if (fig6b_category_bar) {
  # categories ordered by total (largest on top), Other/Uncharacterised last
  cat_levels_6b <- cat_counts_6b %>% distinct(category, total) %>%
    arrange(category == "Other/Uncharacterised", desc(total), category) %>% pull(category)
  cat_counts_6b <- cat_counts_6b %>% mutate(category = factor(category, levels = rev(cat_levels_6b)),
                                            x = if_else(direction == "up", n, -n))
  max_n   <- max(cat_counts_6b$n)
  key_x   <- -max_n * 1.45                              # category colour squares, left of the bars
  p_cat_6b <- ggplot(cat_counts_6b, aes(x = x, y = category)) +
    geom_col(aes(fill = direction), width = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_text(data = ~ filter(.x, n > 0),
              aes(label = n, hjust = if_else(direction == "up", -0.3, 1.3)), size = 2.8) +
    geom_point(data = ~ distinct(.x, category), aes(x = key_x, colour = as.character(category)),
               shape = 15, size = 4.2) +
    scale_fill_manual(values = c(up = "#B2182B", down = "#2166AC"),
                      labels = c(up = paste("Higher in", g2_6b), down = paste("Lower in", g2_6b)), name = NULL) +
    scale_colour_manual(values = category_palette, guide = "none") +
    scale_x_continuous(labels = abs, limits = c(key_x * 1.08, max_n * 1.35),
                       breaks = function(l) { b <- pretty(c(-max_n, max_n), 4); b[abs(b) <= max_n] }) +
    labs(x = "DE genes", y = NULL, title = "DE genes per category") +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 11),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = 8, colour = "black"),
          legend.position = "bottom", legend.text = element_text(size = 8),
          legend.key.size = unit(0.35, "cm"), plot.margin = margin(2, 6, 2, 2))
  bar_h_6b <- 1.3 + 0.24 * length(cat_levels_6b)       # inches

  # legends (the category key is the barplot)
  lg_6b <- c(
    list(Legend(title = "z-score", col_fun = zscore_col_fun, at = -3:3),
         Legend(title = "Status", at = statuses_6b, labels = status_labels[statuses_6b],
                legend_gp = gpar(fill = status_cols[statuses_6b])),
         Legend(title = "Age (days)", col_fun = age_col_fun, at = c(0, 20, 40, 60))),
    extra_legends_6b)
  right_w  <- 4.6                                       # inches, right-hand column
  legend_top_6b <- fig6b_h - 0.35 - bar_h_6b - 0.3
  legends_6b <- packLegend(list = lg_6b, direction = "vertical", gap = unit(5, "mm"),
                           max_height = unit(legend_top_6b - 0.2, "in"))

  png(file.path(out_dir, fig6b_file), width = 10.6 + right_w, height = fig6b_h, units = "in", res = 300)
  grid.newpage()
  pushViewport(viewport(x = 0, y = 0, width = unit(10.6, "in"), height = 1, just = c("left", "bottom")))
  draw(ht_6b, newpage = FALSE, show_heatmap_legend = FALSE, show_annotation_legend = FALSE,
       column_title = fig6b_title, column_title_gp = gpar(fontsize = 13, fontface = "bold"))
  popViewport()
  pushViewport(viewport(x = unit(10.6, "in"), y = unit(fig6b_h - 0.35, "in"),
                        width = unit(right_w - 0.15, "in"), height = unit(bar_h_6b, "in"), just = c("left", "top")))
  grid.draw(ggplotGrob(p_cat_6b))
  popViewport()
  draw(legends_6b, x = unit(10.6 + 0.25, "in"), y = unit(legend_top_6b, "in"), just = c("left", "top"))
  invisible(dev.off())
} else {
  png(file.path(out_dir, fig6b_file), width = 14, height = fig6b_h, units = "in", res = 300)
  draw(ht_6b, merge_legends = TRUE, annotation_legend_list = extra_legends_6b,
       column_title = fig6b_title, column_title_gp = gpar(fontsize = 13, fontface = "bold"))
  invisible(dev.off())
}
cat("Saved:", fig6b_file, "\n")
write_csv(top_genes_6b %>% dplyr::select(gene = gene_clean, gene_label, category, logFC, FDR),
          file.path(out_dir, "heatmap_6B_DE_genes.csv"))

# per gene: mean log2 CPM by status, and Spearman correlation with age
expr_stats_6b <- as.data.frame(expr_6b) %>%
  rownames_to_column("gene_label") %>%
  pivot_longer(-gene_label, names_to = "Sample_id", values_to = "logCPM") %>%
  left_join(meta_6b %>% dplyr::select(Sample_id, Status, Age_days), by = "Sample_id")
age_rho <- expr_stats_6b %>% group_by(gene_label) %>%
  summarise(rho_age = round(cor(logCPM, Age_days, method = "spearman", use = "complete.obs"), 2), .groups = "drop")
status_means <- expr_stats_6b %>% group_by(gene_label, Status) %>%
  summarise(mean_logCPM = round(mean(logCPM), 1), .groups = "drop") %>%
  pivot_wider(names_from = Status, values_from = mean_logCPM)
cat("\nGenes in Figure 6B, in figure order (mean TMM log2 CPM by status; Spearman rho with age):\n")
top_genes_6b %>%
  transmute(gene_label = row_labels_6b, gene = gene_clean, category,
            logFC = round(logFC, 2), FDR = signif(FDR, 2)) %>%
  left_join(status_means, by = "gene_label") %>%
  left_join(age_rho, by = "gene_label") %>%
  as.data.frame() %>% print(row.names = FALSE)


############################################################
# 5 — FIGURE 6C: DE GENES PER TAXON (|log2FC| > fig6c_lfc_cutoff)
############################################################

section(paste0("FIGURE 6C: GENES WITH FDR < ", fdr_cutoff, " AND |log2FC| > ", fig6c_lfc_cutoff,
               " IN ", paste(fig6c_comparisons, collapse = " OR ")))

DE_6c <- DE_taxon %>% filter(comparison %in% fig6c_comparisons)
fig6c_sel <- DE_6c %>% filter(abs(logFC) > fig6c_lfc_cutoff)
if (is.finite(fig6c_max_per_direction)) {
  fig6c_sel <- fig6c_sel %>%
    group_by(taxon, comparison, up = logFC > 0) %>%
    slice_max(abs(logFC), n = fig6c_max_per_direction, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(-up)
  rule_text <- paste0("up to ", fig6c_max_per_direction, " genes per taxon, contrast and direction (largest |log2FC|)")
} else {
  rule_text <- "all such genes"
}
cat("Selection:", rule_text, "| calves:", fig6c_samples, "\n")
if (nrow(fig6c_sel) == 0) fail("No genes for Figure 6C — check fig6c_lfc_cutoff and fig6c_comparisons.")

cat("\nGenes selected per taxon and contrast:\n")
fig6c_sel %>%
  count(taxon, comparison, direction = if_else(logFC > 0, "up", "down")) %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  as.data.frame() %>% print(row.names = FALSE)

# one row per taxon and gene; selected_in = the contrast(s) the gene was selected from
fig6c_list <- fig6c_sel %>%
  group_by(taxon, gene) %>%
  summarise(selected_in = paste(sort(unique(comparison)), collapse = " + "),
            max_logFC = logFC[which.max(abs(logFC))], .groups = "drop") %>%
  arrange(taxon, desc(max_logFC)) %>%
  transmute(gene_host = gene, gene = str_remove(gene, "^gene-"), taxon, selected_in)
cat("\nRows per panel (distinct genes per taxon), by the contrast(s) they were selected from:\n")
fig6c_list %>% count(taxon, selected_in) %>% as.data.frame() %>% print(row.names = FALSE)

# DE support for each gene in its taxon (fig6c_comparisons, FDR < fdr_cutoff, |log2FC| > lfc_cutoff)
fig6c_support <- DE_6c %>%
  filter(abs(logFC) > lfc_cutoff) %>%
  mutate(gene_clean = str_remove(gene, "^gene-")) %>%
  semi_join(fig6c_list, by = c("gene_clean" = "gene", "taxon")) %>%
  arrange(FDR) %>%
  group_by(gene_clean, taxon) %>%
  summarise(DE_tests = paste0(comparison, " ", sprintf("%+.2f", logFC), " (", signif(FDR, 2), ")", collapse = "; "),
            .groups = "drop")

genes_6c <- fig6c_list %>%
  left_join(gene_annot %>% dplyr::select(gene_clean, gene_label, category), by = c("gene" = "gene_clean")) %>%
  mutate(annotated  = !is.na(category),
         gene_label = coalesce(gene_label, gene),
         category   = coalesce(category, "Other/Uncharacterised")) %>%
  left_join(fig6c_support, by = c("gene" = "gene_clean", "taxon")) %>%
  mutate(in_expression_data = gene_host %in% rownames(logcpm_disp),
         tested = gene_host %in% tested_genes)
if (any(!genes_6c$annotated)) cat("WARNING:", sum(!genes_6c$annotated),
                                  "Figure 6C genes not in the 05a annotation (re-run 05a).\n")

cat("\nGenes in Figure 6C and their DE support in the same taxon (comparison logFC (FDR)):\n")
genes_6c %>%
  transmute(line = paste0(taxon, " | ", gene_label, " | ", category, " | selected in: ", selected_in,
                          if_else(in_expression_data, "", " | NOT IN EXPRESSION DATA"),
                          "\n      ", coalesce(DE_tests, "-"))) %>%
  pull(line) %>% cat(sep = "\n")
write_csv(genes_6c, file.path(out_dir, "fig6C_genes.csv"))

meta_6c <- if (fig6c_samples == "DE") {
  meta_host %>% dplyr::select(Sample_id, Status, Age_days)
} else {
  meta_seq %>% filter(Sample_id %in% colnames(logcpm_disp)) %>% dplyr::select(Sample_id, Status, Age_days)
}

heatmap_inputs <- list()
panel_counts   <- list()
for (p in unique(genes_6c$taxon)) {
  genes_p <- genes_6c %>% filter(taxon == p, in_expression_data) %>% distinct(gene_host, .keep_all = TRUE)
  if (nrow(genes_p) == 0) next
  class_df <- pathogen_classified %>%
    filter(taxon == p, class %in% c("negative", "positive")) %>%
    dplyr::select(Sample_id, class) %>%
    inner_join(meta_6c, by = "Sample_id") %>%
    mutate(Status = factor(Status, levels = status_levels),
           class  = factor(class, levels = c("negative", "positive"))) %>%
    filter(!is.na(Status)) %>%
    arrange(class, Status, Age_days) %>%
    # column blocks: taxon-negative / positive, then status (calves ordered by age within block)
    mutate(block = droplevels(factor(paste(class, Status),
                                     levels = as.vector(t(outer(c("negative", "positive"), status_levels, paste))))))
  if (nrow(class_df) < 3) next
  expr <- logcpm_disp[genes_p$gene_host, class_df$Sample_id, drop = FALSE]
  keep <- apply(expr, 1, function(x) var(x) > 0)
  expr <- expr[keep, , drop = FALSE]
  genes_p <- genes_p[keep, ]
  if (nrow(expr) == 0) next
  expr_scaled <- t(scale(t(expr)))
  rownames(expr_scaled) <- make.unique(genes_p$gene_label)
  panel_counts[[p]] <- class_df %>% count(class, Status) %>% mutate(taxon = p)
  heatmap_inputs[[p]] <- list(expr = expr_scaled, class_df = class_df, category = genes_p$category)
}
if (length(heatmap_inputs) == 0) fail("No Figure 6C panels generated — check the gene selection parameters.")

cat("\nCalves per panel (taxon class x clinical status):\n")
bind_rows(panel_counts) %>%
  pivot_wider(names_from = Status, values_from = n, values_fill = 0) %>%
  relocate(taxon) %>% as.data.frame() %>% print(row.names = FALSE)

# row height: 40 px per gene (at res_6c), reduced so the whole figure stays within
# fig6c_max_height_px; row names hidden when rows become too small to label
res_6c      <- 180
n_rows_6c   <- vapply(heatmap_inputs, function(h) nrow(h$expr), integer(1))
overhead_px <- 330 * length(heatmap_inputs) + 60     # title, 3 annotation bars and padding per panel
row_px      <- max(2, min(40, floor((fig6c_max_height_px - overhead_px) / sum(n_rows_6c))))
row_font    <- max(5, floor(13 * row_px / 40))       # 13 pt at 40 px per row
show_names  <- row_px >= 13                          # 5 pt (the smallest font) needs >= 13 px
cat("\nFigure 6C rows:", sum(n_rows_6c), "| row height:", row_px, "px (40 = default)",
    "| row names", if (show_names) paste0("shown (font ", row_font, " pt)") else "HIDDEN (rows too small)", "\n")
if (!show_names) cat("NOTE: set fig6c_max_per_direction to a finite value to label the genes.\n")

heatmap_list <- lapply(names(heatmap_inputs), function(p) {
  hi <- heatmap_inputs[[p]]
  Heatmap(
    hi$expr, name = paste0("z_", p), col = zscore_col_fun,
    cluster_rows = nrow(hi$expr) > 1, cluster_columns = FALSE,
    show_column_names = FALSE, show_row_names = show_names,
    row_names_gp = gpar(fontsize = row_font, fontfamily = "sans"),
    row_names_max_width = unit(8, "cm"),
    top_annotation = HeatmapAnnotation(
      Status = hi$class_df$Status, Detection = hi$class_df$class, Age = hi$class_df$Age_days,
      col = list(Status = status_cols, Detection = detection_cols, Age = age_col_fun),
      na_col = age_na_col,
      annotation_name_gp = gpar(fontsize = 16, fontfamily = "sans"),
      annotation_height = unit(c(6, 6, 6), "mm"), show_legend = FALSE),
    left_annotation = rowAnnotation(
      category = hi$category,
      col = list(category = category_palette[unique(hi$category)]),
      annotation_width = unit(6, "mm"), show_annotation_name = FALSE, show_legend = FALSE),
    column_title = p,
    column_title_gp = gpar(fontsize = 22, fontface = "bold.italic", fontfamily = "sans"),
    height = unit(max(nrow(hi$expr) * row_px / res_6c, 14 / 25.4), "inches"),   # row_px per row, >= 14 mm
    column_split = hi$class_df$block,
    # small gaps between statuses, wide gap between taxon-negative and taxon-positive calves
    column_gap = unit({
      cls <- sub(" .*", "", levels(hi$class_df$block))
      if (length(cls) < 2) 1 else ifelse(head(cls, -1) != tail(cls, -1), 4, 0.8)
    }, "mm"),
    show_heatmap_legend = FALSE)
})

panel_w  <- 3500
legend_w <- 1100
total_w  <- panel_w + legend_w
panel_padding <- function(i) unit(c(6, 1, if (i == 1) 8 else 2, 1), "mm")
# panel heights measured from the drawn heatmaps (body + title + annotation bars + padding),
# so panels neither overlap nor leave blank bands
pdf(NULL)
panel_heights <- vapply(seq_along(heatmap_list), function(i) {
  h <- tryCatch({
    d <- draw(heatmap_list[[i]], show_heatmap_legend = FALSE, show_annotation_legend = FALSE,
              padding = panel_padding(i))
    convertHeight(ComplexHeatmap:::height(d), "inches", valueOnly = TRUE) * res_6c
  }, error = function(err) NA_real_)
  if (is.na(h)) h <- max(n_rows_6c[i] * row_px, 14 / 25.4 * res_6c) + 330 + if (i == 1) 45 else 0
  ceiling(h) + 10
}, numeric(1))
invisible(dev.off())
panels_h <- sum(panel_heights)

statuses_6c   <- intersect(status_levels, unique(unlist(lapply(heatmap_inputs, function(h) as.character(h$class_df$Status)))))
categories_6c <- intersect(names(category_palette), unique(unlist(lapply(heatmap_inputs, `[[`, "category"))))
legend_title_gp <- gpar(fontsize = 18, fontface = "bold")
legend_label_gp <- gpar(fontsize = 16)
legend_6c <- packLegend(
  Legend(title = "z-score", col_fun = zscore_col_fun, at = -3:3,
         title_gp = legend_title_gp, labels_gp = legend_label_gp,
         legend_height = unit(5, "cm"), grid_width = unit(8, "mm")),
  Legend(title = "Status", at = statuses_6c, labels = status_labels[statuses_6c],
         legend_gp = gpar(fill = status_cols[statuses_6c]),
         title_gp = legend_title_gp, labels_gp = legend_label_gp,
         grid_height = unit(7, "mm"), grid_width = unit(7, "mm")),
  Legend(title = "Age (days)", col_fun = age_col_fun, at = c(0, 20, 40, 60),
         title_gp = legend_title_gp, labels_gp = legend_label_gp,
         legend_height = unit(3, "cm"), grid_width = unit(8, "mm")),
  Legend(title = "Age unknown", at = "", legend_gp = gpar(fill = age_na_col),
         title_gp = legend_title_gp, labels_gp = legend_label_gp,
         grid_height = unit(7, "mm"), grid_width = unit(7, "mm")),
  Legend(title = "Detection", at = names(detection_cols), legend_gp = gpar(fill = detection_cols),
         title_gp = legend_title_gp, labels_gp = legend_label_gp,
         grid_height = unit(7, "mm"), grid_width = unit(7, "mm")),
  Legend(title = "Category", at = categories_6c, legend_gp = gpar(fill = category_palette[categories_6c]),
         title_gp = legend_title_gp, labels_gp = legend_label_gp,
         grid_height = unit(7, "mm"), grid_width = unit(7, "mm"), ncol = 1),
  direction = "vertical", gap = unit(6, "mm"))

# panels and legend are rendered to grobs off-screen, then arranged; the figure
# height is the taller of the stacked panels and the legend
pdf(NULL)
legend_h_px <- ceiling((convertHeight(grobHeight(legend_6c@grob), "inches", valueOnly = TRUE) + 0.4) * res_6c)
invisible(dev.off())
total_h <- max(panels_h, legend_h_px)
panel_grobs <- lapply(seq_along(heatmap_list), function(i) {
  grid.grabExpr(draw(heatmap_list[[i]], show_heatmap_legend = FALSE, show_annotation_legend = FALSE,
                     padding = panel_padding(i)),
                width = panel_w / res_6c, height = panel_heights[i] / res_6c)
})
legend_grob <- grid.grabExpr(draw(legend_6c, x = unit(0.05, "npc"), y = unit(0.98, "npc"), just = c("left", "top")),
                             width = legend_w / res_6c, height = total_h / res_6c)
png(file.path(out_dir, "heatmap_combined_panel.png"), width = total_w, height = total_h, res = res_6c)
grid.newpage()
grid.draw(arrangeGrob(
  arrangeGrob(grobs = c(panel_grobs, list(nullGrob())), ncol = 1,
              heights = unit(c(panel_heights, total_h - panels_h), "null")),
  legend_grob, ncol = 2, widths = unit(c(panel_w, legend_w), "null")))
invisible(dev.off())
cat("\nSaved: heatmap_combined_panel.png (", total_w, "x", total_h, "px )\n")


############################################################
# 5b — FIGURE 6C (summary): DE GENES PER FUNCTIONAL CATEGORY AND TAXON
# Bubble matrix: rows = functional categories (05a), columns = taxa, each split into
# up | down; dot size and number = DE genes (FDR < fdr_cutoff, |log2FC| > fig6c_bubble_lfc)
# in fig6c_bubble_comparison. Same genes as the gene-level 6C heatmap when that heatmap
# uses the same contrast and cut-off.
############################################################

fig6c_bubble_comparison <- "positive_vs_negative"
fig6c_bubble_lfc        <- fig6c_lfc_cutoff   # 2 = same genes as the 6C heatmap
fig6c_bubble_show_other <- TRUE               # FALSE = drop the Other/Uncharacterised row

section(paste0("FIGURE 6C (summary): DE GENES PER CATEGORY AND TAXON (", fig6c_bubble_comparison,
               "; FDR < ", fdr_cutoff, ", |log2FC| > ", fig6c_bubble_lfc, ")"))

bub <- DE_taxon %>%
  filter(comparison == fig6c_bubble_comparison, abs(logFC) > fig6c_bubble_lfc) %>%
  annotate() %>%
  mutate(direction = factor(if_else(logFC > 0, "up", "down"), levels = c("up", "down")))
if (nrow(bub) == 0) fail("No genes for the Figure 6C summary — check fig6c_bubble_comparison and fig6c_bubble_lfc.")
if (any(!bub$annotated)) cat("WARNING:", sum(!bub$annotated),
                            "genes not in the 05a annotation (counted as Other/Uncharacterised) — re-run 05a.\n")
if (!fig6c_bubble_show_other) {
  cat("Other/Uncharacterised genes left out:", sum(bub$category == "Other/Uncharacterised"), "\n")
  bub <- bub %>% filter(category != "Other/Uncharacterised")
}

# taxa ordered by number of DE genes; categories by total, Other/Uncharacterised last
taxon_tot <- bub %>%
  count(taxon, direction) %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0, names_expand = TRUE) %>%
  mutate(total = up + down) %>%
  arrange(desc(total), taxon)
cat_order <- bub %>% count(category) %>%
  arrange(category == "Other/Uncharacterised", desc(n), category) %>% pull(category)

cat("DE genes per taxon:\n")
print(as.data.frame(taxon_tot), row.names = FALSE)

# check: genes per taxon here vs in the 6C heatmap (they differ if the selection parameters
# differ, or the heatmap dropped genes without variation / taxa with < 3 calves)
chk_6c <- full_join(
  tibble(taxon = names(heatmap_inputs), heatmap = vapply(heatmap_inputs, function(h) nrow(h$expr), integer(1))),
  taxon_tot %>% transmute(taxon, bubble = as.integer(total)), by = "taxon") %>%
  mutate(across(c(heatmap, bubble), ~ coalesce(.x, 0L))) %>%
  filter(heatmap != bubble)
if (nrow(chk_6c) > 0) {
  cat("\nNOTE: genes per taxon differ between the 6C heatmap and this summary (different contrast, cut-off",
      "or per-direction cap; genes without variation or taxa with < 3 calves dropped from the heatmap;",
      "or Other/Uncharacterised left out here):\n")
  print(as.data.frame(chk_6c), row.names = FALSE)
} else cat("\nCheck: genes per taxon identical to the 6C heatmap.\n")
bub_counts <- bub %>% count(taxon, category, direction, name = "n")
bub_table <- bub_counts %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0, names_expand = TRUE) %>%
  mutate(cell = paste0(up, "/", down)) %>%
  dplyr::select(category, taxon, cell) %>%
  pivot_wider(names_from = taxon, values_from = cell, values_fill = "0/0") %>%
  dplyr::select(category, any_of(taxon_tot$taxon)) %>%
  arrange(match(category, cat_order))
cat("\nDE genes per category (up/down) and taxon:\n")
print(as.data.frame(bub_table), row.names = FALSE)
write_csv(bub_table, file.path(out_dir, "fig6C_category_counts.csv"))

# facet labels: bacterial genera in italics (03a pathogen IDs start with "g__")
italic_tax <- DE_all %>% filter(str_detect(pathogen, "^g__")) %>% distinct(taxon) %>% pull(taxon)
panel_lab  <- setNames(if_else(taxon_tot$taxon %in% italic_tax,
                               paste0("italic('", taxon_tot$taxon, "')"), paste0("'", taxon_tot$taxon, "'")),
                       taxon_tot$taxon)
panel_levels <- unname(panel_lab)

# sub-columns: "up (n)" and "down (n)" under each taxon, up always first
x_df <- taxon_tot %>%
  pivot_longer(c(up, down), names_to = "direction", values_to = "n_dir") %>%
  mutate(x = paste0(direction, "\n(", n_dir, ")"),
         panel = factor(panel_lab[taxon], levels = panel_levels))
x_levels <- c(unique(x_df$x[x_df$direction == "up"]), unique(x_df$x[x_df$direction == "down"]))
x_df <- x_df %>% mutate(x = factor(x, levels = x_levels))

bub_plot <- bub_counts %>%
  mutate(direction = as.character(direction)) %>%
  left_join(x_df %>% dplyr::select(taxon, direction, x, panel), by = c("taxon", "direction"))

p_6c_bub <- ggplot() +
  geom_blank(data = x_df, aes(x = x, y = cat_order[1])) +           # keeps empty up/down sub-columns
  geom_point(data = bub_plot, aes(x = x, y = category, size = n, fill = direction),
             shape = 21, colour = "grey25", stroke = 0.3) +
  geom_text(data = bub_plot, aes(x = x, y = category, label = n), colour = "white", size = 2.6, fontface = "bold") +
  facet_grid(. ~ panel, scales = "free_x", space = "free_x", labeller = label_parsed) +
  scale_y_discrete(limits = rev(cat_order)) +
  scale_size(range = c(4.5, 12), limits = c(1, max(bub_plot$n)), guide = "none") +
  scale_fill_manual(values = c(up = "#B2182B", down = "#2166AC"), breaks = c("up", "down"),
                    labels = c(up = "Higher in taxon-positive calves", down = "Lower in taxon-positive calves"),
                    name = NULL) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  labs(x = NULL, y = NULL,
       title = paste0("Host genes DE with taxon detection (", str_replace_all(fig6c_bubble_comparison, "_", " "),
                      "; FDR < ", fdr_cutoff, ", |log2FC| > ", fig6c_bubble_lfc, ")")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        panel.spacing.x = unit(2, "mm"),
        strip.text = element_text(size = 8, face = "bold"),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.text.x = element_text(size = 8, colour = "grey20", lineheight = 0.9),
        plot.title = element_text(size = 11, face = "bold"),
        legend.position = "bottom", legend.text = element_text(size = 9))

fig6c_bub_w <- 3.8 + 1.4 * nrow(taxon_tot)     # room for the longest genus names
fig6c_bub_h <- 2.2 + 0.32 * length(cat_order)
ggsave(file.path(out_dir, "fig6C_category_bubble.png"), p_6c_bub,
       width = fig6c_bub_w, height = fig6c_bub_h, dpi = 300, bg = "white")
cat("\nSaved: fig6C_category_bubble.png, fig6C_category_counts.csv\n")


############################################################
# 6 — CANDIDATE GENES PER TAXON AND CONTRAST (table)
############################################################

DE_6c %>%
  annotate() %>%
  group_by(taxon, comparison) %>%
  group_modify(~ bind_rows(
    .x %>% filter(logFC > 0) %>% slice_max(logFC, n = top_table_n, with_ties = FALSE),
    .x %>% filter(logFC < 0) %>% slice_min(logFC, n = top_table_n, with_ties = FALSE))) %>%
  ungroup() %>%
  arrange(taxon, comparison, desc(logFC)) %>%
  dplyr::select(taxon, comparison, gene = gene_clean, gene_label, category, logFC, FDR) %>%
  write_csv(file.path(out_dir, "top_up_down_genes_per_pathogen.csv"))
cat("\nSaved: top_up_down_genes_per_pathogen.csv (", top_table_n, "up / down per taxon and contrast in",
    paste(fig6c_comparisons, collapse = " / "), ")\n")


############################################################
# 7 — SUPPLEMENTARY TABLE: ALL DE RESULTS
############################################################

section("SUPPLEMENTARY TABLE: ALL DE RESULTS")

supp_de <- DE_all %>%
  mutate(gene_clean = str_remove(gene, "^gene-"),
         # flags from the unrounded values
         DE      = FDR < fdr_cutoff & abs(logFC) > lfc_cutoff,
         high_FC = FDR < fdr_cutoff & abs(logFC) > lfc_cutoff_high) %>%
  left_join(gene_annot %>% dplyr::select(gene_clean, gene_label, category), by = "gene_clean") %>%
  arrange(factor(comparison, levels = comparison_order), taxon, FDR) %>%
  transmute(comparison, taxon,
            reference_group = group_1_name, n_reference = n_group_1,
            test_group = group_2_name, n_test = n_group_2,
            gene = gene_clean, gene_label = coalesce(gene_label, gene_clean),
            category,                                   # NA = gene not annotated (not DE, not in a network)
            log2FC = round(logFC, 4), average_log2CPM = round(AveExpr, 3), t = round(t, 3),
            P_value = signif(PValue, 4), FDR = signif(FDR, 4), B = round(B, 3),
            DE, high_FC)

write_csv(supp_de, file.path(out_dir, "TableS_DE_results_all.csv.gz"))
write_csv(supp_de %>% filter(DE), file.path(out_dir, "TableS_DE_results_significant.csv"))
cat("log2FC = test_group vs reference_group; FDR = Benjamini-Hochberg within each test (taxon x contrast)\n",
    "DE: FDR <", fdr_cutoff, "and |log2FC| >", lfc_cutoff, "| high_FC: FDR <", fdr_cutoff, "and |log2FC| >", lfc_cutoff_high, "\n")
cat("All results:", nrow(supp_de), "rows (", n_distinct(paste(supp_de$comparison, supp_de$taxon)), "tests ) |",
    "significant:", sum(supp_de$DE), "rows\n")
cat("Saved: TableS_DE_results_all.csv.gz, TableS_DE_results_significant.csv\n")


############################################################
# 8 — FIGURE 7: HOST-TAXON CORRELATION NETWORKS (04a)
# Edges from 04a (all genes); gene labels and categories from 05a. Nothing is re-computed.
############################################################

section("FIGURE 7: CORRELATION NETWORKS (04a, all genes)")

# viral node names (abbreviations) from taxa_long.csv
virus_abbrev <- read_csv(file.path(in_dir, "taxa_long.csv"), col_select = c(taxon_label, type),
                         col_types = cols(.default = col_character())) %>%
  filter(type == "virus") %>% pull(taxon_label) %>% unique()

shape_values <- c(
  "Interferon/Antiviral"          = 16, "Innate Sensing/Inflammasome"   = 17,
  "Cytokine/Inflammation"         = 15, "Innate Humoral Defense"        = 18,
  "Adaptive/Leukocyte"            = 8,  "Autophagy/Lysosomal Defense"   = 7,
  "Ciliary/Mucociliary"           = 9,  "Epithelial Barrier"            = 10,
  "ECM/Matrix/Proteolysis"        = 11, "Ion/Solute Transport"          = 12,
  "Metabolism"                    = 13, "Signaling/Trafficking/Cytoskeleton/Adhesion" = 14,
  "Apoptosis/Stress/Remodeling"   = 3,  "Neuronal/Sensory/Reproductive" = 4,
  "Transcription/RNA/Translation" = 5,  "Other/Uncharacterised"         = 1,
  "Innate Immunity (general)"     = 6,  "Immune Response (general)"     = 0, "Coagulation/Haemostasis" = 2,
  "Taxon"                         = 19
)
missing_shapes <- setdiff(unique(gene_annot$category), names(shape_values))
if (length(missing_shapes) > 0) fail("No node shape for categories: ", paste(missing_shapes, collapse = ", "))

# genes DE for the same taxon in positive_vs_negative: FDR < fdr_cutoff, no fold-change threshold
# (the same definition as the edge-level overlap in 04a); logFC kept for the direction check
fig7_DE_same_taxon <- DE_taxon %>% filter(comparison == "positive_vs_negative") %>%
  distinct(taxon, gene, .keep_all = TRUE) %>%
  transmute(taxon, gene, logFC_DE = logFC, DE_same_taxon = TRUE)

# thresholds used by 04a (for the figure title)
net_params_file <- file.path(net_dir, "network_parameters.csv")
net_params <- if (file.exists(net_params_file))
  deframe(read_csv(net_params_file, col_types = cols(.default = col_character()))) else character()
np <- function(x) if (x %in% names(net_params)) net_params[[x]] else "?"

# layout weights: full-colour taxon-gene edges of taxa with >= fig7_bold_min_genes such edges pull
# with fig7_bold_edge_weight (edge attribute push), all others with 1
fig7_edge_weights <- function(sg) {
  ea <- igraph::edge_attr(sg)
  if (igraph::ecount(sg) == 0) return(NULL)
  ifelse(ea$push, fig7_bold_edge_weight, 1)
}

# full-colour genes of crowded taxa (edge attribute push), linked to that taxon only: spread evenly on a
# ring around it (radius fig7_bold_ring x its pastel genes' 90th-centile distance), over the largest arc
# free of the other taxa and of genes shared with other taxa (15 degrees clear of them), at most
# fig7_bold_arc degrees; all around if nothing else is attached; ordered by rho
fig7_spread_bold <- function(sg, xy) {
  if (fig7_bold_arc <= 0 || igraph::ecount(sg) == 0) return(xy)
  vt  <- igraph::V(sg)$type
  el  <- igraph::as_edgelist(sg, names = FALSE)
  ea  <- igraph::edge_attr(sg)
  deg <- igraph::degree(sg)
  taxa <- which(vt != "host_gene")
  for (t in taxa) {
    at_t  <- el[, 1] == t | el[, 2] == t
    other <- ifelse(el[, 1] == t, el[, 2], el[, 1])
    sel   <- at_t & ea$push & deg[other] == 1
    ring  <- other[sel]
    if (length(ring) < 2) next
    c0   <- xy[t, ]
    dist <- function(i) sqrt((xy[i, 1] - c0[1])^2 + (xy[i, 2] - c0[2])^2)
    single <- other[at_t & !ea$push & deg[other] == 1]
    R <- fig7_bold_ring * if (length(single) >= 3) unname(quantile(dist(single), 0.9)) else median(dist(ring)) / fig7_bold_ring
    obst <- unique(c(setdiff(taxa, t), other[at_t & deg[other] > 1]))     # other taxa, shared genes
    n <- length(ring); margin <- 15 * pi / 180; amax <- min(fig7_bold_arc, 360) * pi / 180
    if (length(obst) == 0) {
      ang <- 2 * pi * (seq_len(n) - 1) / n + pi / 2
    } else {
      a   <- sort(atan2(xy[obst, 2] - c0[2], xy[obst, 1] - c0[1]) %% (2 * pi))
      gap <- diff(c(a, a[1] + 2 * pi))
      k   <- which.max(gap)
      w   <- min(amax, max(gap[k] - 2 * margin, gap[k] / 2))
      mid <- a[k] + gap[k] / 2
      ang <- mid - w / 2 + w * (seq_len(n) - 0.5) / n
    }
    ord <- order(-ea$corr[sel])
    xy[ring[ord], 1] <- c0[1] + R * cos(ang)
    xy[ring[ord], 2] <- c0[2] + R * sin(ang)
  }
  xy
}

# "packed_right" / "packed_bottom" layout: each connected part laid out with FR, centred and
# scaled to radius sqrt(number of nodes); the largest part placed first, the others (largest
# first) in a column to its right / a row below it, wrapping when the column / row is full
fig7_packed_layout <- function(g, side = "right") {
  comp  <- igraph::components(g)
  parts <- order(comp$csize, decreasing = TRUE)
  lays <- lapply(parts, function(k) {
    v  <- which(comp$membership == k)
    sg <- igraph::induced_subgraph(g, v)
    xy <- if (length(v) == 1) matrix(0, 1, 2) else fig7_spread_bold(sg, igraph::layout_with_fr(sg, weights = fig7_edge_weights(sg)))
    xy <- sweep(xy, 2, (apply(xy, 2, max) + apply(xy, 2, min)) / 2)
    r  <- max(sqrt(rowSums(xy^2)))
    if (r > 0) xy <- xy / r * sqrt(length(v))
    list(v = v, xy = xy)
  })
  bb <- function(xy) c(xmin = min(xy[, 1]), xmax = max(xy[, 1]), ymin = min(xy[, 2]), ymax = max(xy[, 2]))
  if (fig7_mirror) lays[[1]]$xy[, 1] <- -lays[[1]]$xy[, 1]

  # place the other parts beside ("right") / below ("bottom") the largest part
  pack <- function(lays, side) {
    b0  <- bb(lays[[1]]$xy)
    gap <- fig7_pack_gap * max(b0["xmax"] - b0["xmin"], b0["ymax"] - b0["ymin"])
    out <- matrix(NA_real_, igraph::vcount(g), 2)
    out[lays[[1]]$v, ] <- lays[[1]]$xy
    if (length(lays) > 1) {
      if (side == "right") {
        col_x <- b0["xmax"] + 2 * gap; cursor <- b0["ymax"]; col_w <- 0
        for (p in lays[-1]) {
          b <- bb(p$xy); w <- b["xmax"] - b["xmin"]; h <- b["ymax"] - b["ymin"]
          if (cursor - h < b0["ymin"] && cursor < b0["ymax"]) { col_x <- col_x + col_w + gap; cursor <- b0["ymax"]; col_w <- 0 }
          out[p$v, 1] <- p$xy[, 1] - b["xmin"] + col_x
          out[p$v, 2] <- p$xy[, 2] - b["ymax"] + cursor
          cursor <- cursor - h - gap; col_w <- max(col_w, w)
        }
      } else {
        row_y <- b0["ymin"] - 2 * gap; cursor <- b0["xmin"]; row_h <- 0
        for (p in lays[-1]) {
          b <- bb(p$xy); w <- b["xmax"] - b["xmin"]; h <- b["ymax"] - b["ymin"]
          if (cursor + w > b0["xmax"] && cursor > b0["xmin"]) { row_y <- row_y - row_h - gap; cursor <- b0["xmin"]; row_h <- 0 }
          out[p$v, 1] <- p$xy[, 1] - b["xmin"] + cursor
          out[p$v, 2] <- p$xy[, 2] - b["ymax"] + row_y
          cursor <- cursor + w + gap; row_h <- max(row_h, h)
        }
      }
    }
    stopifnot(!anyNA(out))
    out
  }
  # candidate arrangements: the largest part as laid out and, with fig7_rotate_largest, turned (long axis
  # horizontal / vertical, or 90 degrees); the other parts to its right and/or below ("packed_auto": both).
  # Kept: the one that fills a panel of height / width fig7_target_aspect best (largest drawing scale,
  # i.e. smallest max(width, height / aspect)); the first (as laid out, parts to the right) unless
  # another is > 10% better
  orient <- list("as laid out" = lays[[1]]$xy)
  if (fig7_rotate_largest && nrow(lays[[1]]$xy) > 2) {
    pcs <- prcomp(lays[[1]]$xy, center = TRUE, scale. = FALSE)$x[, 1:2, drop = FALSE]
    orient <- c(orient, list("long axis horizontal" = pcs, "long axis vertical" = cbind(-pcs[, 2], pcs[, 1]),
                             "turned 90 degrees" = cbind(-lays[[1]]$xy[, 2], lays[[1]]$xy[, 1])))
  }
  sides <- if (side == "auto") c("right", "bottom") else side
  cand  <- expand.grid(o = names(orient), s = sides, stringsAsFactors = FALSE)
  cand  <- cand[order(cand$o != "as laid out", cand$s != sides[1]), ]
  score <- function(o) { b <- bb(o); max(b["xmax"] - b["xmin"], (b["ymax"] - b["ymin"]) / fig7_target_aspect) }
  outs  <- lapply(seq_len(nrow(cand)), function(i) { l <- lays; l[[1]]$xy <- orient[[cand$o[i]]]; pack(l, cand$s[i]) })
  sc    <- vapply(outs, score, numeric(1))
  best  <- if (min(sc) < 0.9 * sc[1]) which.min(sc) else 1
  if (nrow(cand) > 1)
    cat("Layout: largest part ", cand$o[best], ", other parts ", if (cand$s[best] == "right") "to its right" else "below",
        if (best > 1) paste0(" (drawing scale x", round(sc[1] / sc[best], 2), " vs as laid out)") else "", "\n", sep = "")
  outs[[best]]
}

fig7_genes <- list()
fig7_data  <- list()     # per status: layout, edges, nodes (plotted together after the loop)
for (st in fig7_statuses) {
  f <- file.path(net_dir, paste0("network_", st, "_stable_bootstrap_edges.csv"))
  if (!file.exists(f)) { cat("\n", st, ": edge file not found (", f, ") — skipped\n", sep = ""); next }
  e <- read_csv(f, show_col_types = FALSE)
  if (!"stable_strict" %in% names(e)) {
    if (fig7_highlight_rule == "strict")
      fail("Column stable_strict missing in ", f, " — re-run the current 04a_networks_DE_genes.R.")
    e <- e %>% mutate(stability_strict = NA_real_, stable_strict = NA)
  }
  if (nrow(e) == 0) { cat("\n", st, ": no stable edges — skipped\n", sep = ""); next }

  # taxon first; taxon-taxon edges (bacteria-bacteria, virus-virus, bacteria-virus) kept as they are
  e <- e %>%
    mutate(is_gene_from = str_detect(from, "^gene-"), is_gene_to = str_detect(to, "^gene-"),
           taxon   = if_else(is_gene_from, to, from),
           partner = if_else(is_gene_from, from, to),
           edge_kind = if_else(!is_gene_from & !is_gene_to, "taxon-taxon", "taxon-gene"))
  n_all <- nrow(e)
  n_tt  <- sum(e$edge_kind == "taxon-taxon")
  # display filter on |rho| (all edge kinds)
  e <- e %>% filter(abs(corr) >= fig7_min_abs_rho)
  cat("\n---", st, "--- stable edges:", n_all, "(taxon-taxon:", n_tt, ") | with |rho| >=", fig7_min_abs_rho, ":",
      nrow(e), "(taxon-taxon:", sum(e$edge_kind == "taxon-taxon"), ")\n")
  if (nrow(e) == 0) { cat(st, ": no edges with |rho| >=", fig7_min_abs_rho, "— skipped\n"); next }
  cat("Taxon-gene edges per taxon at |rho| >=", fig7_min_abs_rho, ":\n")
  e %>% filter(edge_kind == "taxon-gene") %>% count(taxon, sort = TRUE) %>% as.data.frame() %>% print(row.names = FALSE)
  n_rho <- nrow(e)
  if (is.finite(fig7_max_genes_per_taxon)) {
    e <- bind_rows(
      e %>% filter(edge_kind == "taxon-taxon"),                       # always all drawn
      e %>% filter(edge_kind == "taxon-gene") %>% group_by(taxon) %>%
        slice_max(abs(corr), n = fig7_max_genes_per_taxon, with_ties = FALSE) %>% ungroup())
  }

  # ---- full colour (bold) vs pastel (fig7_highlight_rule) ----
  e <- e %>%
    mutate(edge_id = row_number(),
           gene_clean = if_else(edge_kind == "taxon-gene", str_remove(partner, "^gene-"), NA_character_)) %>%
    left_join(gene_annot %>% dplyr::select(gene_clean, gene_label, category), by = "gene_clean") %>%
    left_join(fig7_DE_same_taxon, by = c("taxon", "partner" = "gene")) %>%
    mutate(DE_same_taxon = coalesce(DE_same_taxon, FALSE),
           # DE for the same taxon AND sign of log2FC (positive vs negative) = sign of rho
           DE_same_direction = DE_same_taxon & coalesce(sign(logFC_DE) == sign(corr), FALSE),
           DE_ok = edge_kind == "taxon-taxon" | !fig7_highlight_DE_only | DE_same_taxon)
  if (fig7_highlight_rule == "strict") {
    # edges also stable at the stricter 04a level (COR_THRESH_STRICT)
    e <- e %>% mutate(highlight = coalesce(stable_strict, FALSE) & DE_ok)
  } else if (fig7_highlight_rule == "rho") {
    # every edge (taxon-gene and taxon-taxon) with |rho| >= fig7_highlight_min_abs_rho
    e <- e %>% mutate(highlight = abs(corr) >= fig7_highlight_min_abs_rho & DE_ok)
  } else {
    # all taxon-taxon edges + per taxon, the strongest annotated genes (all of them if
    # fig7_highlight_per_taxon = Inf); categories in fig7_highlight_exclude never highlighted
    hl_ids <- e %>%
      filter(edge_kind == "taxon-gene", !is.na(category), !category %in% fig7_highlight_exclude, DE_ok) %>%
      group_by(taxon) %>%
      slice_max(abs(corr), n = min(fig7_highlight_per_taxon, nrow(e)), with_ties = FALSE) %>% ungroup() %>% pull(edge_id)
    e <- e %>% mutate(highlight = edge_kind == "taxon-taxon" | edge_id %in% hl_ids)
  }
  # gene labels: genes on full-colour edges, at most fig7_label_per_taxon per taxon
  lab_e <- e %>% filter(highlight, edge_kind == "taxon-gene")
  if (is.finite(fig7_label_per_taxon))
    lab_e <- lab_e %>% group_by(taxon) %>% slice_max(abs(corr), n = fig7_label_per_taxon, with_ties = FALSE) %>% ungroup()
  lab_genes <- unique(lab_e$partner)
  cat("Drawn:", nrow(e),
      if (nrow(e) < n_rho) paste0("(all taxon-taxon + strongest ", fig7_max_genes_per_taxon, " gene edges per taxon)") else "",
      "| full colour:", sum(e$highlight), "| genes labelled:", length(lab_genes), "\n")

  nodes <- tibble(name = unique(c(e$from, e$to))) %>%
    mutate(type = case_when(str_detect(name, "^gene-") ~ "host_gene",
                            name %in% virus_abbrev       ~ "virus",
                            TRUE                         ~ "bacteria"),
           gene_clean = str_remove(name, "^gene-")) %>%
    left_join(gene_annot %>% dplyr::select(gene_clean, gene_label, category), by = "gene_clean") %>%
    mutate(category = case_when(type != "host_gene" ~ "Taxon",
                                TRUE ~ coalesce(category, "Other/Uncharacterised")),
           label = if_else(type == "host_gene", coalesce(gene_label, gene_clean), name),
           highlight = type != "host_gene" | name %in% c(e$from[e$highlight], e$to[e$highlight]),
           labelled  = type != "host_gene" | name %in% lab_genes,
           # labels: taxa larger and bold (bacterial genera in italics)
           label_face = case_when(type == "bacteria" ~ "bold.italic", type == "virus" ~ "bold", TRUE ~ "plain"))
  # stable taxon-gene edges drawn per taxon (node size)
  tax_n <- e %>% filter(edge_kind == "taxon-gene") %>% group_by(taxon) %>%
    summarise(n_edges = n(), n_strict = sum(coalesce(stable_strict, FALSE)), .groups = "drop")
  nodes <- nodes %>% left_join(tax_n, by = c("name" = "taxon")) %>%
    mutate(n_edges  = if_else(type == "host_gene", NA_integer_, coalesce(n_edges, 0L)),
           n_strict = if_else(type == "host_gene", NA_integer_, coalesce(n_strict, 0L)),
           size_val = if (fig7_taxon_size_by == "strict") n_strict else n_edges)
  # node types must agree with those written by 04a (guards against a gene ID without
  # the "gene-" prefix being drawn as a taxon)
  nf <- file.path(net_dir, paste0("network_", st, "_stable_bootstrap_nodes.csv"))
  if (file.exists(nf)) {
    n04 <- read_csv(nf, show_col_types = FALSE, col_types = cols(.default = col_character()))
    type_check <- nodes %>% dplyr::select(name, type) %>%
      left_join(n04 %>% dplyr::select(name, type_04a = type), by = "name") %>%
      filter(is.na(type_04a) | type != type_04a)
    if (nrow(type_check) > 0) fail("Node types differ from 04a (", st, "): ",
                                   paste0(type_check$name, " ", type_check$type, "/", type_check$type_04a, collapse = ", "))
  }
  n_unannot <- sum(nodes$type == "host_gene" & !nodes$gene_clean %in% gene_annot$gene_clean)
  if (n_unannot > 0) fail(n_unannot, " network genes (", st, ") not in the 05a annotation (re-run 05a): ",
                          paste(setdiff(nodes$gene_clean[nodes$type == "host_gene"], gene_annot$gene_clean), collapse = ", "))
  cat("Nodes:", nrow(nodes), "| taxa:", sum(nodes$type != "host_gene"), "| genes:", sum(nodes$type == "host_gene"),
      "| genes in full colour:", sum(nodes$type == "host_gene" & nodes$highlight), "\n")
  cat("Stable taxon-gene edges per taxon (node size):",
      paste0(nodes$name[nodes$type != "host_gene"], " ", nodes$n_edges[nodes$type != "host_gene"],
             " (", nodes$n_strict[nodes$type != "host_gene"], " also strict)", collapse = "; "), "\n")
  cat("Genes per category (all drawn | in full colour):\n")
  nodes %>% filter(type == "host_gene") %>%
    group_by(category) %>% summarise(drawn = n(), full_colour = sum(highlight), .groups = "drop") %>%
    arrange(desc(drawn)) %>% as.data.frame() %>% print(row.names = FALSE)
  cat("\nEdges in full colour (taxon | partner | category | rho | stability | DE same taxon | DE same direction;",
      "DE: FDR <", fdr_cutoff, "in positive_vs_negative, no fold-change threshold):\n")
  e %>% filter(highlight) %>% arrange(taxon, desc(abs(corr))) %>%
    transmute(taxon, partner = coalesce(gene_label, partner), category, rho = round(corr, 2), stability,
              DE_same_taxon = if_else(edge_kind == "taxon-gene", as.character(DE_same_taxon), "-"),
              DE_same_direction = if_else(edge_kind == "taxon-gene", as.character(DE_same_direction), "-")) %>%
    as.data.frame() %>% print(row.names = FALSE)
  no_hl <- setdiff(unique(e$taxon[e$edge_kind == "taxon-gene"]), unique(e$taxon[e$highlight & e$edge_kind == "taxon-gene"]))
  if (length(no_hl) > 0) cat("Taxa with gene edges but none in full colour:", paste(no_hl, collapse = ", "), "\n")

  # layout only: full-colour gene edges of taxa with several of them are drawn longer (fig7_bold_edge_weight)
  e <- e %>% group_by(taxon) %>%
    mutate(push = highlight & edge_kind == "taxon-gene" &
             sum(highlight & edge_kind == "taxon-gene") >= fig7_bold_min_genes) %>% ungroup()
  if (fig7_bold_edge_weight != 1 && any(e$push))
    cat("Full-colour gene edges drawn longer for:", paste(unique(e$taxon[e$push]), collapse = ", "), "\n")
  g <- igraph::graph_from_data_frame(e %>% dplyr::select(from, to, corr, highlight, edge_kind, push),
                                     vertices = nodes, directed = FALSE)
  set.seed(fig7_seed)
  layout_net <- if (fig7_layout %in% c("packed_auto", "packed_right", "packed_bottom")) {
    xy <- fig7_packed_layout(g, side = sub("packed_", "", fig7_layout))
    ggraph::create_layout(g, layout = "manual", x = xy[, 1], y = xy[, 2])
  } else if (fig7_layout == "fr_packed") {
    # each connected part laid out separately, then packed at random (igraph::merge_coords)
    comp <- igraph::components(g)
    subs <- lapply(seq_len(comp$no), function(k) igraph::induced_subgraph(g, which(comp$membership == k)))
    lays <- lapply(subs, function(sg) if (igraph::vcount(sg) == 1) matrix(0, 1, 2) else
      fig7_spread_bold(sg, igraph::layout_with_fr(sg, weights = fig7_edge_weights(sg))))
    xy   <- igraph::merge_coords(subs, lays)
    xy   <- xy[match(igraph::V(g)$name, unlist(lapply(subs, function(sg) igraph::V(sg)$name))), , drop = FALSE]
    ggraph::create_layout(g, layout = "manual", x = xy[, 1], y = xy[, 2])
  } else {
    tryCatch(ggraph::create_layout(g, layout = fig7_layout),
             error = function(err) { message(fig7_layout, " layout failed — using KK"); ggraph::create_layout(g, layout = "kk") })
  }
  # check: node positions belong to the right nodes (edge endpoints match node names)
  stopifnot(identical(as.character(layout_net$name), igraph::V(g)$name))
  fig7_data[[st]] <- list(layout = layout_net, e = e, nodes = nodes, n_all = n_all)

  e_out <- e %>%
    transmute(status = st, taxon, partner = str_remove(partner, "^gene-"), edge_kind, gene_label, category,
              rho = corr, p_adj, stability, stability_strict, stable_strict,
              DE_same_taxon = if_else(edge_kind == "taxon-gene", DE_same_taxon, NA),
              DE_same_direction = if_else(edge_kind == "taxon-gene", DE_same_direction, NA),
              full_colour = highlight,
              labelled = edge_kind == "taxon-gene" & paste0("gene-", partner) %in% lab_genes)
  write_csv(e_out, file.path(out_dir, paste0("fig7_edges_", st, ".csv")))
  cat("Saved: fig7_edges_", st, ".csv\n", sep = "")
  # genes named in the figure (labelled) go to the annotation check
  fig7_genes[[st]] <- e_out %>% filter(edge_kind == "taxon-gene", labelled) %>%
    transmute(gene = partner, figure = "7", context = paste0(st, ": ", taxon))
}
fig7_genes <- bind_rows(fig7_genes)

# ---- Figure 7: one file per status, or panels combined; same rho / |rho| scales throughout ----
if (length(fig7_data) > 0) {
  hl_all <- bind_rows(lapply(fig7_data, function(d) d$e %>% filter(highlight)))
  rho_lim   <- if (nrow(hl_all) > 0) max(abs(hl_all$corr)) else 1
  width_lim <- if (nrow(hl_all) > 0) range(abs(hl_all$corr)) else c(0, 1)
  # legend entries: categories of full-colour genes and node types, over all panels
  cats_7  <- intersect(setdiff(names(shape_values), "Taxon"),
                       unique(unlist(lapply(fig7_data, function(d) d$nodes$category[d$nodes$type == "host_gene" & d$nodes$highlight]))))
  types_7 <- intersect(c("bacteria", "virus", "host_gene"),
                       unique(unlist(lapply(fig7_data, function(d) d$nodes$type[d$nodes$type != "host_gene" | d$nodes$highlight]))))
  faded_any <- any(vapply(fig7_data, function(d) !all(d$e$highlight), logical(1)))
  # taxon node size: same scale in every figure (area ~ number of edges)
  size_max <- max(1, unlist(lapply(fig7_data, function(d) d$nodes$size_val)), na.rm = TRUE)
  size_breaks <- pretty(c(0, size_max), n = 4)
  size_breaks <- size_breaks[size_breaks > 0 & size_breaks <= size_max]
  if (length(size_breaks) == 0) size_breaks <- size_max
  size_title <- if (fig7_taxon_size_by == "strict") paste0("Edges stable at\n|rho| >= ", np("COR_THRESH_STRICT")) else "Stable edges"
  # drawn size of each taxon (same formula as scale_size with limits c(0, size_max)), used to keep
  # taxon names clear of their circles
  taxon_pt <- function(v) if (fig7_taxon_size_by == "none") rep(7, length(v)) else
    fig7_taxon_size_range[1] + diff(fig7_taxon_size_range) * sqrt(pmax(v, 0) / size_max)
  fig7_one_label_layer <- utils::packageVersion("ggplot2") >= "3.5.0"

  # categories / node types of the full-colour nodes of one panel (for "own" legends)
  cats_of  <- function(d) intersect(setdiff(names(shape_values), "Taxon"),
                                    unique(d$nodes$category[d$nodes$type == "host_gene" & d$nodes$highlight]))
  types_of <- function(d) intersect(c("bacteria", "virus", "host_gene"),
                                    unique(d$nodes$type[d$nodes$type != "host_gene" | d$nodes$highlight]))

  plot_network <- function(d, title, cats = cats_7, types = types_7) {
    ggraph::ggraph(d$layout) +
      # faded: pastel red (positive) / blue (negative)
      ggraph::geom_edge_link(aes(filter = !highlight & corr > 0), colour = fig7_faded_cols[["positive"]], width = 0.3, alpha = 0.6) +
      ggraph::geom_edge_link(aes(filter = !highlight & corr < 0), colour = fig7_faded_cols[["negative"]], width = 0.3, alpha = 0.6) +
      ggraph::geom_edge_link(aes(filter = highlight, width = abs(corr), colour = corr), alpha = 0.9) +
      ggraph::scale_edge_colour_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                                          limits = c(-rho_lim, rho_lim), name = "rho") +
      ggraph::scale_edge_width(range = c(0.6, 2.2), limits = width_lim, name = "|rho|") +
      # taxa first (underneath), so genes next to a large hub stay visible;
      # size (area) ~ number of stable edges (fig7_taxon_size_by), or all the same size
      (if (fig7_taxon_size_by == "none")
         ggraph::geom_node_point(aes(filter = type != "host_gene", colour = type), shape = 19, size = 7)
       else ggraph::geom_node_point(aes(filter = type != "host_gene", colour = type, size = size_val), shape = 19)) +
      ggraph::geom_node_point(aes(filter = type == "host_gene" & !highlight), colour = fig7_faded_cols[["gene"]], size = 1) +
      # show.legend = TRUE: legend keys for every category, including those absent from this panel
      ggraph::geom_node_point(aes(filter = type == "host_gene" & highlight, colour = type, shape = category),
                              size = 3.5, show.legend = c(shape = TRUE, colour = TRUE)) +
      (if (fig7_one_label_layer) {
        # gene and taxon labels in ONE repel layer, so neither overlaps the other; text size given
        # with I() so it bypasses the node size scale (needs ggplot2 >= 3.5.0)
        ggraph::geom_node_text(aes(filter = labelled, label = label, fontface = label_face,
                                   size = I(if_else(type == "host_gene", 2.7, 4.2)),
                                   point.size = if_else(type == "host_gene", 3.5, taxon_pt(size_val))),
                               repel = TRUE, max.overlaps = Inf, segment.colour = "grey60", seed = fig7_seed,
                               min.segment.length = 0.2, box.padding = 0.3, point.padding = 0.05)
      } else list(
        # ggplot2 < 3.5: two layers. Gene labels; taxa are included with an empty label so the labels
        # are pushed off the hubs
        ggraph::geom_node_text(aes(filter = (labelled & type == "host_gene") | type != "host_gene",
                                   label = if_else(type == "host_gene", label, ""),
                                   point.size = if_else(type == "host_gene", 3.5, taxon_pt(size_val))), size = 2.7,
                               repel = TRUE, max.overlaps = Inf, segment.colour = "grey60", seed = fig7_seed,
                               min.segment.length = 0.2, box.padding = 0.3, point.padding = 0.05),
        # taxon labels, kept clear of their (sized) circles
        ggraph::geom_node_text(aes(filter = type != "host_gene", label = label, fontface = label_face,
                                   point.size = taxon_pt(size_val)), size = 4.2,
                               repel = TRUE, max.overlaps = Inf, segment.colour = "grey60", seed = fig7_seed,
                               min.segment.length = 0.2, box.padding = 0.25, point.padding = 0.1))) +
      (if (fig7_taxon_size_by == "none") NULL
       else scale_size(range = fig7_taxon_size_range, limits = c(0, size_max), breaks = size_breaks,
                       name = size_title)) +
      scale_colour_manual(values = c(virus = "#E1A730", bacteria = "#2A9D8F", host_gene = "#6A4C93"),
                          limits = types,
                          labels = c(virus = "virus", bacteria = "bacterium", host_gene = "host gene"), name = NULL) +
      scale_shape_manual(values = shape_values, limits = cats, name = "Category") +
      guides(colour = guide_legend(order = 1, override.aes = list(size = 5, shape = 16)),
             size   = guide_legend(order = 2, override.aes = list(colour = "grey45", shape = 19)),
             shape  = guide_legend(order = 3, override.aes = list(size = 4, colour = "#6A4C93")),
             edge_colour = guide_edge_colourbar(order = 4), edge_width = guide_legend(order = 5)) +
      scale_x_continuous(expand = expansion(mult = fig7_expand)) +
      scale_y_continuous(expand = expansion(mult = fig7_expand)) +
      (if (fig7_fixed_aspect) coord_fixed(clip = "off") else NULL) +
      labs(title = title) +
      theme_void() +
      theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            legend.text = element_text(size = 9), legend.key.height = grid::unit(0.45, "cm"),
            plot.margin = margin(5, 5, 5, 5))
  }
  fig7_plots <- lapply(names(fig7_data), function(st) plot_network(fig7_data[[st]], str_to_sentence(st)))

  # one legend, taken from the first panel (all panels share the same scale limits)
  g1 <- ggplotGrob(fig7_plots[[1]])
  leg_idx <- which(grepl("^guide-box", vapply(g1$grobs, function(x) x$name, character(1))))
  leg_idx <- leg_idx[!vapply(g1$grobs[leg_idx], inherits, logical(1), "zeroGrob")]
  fig7_legend <- g1$grobs[[leg_idx[1]]]
  fig7_panels <- lapply(fig7_plots, function(p) ggplotGrob(p + theme(legend.position = "none")))

  # figure title with the display parameters (fig7_show_title <- FALSE to drop it)
  fig7_title <- paste0(
    "Host-taxon correlation networks: stable edges (|rho| >= ", np("COR_THRESH"), ", BH-adjusted P <= ", np("P_THRESH"),
    ", same sign and |rho| >= ", np("COR_THRESH"), " in >= ", as.numeric(np("STABILITY_THRESHOLD")) * 100, "% of ",
    np("N_BOOTSTRAP_ITER"), " bootstrap resamples)",
    if (isTRUE(fig7_min_abs_rho > suppressWarnings(as.numeric(np("COR_THRESH"))))) paste0(" with |rho| >= ", fig7_min_abs_rho) else "",
    if (faded_any) " (pastel red / blue = positive / negative rho)" else "", ";\n",
    if (fig7_highlight_rule == "strict")
      paste0("full colour: also stable at |rho| >= ", np("COR_THRESH_STRICT"), " (same sign and |rho| >= ",
             np("COR_THRESH_STRICT"), " in >= ", as.numeric(np("STABILITY_THRESHOLD")) * 100, "% of the same resamples)")
    else if (fig7_highlight_rule == "rho") paste0("full colour: |rho| >= ", fig7_highlight_min_abs_rho)
    else if (is.finite(fig7_highlight_per_taxon))
      paste0("full colour: taxon-taxon edges and the ", fig7_highlight_per_taxon, " annotated genes per taxon with the strongest |rho|")
    else "full colour: taxon-taxon edges and all annotated genes",
    if (fig7_highlight_DE_only) paste0(" (genes DE for the same taxon, FDR < ", fdr_cutoff, ")") else "",
    if (is.finite(fig7_label_per_taxon)) paste0("; labels: up to ", fig7_label_per_taxon, " genes per taxon") else "",
    "; node positions and edge lengths are not informative",
    if (fig7_bold_edge_weight != 1 && fig7_layout %in% c("packed_auto", "packed_right", "packed_bottom", "fr_packed"))
      paste0(" (full-colour genes of taxa with >= ", fig7_bold_min_genes, " of them drawn further out",
             if (fig7_bold_arc > 0) ", on a ring around the taxon ordered by rho," else "",
             " so that the gene labels do not overlap)") else "")
  legend_w_in <- 3.2
  size_in <- if (!is.null(fig7_size_in)) fig7_size_in else
    switch(fig7_arrangement, side = c(width = 19, height = 9), stack = c(width = 12, height = 16),
           separate = c(width = 11, height = 10))
  # title wrapped to the figure width (about 13 characters per inch at 11 pt)
  title_wrapped <- paste(unlist(lapply(strsplit(fig7_title, "\n")[[1]], strwrap, width = floor(size_in[["width"]] * 13))),
                         collapse = "\n")
  title_grob <- if (fig7_show_title) textGrob(title_wrapped, gp = gpar(fontsize = 11)) else NULL
  n_nodes <- vapply(fig7_data, function(d) nrow(d$nodes), numeric(1))
  if (fig7_arrangement == "separate") {
    # one file per status, each with its own legend (fig7_legend_categories)
    for (st in names(fig7_data)) {
      d <- fig7_data[[st]]
      p_st <- if (fig7_legend_categories == "own") plot_network(d, str_to_sentence(st), cats_of(d), types_of(d)) else
        plot_network(d, str_to_sentence(st))
      # extra right margin: text is measured on the default device but drawn in the PNG, whose
      # font metrics can be wider (long category names would otherwise be cut at the edge)
      gp_i <- ggplotGrob(p_st + theme(plot.margin = margin(5, 45, 5, 5)))
      g_i  <- arrangeGrob(gp_i, top = title_grob)
      size_st <- size_in
      if (is.null(fig7_size_in) && fig7_fixed_aspect) {
        # height fitted to the network: panel width (figure width minus legend and margins) x the
        # network's height / width, plus the fixed rows (header, margins, title); never less than the legend
        pp    <- ggplot_build(p_st)$layout$panel_params[[1]]
        asp   <- diff(pp$y.range) / diff(pp$x.range)
        abs_w <- sum(grid::convertWidth(gp_i$widths[grid::unitType(gp_i$widths) != "null"], "in", valueOnly = TRUE))
        abs_h <- sum(grid::convertHeight(gp_i$heights[grid::unitType(gp_i$heights) != "null"], "in", valueOnly = TRUE))
        leg   <- gp_i$grobs[grepl("^guide-box", gp_i$layout$name) & !vapply(gp_i$grobs, inherits, logical(1), "zeroGrob")]
        leg_h <- if (length(leg) > 0) grid::convertHeight(grid::grobHeight(leg[[1]]), "in", valueOnly = TRUE) else 0
        top_h <- if (is.null(title_grob)) 0 else grid::convertHeight(grid::grobHeight(title_grob), "in", valueOnly = TRUE) + 0.2
        # (legend measured on the default device; 5% + 0.2 in added for the PNG's wider font metrics)
        size_st[["height"]] <- max((size_in[["width"]] - abs_w) * asp + abs_h, leg_h * 1.05 + abs_h + 0.2) + top_h
      }
      ggsave(file.path(out_dir, paste0("network_", st, ".png")), g_i,
             width = size_st[["width"]], height = size_st[["height"]], dpi = 300, bg = "white")
      cat("Saved: network_", st, ".png (", round(size_st[["width"]], 1), " x ", round(size_st[["height"]], 1), " in)\n", sep = "")
    }
  } else {
    if (fig7_arrangement == "side") {
      # panel widths follow network size (square root of the number of nodes)
      w <- sqrt(n_nodes) / sum(sqrt(n_nodes)) * (size_in[["width"]] - legend_w_in)
      body <- arrangeGrob(grobs = fig7_panels, nrow = 1, widths = unit(w, "in"))
    } else {
      # stacked: heights follow network size
      h <- sqrt(n_nodes) / sum(sqrt(n_nodes))
      body <- arrangeGrob(grobs = fig7_panels, ncol = 1, heights = unit(h, "null"))
    }
    fig7_grob <- arrangeGrob(body, fig7_legend, nrow = 1,
                             widths = unit(c(size_in[["width"]] - legend_w_in, legend_w_in), "in"), top = title_grob)
    fig7_file <- paste0("network_", paste(names(fig7_data), collapse = "_"), ".png")
    ggsave(file.path(out_dir, fig7_file), fig7_grob,
           width = size_in[["width"]], height = size_in[["height"]], dpi = 300, bg = "white")
    cat("Saved:", fig7_file, "(", fig7_arrangement, ")\n")
  }
  cat("\nFigure 7 title (also for the figure legend):", fig7_title, "\n")
}
if (nrow(fig7_genes) == 0) fig7_genes <- tibble(gene = character(), figure = character(), context = character())


############################################################
# 9 — ANNOTATION CHECK TABLE: every gene shown in a figure
############################################################

section("ANNOTATION CHECK: GENES SHOWN IN THE FIGURES")

# which figure(s) each gene appears in, and in which context
fig_genes <- bind_rows(
  # 6A: genes in the Euler sets drawn (taxon and direction)
  enframe(sets_up,   name = "taxon", value = "gene") %>% unnest(gene) %>% mutate(figure = "6A", context = paste(taxon, "up")),
  enframe(sets_down, name = "taxon", value = "gene") %>% unnest(gene) %>% mutate(figure = "6A", context = paste(taxon, "down")),
  # 6B: respiratory vs healthy
  top_genes_6b %>% transmute(gene, figure = "6B", context = fig6b_comparison),
  # 6C: genes per taxon panel
  genes_6c %>% filter(in_expression_data) %>% transmute(gene = gene_host, figure = "6C", context = taxon),
  # 7: genes in the network figures (status: taxon)
  fig7_genes
) %>%
  mutate(gene_clean = str_remove(gene, "^gene-")) %>%
  distinct(gene_clean, figure, context)

fig_summary <- fig_genes %>%
  group_by(gene_clean) %>%
  summarise(figures = paste(sort(unique(figure)), collapse = ", "),
            figure_context = paste0(figure, " (", context, ")", collapse = "; "),
            .groups = "drop")

# annotation evidence from 05a (route, evidence, GO-supported categories, conflicts)
evidence <- read_csv(file.path(annot_dir, "gene_annotation_evidence.csv"), show_col_types = FALSE,
                     col_types = cols(.default = col_character()))

annotation_check <- fig_summary %>%
  left_join(evidence, by = c("gene_clean" = "gene")) %>%
  transmute(gene = gene_clean, gene_label, figures, figure_context,
            category, route, evidence,
            go_supported_categories, family_rule, manual_category, conflict,
            description, human_orthologues, ENTREZID,
            ncbi_link = if_else(is.na(ENTREZID), NA_character_,
                                paste0("https://www.ncbi.nlm.nih.gov/gene/", ENTREZID))) %>%
  arrange(category, gene)

write_csv(annotation_check, file.path(out_dir, "annotation_check_figure_genes.csv"), na = "")

cat("Genes shown in at least one figure:", nrow(annotation_check), "\n")
print(as.data.frame(count(fig_genes %>% distinct(gene_clean, figure), figure)), row.names = FALSE)
cat("\nBy category and annotation route:\n")
annotation_check %>% count(category, route) %>%
  pivot_wider(names_from = route, values_from = n, values_fill = 0) %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\nFlagged for review (manual or family-rule category not supported by the gene's GO terms):",
    sum(annotation_check$conflict == "TRUE", na.rm = TRUE), "\n")
cat("Other/Uncharacterised:", sum(annotation_check$category == "Other/Uncharacterised", na.rm = TRUE),
    "| missing from 05a evidence table:", sum(is.na(annotation_check$route)), "\n")
cat("Saved: annotation_check_figure_genes.csv\n")


section("SESSION INFORMATION")
print(sessionInfo())
cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

sink(type = "message")
sink()
close(log_con)
