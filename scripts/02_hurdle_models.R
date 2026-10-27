############################################################
# 02_hurdle_models.R
#
# Microbial associations with respiratory and unspecific clinical
# status across age groups (<25 / >=25 days).
#
# Two-part (hurdle) mixed models per taxon, farm as random intercept:
#   Detection          : glmer, present ~ Status + (1 | Farm_Number),
#                        binomial (Wald CIs and p-values)
#   Positive abundance : lmer (lmerTest), log1p(value) ~ Status + (1 | Farm_Number),
#                        maximum likelihood, fitted to samples above the
#                        detection cutoff (Satterthwaite p-values)
# Comparisons, fitted separately within each age group:
#   respiratory vs healthy
#   unspecific  vs healthy
#
# Samples
#   Nasal swabs ("_N_" in the sample ID) from healthy, respiratory and
#   unspecific calves with a single-value age in days (Age_days in
#   metadata.csv; ranges such as "15-42" are NA there and not used).
#
# Taxa
#   Viruses : target viral species, virus_targets (normalised abundance;
#             present = value > 5)
#   Bacteria: genera detected (> 3 reads per 6e7 paired reads) in
#             >= 18 sequenced calves (the genera in taxa_long.csv; the script
#             checks each target against this rule) (present = value > 3)
#
# Eligibility per taxon x age group x comparison
#   >= 10 samples, >= 5 samples per group, >= 5 positives in total,
#   >= 2 positives per group (detection and abundance models).
#
# Multiple testing
#   BH-adjusted p-values within each layer x model part x comparison
#   (across taxa and age groups), column p_adj_BH.
#
# Model QC (columns in the result tables)
#   singular, converged, fit_warnings, ci_width, wide_ci (CI width >= 20)
#
# Figure 3: forest plots of both layers. Bovine picornavirus and
#   estimates with CI width >= 20 are not shown (see legend).
#
# Inputs (input_files/)
#   metadata.csv    Status, Age_days, Farm_Number, Total_paired_reads
#   taxa_long.csv   viral normalised abundance, bacterial read counts
#
# Outputs (02_hurdle_models/)
#   log.txt                                   console log + session information
#   {virus,bacteria}_age_stratified_detection_results.csv
#   {virus,bacteria}_age_stratified_abundance_results.csv
#   {virus,bacteria}_age_stratified_feature_qc.csv
#   {virus,bacteria}_age_stratified_group_qc.csv
#   {virus,bacteria}_age_stratified_fit_log.csv   every fit / skip with reason
#   bacteria_target_check.csv                 prevalence of each target genus
#   Fig3_hurdle_{respiratory,unspecific}_vs_healthy.png
############################################################


############################################################
# 0 — SETUP
############################################################

# Run from the repository root: all paths below are relative to it.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(tidytext)
})

out_dir <- "02_hurdle_models"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

in_dir <- "input_files"

# ---- parameters ----
virus_targets <- c(
  "Bovine adenovirus 3", "Bovine coronavirus", "Bovine picornavirus",
  "Bovine rhinitis A", "Bovine rhinitis B", "Bovine torovirus", "BVDV 1",
  "Influenza D virus"
)

# genera detected (> 3 reads per 6e7 paired reads) in >= 18 sequenced calves
bacteria_targets <- c(
  "g__Mesomycoplasma", "g__Moraxella", "g__Mannheimia",
  "g__Mycoplasmopsis_A", "g__Moraxella_A", "g__Pasteurella", "g__Alloprevotella",
  "g__Streptococcus", "g__Prevotella", "g__Burkholderia", "g__Psychrobacter"
)
bact_min_calves <- 18

virus_present_cutoff <- 5     # normalised viral abundance
bact_present_cutoff  <- 3     # bacterial reads per 6e7 paired reads
bact_reads_scale     <- 6e7

min_total_samples        <- 10
min_group_n              <- 5
min_total_positive       <- 5
min_group_positive       <- 2
min_group_positive_abund <- 2

wide_ci_limit <- 20           # CI width flagged (and hidden in Figure 3)

age_groups  <- c("<25", ">=25")
comparisons <- list(c("healthy", "respiratory"),
                    c("healthy", "unspecific"))

# ---- log: console output and messages/warnings go to log.txt ----
# If the script stops with an error, run  sink(); sink(type = "message")
log_con <- file(file.path(out_dir, "log.txt"), open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

options(width = 250, readr.show_progress = FALSE)
cat("02_hurdle_models.R\n")
cat("Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Virus targets:", paste(virus_targets, collapse = ", "), "\n")
cat("Bacteria targets:", paste(bacteria_targets, collapse = ", "), "\n")
cat("Present: virus >", virus_present_cutoff, "| bacteria >", bact_present_cutoff, "reads per",
    format(bact_reads_scale, scientific = TRUE), "paired reads\n")
cat("Eligibility: >=", min_total_samples, "samples, >=", min_group_n, "per group, >=", min_total_positive,
    "positives, >=", min_group_positive, "positives per group\n")


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

# taxa x samples matrix from taxa_long.csv (bacterial genera given with "g__")
taxa_matrix <- function(taxa, layer) {
  d <- taxa_long %>% filter(type == layer, taxon %in% str_remove(taxa, "^g__"))
  missing <- setdiff(str_remove(taxa, "^g__"), d$taxon)
  if (length(missing) > 0) fail("Taxa absent from taxa_long.csv (", layer, "): ", paste(missing, collapse = ", "))
  m <- d %>% dplyr::select(taxon, Sample_id, value) %>%
    pivot_wider(names_from = Sample_id, values_from = value) %>%
    column_to_rownames("taxon") %>% as.matrix()
  m <- m[str_remove(taxa, "^g__"), , drop = FALSE]
  rownames(m) <- taxa
  m
}

safe_bind_rows <- function(x) {
  if (length(x) == 0) return(tibble())
  bind_rows(x)
}

# evaluate a model call, capturing warnings instead of printing them
fit_capture <- function(expr) {
  warns <- character(0)
  fit <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warnings = warns)
}

# convergence / singularity summary of a captured fit
model_qc <- function(res) {
  if (inherits(res$fit, "error")) {
    return(list(error = TRUE, message = conditionMessage(res$fit)))
  }
  conv <- res$fit@optinfo$conv$lme4$messages
  msgs <- unique(c(res$warnings, conv))
  list(
    error     = FALSE,
    singular  = lme4::isSingular(res$fit),
    converged = !any(grepl("converge", msgs, ignore.case = TRUE)),
    warnings  = if (length(msgs) > 0) paste(msgs, collapse = "; ") else NA_character_
  )
}

add_bh <- function(df) {
  if (nrow(df) == 0) return(df)
  df %>%
    group_by(Comparison) %>%
    mutate(p_adj_BH = p.adjust(p.value, method = "BH")) %>%
    ungroup()
}


############################################################
# 2 — HURDLE MODELS FOR ONE LAYER
############################################################

run_hurdle_models <- function(long_df, meta2, feature_col, present_cutoff) {
  
  results_detection <- list()
  results_abundance <- list()
  qc_feature_checks <- list()
  qc_group_counts   <- list()
  fit_log           <- list()
  
  log_fit <- function(feat, age_cat, comparison_label, model, outcome) {
    fit_log[[length(fit_log) + 1]] <<- tibble(
      feature = feat, Age_group = age_cat, Comparison = comparison_label,
      Model = model, outcome = outcome
    )
  }
  
  for (age_cat in age_groups) {
    for (comp in comparisons) {
      
      reference_group  <- comp[1]
      case_group       <- comp[2]
      comparison_label <- paste(case_group, "vs", reference_group, sep = "_")
      target_term      <- paste0("Status", case_group)
      
      # ---- comparison-level sample check ----
      meta_sub <- meta2 %>%
        filter(Age_category == age_cat,
               Status %in% c(reference_group, case_group)) %>%
        mutate(Status = factor(Status, levels = c(reference_group, case_group)))
      
      group_counts_meta <- meta_sub %>%
        count(Status, name = "n_samples") %>%
        mutate(Age_group = age_cat, Comparison = comparison_label)
      
      qc_group_counts[[length(qc_group_counts) + 1]] <- group_counts_meta
      
      if (nrow(group_counts_meta) < 2 || any(group_counts_meta$n_samples < min_group_n)) {
        log_fit("all", age_cat, comparison_label, "all", "skipped: insufficient samples for comparison")
        next
      }
      
      df_sub <- long_df %>%
        filter(Age_category == age_cat,
               Status %in% c(reference_group, case_group)) %>%
        mutate(Status = factor(Status, levels = c(reference_group, case_group)),
               Farm_Number = factor(Farm_Number))
      
      for (feat in unique(df_sub[[feature_col]])) {
        
        df <- df_sub %>%
          filter(.data[[feature_col]] == feat) %>%
          droplevels()
        
        # ---- per-feature QC ----
        sample_counts <- df %>% count(Status, name = "n_samples")
        
        pos_counts <- df %>%
          group_by(Status) %>%
          summarise(n_positive = sum(present, na.rm = TRUE),
                    n_zero     = sum(present == 0, na.rm = TRUE),
                    .groups = "drop")
        
        qc_feature_checks[[length(qc_feature_checks) + 1]] <- sample_counts %>%
          full_join(pos_counts, by = "Status") %>%
          mutate(!!feature_col := feat, Age_group = age_cat, Comparison = comparison_label) %>%
          dplyr::select(all_of(c(feature_col, "Age_group", "Comparison", "Status",
                                 "n_samples", "n_positive", "n_zero")))
        
        # ---- eligibility ----
        reason <- NULL
        if (nrow(sample_counts) < 2) {
          reason <- "one status group absent"
        } else if (any(sample_counts$n_samples < min_group_n)) {
          reason <- paste("fewer than", min_group_n, "samples in a group")
        } else if (nrow(df) < min_total_samples) {
          reason <- paste("fewer than", min_total_samples, "samples")
        } else if (sum(df$present, na.rm = TRUE) < min_total_positive) {
          reason <- paste("fewer than", min_total_positive, "positives")
        } else if (nrow(pos_counts) < 2 || any(pos_counts$n_positive < min_group_positive)) {
          reason <- paste("fewer than", min_group_positive, "positives in a group")
        }
        if (!is.null(reason)) {
          log_fit(feat, age_cat, comparison_label, "Detection + abundance", paste("skipped:", reason))
          next
        }
        
        # ---- detection model ----
        det <- fit_capture(
          glmer(present ~ Status + (1 | Farm_Number),
                data = df,
                family = binomial(),
                control = glmerControl(optimizer = "bobyqa",
                                       optCtrl = list(maxfun = 2e5)))
        )
        det_qc <- model_qc(det)
        
        if (det_qc$error) {
          log_fit(feat, age_cat, comparison_label, "Detection", paste("failed:", det_qc$message))
        } else {
          det_res <- broom.mixed::tidy(det$fit, effects = "fixed", conf.int = TRUE) %>%
            filter(term == target_term) %>%
            mutate(
              !!feature_col := feat,
              Model = "Detection",
              Age_group = age_cat,
              Comparison = comparison_label,
              Reference_group = reference_group,
              Case_group = case_group,
              OR = exp(estimate),
              OR_low = exp(conf.low),
              OR_high = exp(conf.high),
              Direction = ifelse(estimate > 0,
                                 paste0(case_group, "-associated"),
                                 paste0(reference_group, "-associated")),
              present_cutoff = present_cutoff,
              min_group_n = min_group_n,
              min_total_positive = min_total_positive,
              min_group_positive = min_group_positive,
              singular = det_qc$singular,
              converged = det_qc$converged,
              fit_warnings = det_qc$warnings,
              ci_width = conf.high - conf.low,
              wide_ci = ci_width >= wide_ci_limit
            )
          if (nrow(det_res) == 1) {
            results_detection[[length(results_detection) + 1]] <- det_res
            log_fit(feat, age_cat, comparison_label, "Detection",
                    ifelse(is.na(det_qc$warnings), "fitted", "fitted with warnings"))
          } else {
            log_fit(feat, age_cat, comparison_label, "Detection", "failed: target term not found")
          }
        }
        
        # ---- positive-abundance model (samples above cutoff) ----
        df_pos <- df %>%
          filter(present == 1) %>%
          droplevels()
        
        pos_group_counts <- df_pos %>% count(Status, name = "n_positive_samples")
        
        if (nrow(pos_group_counts) < 2 ||
            any(pos_group_counts$n_positive_samples < min_group_positive_abund)) {
          log_fit(feat, age_cat, comparison_label, "Positive_abundance",
                  paste("skipped: fewer than", min_group_positive_abund, "positives in a group"))
          next
        }
        
        pos <- fit_capture(
          lmer(log_value ~ Status + (1 | Farm_Number),
               data = df_pos,
               REML = FALSE)
        )
        pos_qc <- model_qc(pos)
        
        if (pos_qc$error) {
          log_fit(feat, age_cat, comparison_label, "Positive_abundance", paste("failed:", pos_qc$message))
        } else {
          pos_res <- broom.mixed::tidy(pos$fit, effects = "fixed", conf.int = TRUE) %>%
            filter(term == target_term) %>%
            mutate(
              !!feature_col := feat,
              Model = "Positive_abundance",
              Age_group = age_cat,
              Comparison = comparison_label,
              Reference_group = reference_group,
              Case_group = case_group,
              Fold_change = exp(estimate),
              FC_low = exp(conf.low),
              FC_high = exp(conf.high),
              Direction = ifelse(estimate > 0,
                                 paste0(case_group, "-associated"),
                                 paste0(reference_group, "-associated")),
              present_cutoff = present_cutoff,
              min_group_n = min_group_n,
              min_total_positive = min_total_positive,
              min_group_positive_abund = min_group_positive_abund,
              singular = pos_qc$singular,
              converged = pos_qc$converged,
              fit_warnings = pos_qc$warnings,
              ci_width = conf.high - conf.low,
              wide_ci = ci_width >= wide_ci_limit
            )
          if (nrow(pos_res) == 1) {
            results_abundance[[length(results_abundance) + 1]] <- pos_res
            log_fit(feat, age_cat, comparison_label, "Positive_abundance",
                    ifelse(is.na(pos_qc$warnings), "fitted", "fitted with warnings"))
          } else {
            log_fit(feat, age_cat, comparison_label, "Positive_abundance", "failed: target term not found")
          }
        }
      }
    }
  }
  
  list(
    detection  = add_bh(safe_bind_rows(results_detection)),
    abundance  = add_bh(safe_bind_rows(results_abundance)),
    feature_qc = safe_bind_rows(qc_feature_checks),
    group_qc   = safe_bind_rows(qc_group_counts),
    fit_log    = safe_bind_rows(fit_log)
  )
}

save_hurdle_results <- function(res, output_tag, feature_col) {
  write.csv(res$detection,  file.path(out_dir, paste0(output_tag, "_detection_results.csv")), row.names = FALSE)
  write.csv(res$abundance,  file.path(out_dir, paste0(output_tag, "_abundance_results.csv")), row.names = FALSE)
  write.csv(res$feature_qc, file.path(out_dir, paste0(output_tag, "_feature_qc.csv")),        row.names = FALSE)
  write.csv(res$group_qc,   file.path(out_dir, paste0(output_tag, "_group_qc.csv")),          row.names = FALSE)
  write.csv(res$fit_log,    file.path(out_dir, paste0(output_tag, "_fit_log.csv")),           row.names = FALSE)
  
  cat("\n====================================\n")
  cat(" MODEL SUMMARY —", output_tag, "\n")
  cat("====================================\n")
  cat("Detection results:", nrow(res$detection),
      "| abundance results:", nrow(res$abundance), "\n")
  cat("\nFit outcomes:\n")
  print(res$fit_log %>% count(Model, outcome), n = Inf)
  
  both <- bind_rows(res$detection, res$abundance)
  if (nrow(both) > 0) {
    cat("\nQC flags among reported estimates:\n")
    print(both %>% group_by(Model) %>%
            summarise(n = n(), singular = sum(singular), not_converged = sum(!converged),
                      wide_ci = sum(wide_ci), .groups = "drop"))
    cat("\nEstimates with p < 0.05 (with BH-adjusted p and QC flags):\n")
    print(both %>%
            filter(p.value < 0.05) %>%
            dplyr::select(all_of(feature_col), Model, Age_group, Comparison,
                          estimate, conf.low, conf.high, p.value, p_adj_BH,
                          singular, converged, wide_ci) %>%
            arrange(Comparison, Age_group, Model, p.value),
          n = Inf, width = Inf)
  }
}


############################################################
# 3 — METADATA
# Healthy, respiratory and unspecific calves with a single-value age
# in days.
############################################################

metadata_seq <- read_csv(file.path(in_dir, "metadata.csv"), col_types = cols(
  Sample_id = col_character(), Status = col_character(), Age_days = col_double(),
  Age_group = col_character(), Farm_Number = col_character(), Total_paired_reads = col_double(),
  has_virus = col_logical(), has_bacteria = col_logical(), has_host = col_logical()))
taxa_long <- read_csv(file.path(in_dir, "taxa_long.csv"), col_types = cols(
  Sample_id = col_character(), taxon = col_character(), taxon_label = col_character(),
  type = col_character(), value = col_double()))

metadata2 <- metadata_seq %>%
  filter(str_detect(Sample_id, "_N_")) %>%
  mutate(Age_category = if_else(Age_days < 25, "<25", ">=25"))

cat("\nSequenced samples with metadata:", nrow(metadata_seq),
    "| nasal-labelled (_N_ in the sample ID):", nrow(metadata2), "\n")
cat("\nClinical status, nasal samples:\n")
metadata2 %>% count(Status) %>% print(n = Inf)
cat("\nNasal samples without a single-value age in days:", sum(is.na(metadata2$Age_days)), "\n")

metadata_filtered <- metadata2 %>%
  filter(Status %in% c("healthy", "respiratory", "unspecific"),
         !is.na(Age_days))

cat("\nCalves eligible (healthy/respiratory/unspecific with age in days):",
    nrow(metadata_filtered), "\n")
cat("Of these, farm missing (dropped by the models):",
    sum(is.na(metadata_filtered$Farm_Number)), "\n")


############################################################
# SECTION A — VIRAL SPECIES
############################################################

virus_mat <- taxa_matrix(virus_targets, "virus")

common_v  <- intersect(colnames(virus_mat), metadata_filtered$Sample_id)
virus_mat <- virus_mat[, common_v, drop = FALSE]

meta_virus <- metadata_filtered %>%
  filter(Sample_id %in% common_v) %>%
  arrange(match(Sample_id, colnames(virus_mat)))
stopifnot(identical(meta_virus$Sample_id, colnames(virus_mat)))

cat("\nEligible calves missing from the viral matrix:",
    length(setdiff(metadata_filtered$Sample_id, colnames(virus_mat))), "\n")

virus_long <- as.data.frame(virus_mat) %>%
  rownames_to_column("Virus") %>%
  pivot_longer(-Virus, names_to = "Sample_id", values_to = "value") %>%
  left_join(meta_virus, by = "Sample_id") %>%
  mutate(present     = as.integer(value > virus_present_cutoff),
         log_value   = log1p(value),
         Farm_Number = factor(Farm_Number))

cat("\nViruses:", nrow(virus_mat), "| calves:", ncol(virus_mat), "\n")
print(meta_virus %>% count(Age_category, Status))

res_virus <- run_hurdle_models(virus_long, meta_virus,
                               feature_col    = "Virus",
                               present_cutoff = virus_present_cutoff)
save_hurdle_results(res_virus, "virus_age_stratified", "Virus")


############################################################
# SECTION B — BACTERIAL GENERA
############################################################

bact_raw <- taxa_matrix(bacteria_targets, "bacteria")   # read counts, taxa x samples

# ---- check: each target genus is detected in >= bact_min_calves sequenced calves ----
seq_cols  <- intersect(colnames(bact_raw), metadata_seq$Sample_id)
reads_all <- metadata_seq$Total_paired_reads[match(seq_cols, metadata_seq$Sample_id)]
norm_all  <- sweep(bact_raw[, seq_cols, drop = FALSE], 2, reads_all, "/") * bact_reads_scale

target_check <- tibble(
  genus    = rownames(norm_all),
  n_calves = rowSums(norm_all > bact_present_cutoff)
) %>%
  mutate(pct_calves = round(100 * n_calves / length(seq_cols), 1),
         meets_rule = n_calves >= bact_min_calves) %>%
  arrange(desc(n_calves))
write_csv(target_check, file.path(out_dir, "bacteria_target_check.csv"))

cat("\nTarget genera: calves (of", length(seq_cols), "sequenced) with >", bact_present_cutoff,
    "reads per", format(bact_reads_scale, scientific = TRUE), "paired reads:\n")
print(target_check, n = Inf)
if (!all(target_check$meets_rule))
  fail("Target genera below the prevalence rule (>= ", bact_min_calves, " calves): ",
       paste(target_check$genus[!target_check$meets_rule], collapse = ", "))

# ---- model matrix: target genera, eligible calves, reads per 6e7 ----
common_b <- intersect(colnames(bact_raw), metadata_filtered$Sample_id)
reads_b  <- metadata_seq$Total_paired_reads[match(common_b, metadata_seq$Sample_id)]
bact_mat <- bact_raw[, common_b, drop = FALSE]

meta_bact <- metadata_filtered %>%
  filter(Sample_id %in% common_b) %>%
  arrange(match(Sample_id, colnames(bact_mat)))
stopifnot(identical(meta_bact$Sample_id, colnames(bact_mat)))

cat("Eligible calves missing from the bacterial matrix:",
    length(setdiff(metadata_filtered$Sample_id, colnames(bact_mat))), "\n")

bact_mat <- sweep(bact_mat, 2, reads_b, "/") * bact_reads_scale

bact_long <- as.data.frame(bact_mat) %>%
  rownames_to_column("BACT") %>%
  pivot_longer(-BACT, names_to = "Sample_id", values_to = "value") %>%
  left_join(meta_bact, by = "Sample_id") %>%
  mutate(present     = as.integer(value > bact_present_cutoff),
         log_value   = log1p(value),
         Farm_Number = factor(Farm_Number))

cat("\nGenera:", nrow(bact_mat), "| calves:", ncol(bact_mat), "\n")
print(meta_bact %>% count(Age_category, Status))

res_bact <- run_hurdle_models(bact_long, meta_bact,
                              feature_col    = "BACT",
                              present_cutoff = bact_present_cutoff)
save_hurdle_results(res_bact, "bacteria_age_stratified", "BACT")


############################################################
# SECTION C — FIGURE 3
# Forest plot of detection (log-odds) and positive-abundance
# (log fold change) estimates, by age group, both layers.
############################################################

features_to_exclude <- c("Bovine picornavirus")  # excluded from the figure only

load_hurdle <- function(file, type, model_part) {
  df <- read_csv(file, show_col_types = FALSE) %>%
    mutate(Layer = type, Part = model_part)
  if ("Virus" %in% colnames(df)) df <- dplyr::rename(df, Feature = Virus)
  if ("BACT"  %in% colnames(df)) df <- dplyr::rename(df, Feature = BACT)
  df %>% mutate(Feature = str_remove(Feature, "^g__"))
}

all_res <- bind_rows(
  load_hurdle(file.path(out_dir, "virus_age_stratified_detection_results.csv"),    "Virus",    "Detection"),
  load_hurdle(file.path(out_dir, "virus_age_stratified_abundance_results.csv"),    "Virus",    "Abundance"),
  load_hurdle(file.path(out_dir, "bacteria_age_stratified_detection_results.csv"), "Bacteria", "Detection"),
  load_hurdle(file.path(out_dir, "bacteria_age_stratified_abundance_results.csv"), "Bacteria", "Abundance")
)

for (target_comparison in c("respiratory_vs_healthy", "unspecific_vs_healthy")) {
  
  plot_data <- all_res %>%
    filter(Comparison == target_comparison) %>%
    filter(!Feature %in% features_to_exclude) %>%
    filter(!is.na(estimate) & !is.na(std.error)) %>%
    filter((conf.high - conf.low) < wide_ci_limit) %>%
    mutate(
      Part = case_when(
        Part %in% c("Abundance", "Positive_abundance") ~ "Positive Abundance (Log-FC)",
        Part == "Detection" ~ "Detection (Log-Odds)",
        TRUE ~ as.character(Part)
      ),
      Part = factor(Part, levels = c("Detection (Log-Odds)", "Positive Abundance (Log-FC)")),
      Feature_panel = reorder_within(Feature, estimate, interaction(Age_group, Part))
    )
  
  p <- ggplot(plot_data, aes(x = estimate, y = Feature_panel, color = Layer)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.2, linewidth = 0.5, orientation = "y") +
    geom_point(size = 2.5) +
    facet_wrap(Age_group ~ Part, scales = "free", ncol = 2) +
    scale_y_reordered() +
    coord_cartesian(xlim = c(-6, 6)) +
    scale_color_manual(values = c("Bacteria" = "#2A9D8F", "Virus" = "#E1A730")) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold"),
      legend.position = "top",
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(margin = ggplot2::margin(t = 10)),
      axis.text.y = element_text(size = 9, face = "italic")
    ) +
    labs(
      title = paste("Association with", str_replace(target_comparison, "_vs_healthy", "")),
      x = "Coefficient Estimate (Log-Odds / Log-FC)",
      y = NULL,
      color = "Microbial Layer"
    )
  
  print(p)
  ggsave(file.path(out_dir, paste0("Fig3_hurdle_", target_comparison, ".png")),
         p, width = 10, height = 8, dpi = 300)
}

cat("\nEstimates hidden from Figure 3 (CI width >=", wide_ci_limit, "or excluded taxon):\n")
all_res %>%
  filter(!is.na(estimate), (conf.high - conf.low) >= wide_ci_limit | Feature %in% features_to_exclude) %>%
  dplyr::select(Feature, Layer, Model, Age_group, Comparison, estimate, conf.low, conf.high) %>%
  print(n = Inf, width = Inf)


############################################################
# 4 — SESSION INFORMATION AND CLOSE LOG
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