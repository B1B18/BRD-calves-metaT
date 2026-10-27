############################################################
# 05a_gene_annotation.R
#
# Transparent, rule-based annotation of every host gene reported in the results.
# Descriptive only: no enrichment testing.
#
# Pipeline position: run after 03a (DE) and 04a (networks).
# 05b_figures.R reads its outputs.
#
# Genes annotated
#   all DE genes from 03a (FDR < fdr_cutoff, |log2FC| > lfc_cutoff, any test)
#   + the annot_resp_top_n genes with the lowest FDR in respiratory_vs_healthy
#     (supplementary heatmap in 05b; may include genes with |log2FC| <= lfc_cutoff)
#   + every host gene in a stable edge of the correlation networks (04a 
#     network_edge_files)
#
# Sources
#   Entrez Gene ID: org.Bt.eg.db, matched by the RefSeq symbol used for alignment
#     as (a) current symbol, (b) previous/alternative symbol (alias, unique
#     matches only), (c) LOC genes: the number after "LOC" (NCBI Gene ID);
#     (d) otherwise from BioMart. NCBI description and current symbol: org.Bt.eg.db
#   GO biological process: org.Bt.eg.db PLUS the GO BP of the human gene
#     (org.Hs.eg.db): the single BioMart orthologue or, if there is none, the human
#     gene with the identical symbol; the source of each term is kept in the
#     evidence. Term names and hierarchy: GO.db
#   Human orthologues: Ensembl BioMart release 115 (ARS-UCD1.3, last release on
#     the ARS-UCD1.x assembly), read from the deposited cache
#     input_files/input_annotations.rds; only genes not in the cache are queried
#     (query_ensembl = TRUE); if Ensembl is unreachable the cache is used and the
#     log says how many genes lack BioMart data.
#
# Gene labels (figures)
#   The NCBI (RefSeq) symbol. LOC genes with exactly one human orthologue are
#   shown as "LOC... (HUMAN)". Symbols are never replaced by inferred names.
#
# Functional category: ONE per gene, from the first route that applies, in
# this fixed order:
#   1 manual         manual_categories (expert curation; listed in the log,
#                    flagged when GO supports a different category)
#   2 gene family    family_rules: explicit patterns on the NCBI description,
#                    mainly for LOC genes of expanded families without GO
#                    annotation (e.g. "immunoglobulin lambda-1 light chain-like")
#   3 GO BP          the gene's GO BP annotation falls within the GO terms
#                    defining a category (category_go; descendant terms
#                    included); categories checked in the listed priority order
#   4 GO BP, other   GO BP annotation, none within the category terms
#                    -> "Other/Uncharacterised"
#   5 unannotated    no GO BP annotation and no family match
#                    -> "Other/Uncharacterised" (route 4 vs 5 kept in the evidence table)
#
# Inputs
#   03a_host_DE/03a_objects.rds
#   04a_networks_DE_genes/network_<status>_stable_bootstrap_edges.csv
#   input_files/input_annotations.rds   BioMart cache (Ensembl 115): human orthologues,
#                                       Entrez IDs; list(basic, entrez, genes, ...)
#
# Outputs (05a_gene_annotation/)
#   log.txt                        parameters, GO terms used (IDs and names),
#                                  genes per route/category, manual calls, conflicts
#   05a_annotation.rds             hand-off to 05b
#   gene_annotation_evidence.csv   every annotated gene: symbol, label, NCBI
#                                  description, Entrez ID, human orthologue(s),
#                                  category, route, evidence, all GO-supported
#                                  categories, conflict flag, DE / network membership
#                                  (supplementary table)
#   gene_annotation_table.csv      one row per annotated gene (label, category, lowest-FDR
#                                  DE test if DE, network membership)
#   host_DE_all_results_entrezid.csv
#   biomart_btaurus_ensembl115.rds copy of the updated BioMart cache, written ONLY if genes were
#                                  added; input_files/input_annotations.rds is updated at the same time
############################################################


############################################################
# 0 — PARAMETERS (see README)
############################################################

# Run from the repository root: all paths below are relative to it.

in_dir  <- "03a_host_DE"
out_dir <- "05a_gene_annotation"

# stable edges of the correlation networks (04a); every host gene in them is annotated.
# A missing file is reported in the log (run 04a first).
network_edge_files <- c(
  "04a" = "04a_networks_DE_genes/network_{status}_stable_bootstrap_edges.csv"
)
network_statuses <- c("healthy", "respiratory")

# extra genes from respiratory_vs_healthy (lowest FDR, any |log2FC|); 0 = none
# (Fig. 6A shows DE genes only, which are always annotated)
annot_resp_top_n <- 0

# BioMart
ensembl_host  <- "https://sep2025.ARCHIVE.ensembl.org"   # upper-case ARCHIVE: biomaRt 2.66 workaround
biomart_input <- file.path("input_files", "input_annotations.rds")         # deposited cache
biomart_cache <- file.path(out_dir, "biomart_btaurus_ensembl115.rds")       # record of the cache, written only if genes are added
query_ensembl <- TRUE      # FALSE = use the cache only (no internet)

# Route 1: manual curation (gene symbol = category)
manual_categories <- c(
  "PAPPA2"="ECM/Matrix/Proteolysis",
  "SULT1C2"="Metabolism",
  "ITLN2"="Innate Humoral Defense",
  "INHBA"="Cytokine/Inflammation",
  "PMP22"="Neuronal/Sensory/Reproductive",
  "GDF15"="Cytokine/Inflammation",
  "UNC5C"="Neuronal/Sensory/Reproductive",
  "MME"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "ACE2"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "VIPR1"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "PTGER2"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "LGALS1"="Innate Humoral Defense",
  "SPARC"="ECM/Matrix/Proteolysis",
  "ENO2"="Metabolism",
  "SUCNR1"="Innate Sensing/Inflammasome",
  "ALB"="Other/Uncharacterised",
  "SPARCL1"="ECM/Matrix/Proteolysis",
  "FUT5"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "MMEL1"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "MARCH10"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "RGS4"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "F2"="Coagulation/Haemostasis",
  "SELENOM"="Metabolism",
  "SELENBP1"="Metabolism",
  "AHNAK2"="Epithelial Barrier",
  "ANXA9"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "CYP4F2"="Metabolism",
  "BOSTAUV1R414"="Neuronal/Sensory/Reproductive",
  "F2R"="Coagulation/Haemostasis",
  "PCSK5"="ECM/Matrix/Proteolysis",
  "SPINK5"="Epithelial Barrier",
  "CES2"="Metabolism",
  "FBXO15"="Apoptosis/Stress/Remodeling",
  "UNC5CL"="Innate Sensing/Inflammasome",
  "CES4A"="Metabolism",
  "EMILIN2"="ECM/Matrix/Proteolysis",
  "SELENOK"="Apoptosis/Stress/Remodeling",
  "CAPN14"="Epithelial Barrier",
  "DMBT1"="Innate Humoral Defense",
  "SEPTIN4"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "GGT6"="Metabolism",
  "SESN2"="Apoptosis/Stress/Remodeling",
  "CYR61"="ECM/Matrix/Proteolysis",
  "TSPAN9"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "CPNE2"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "ARRB2"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "F13A1"="Coagulation/Haemostasis",
  "PDE4B"="Cytokine/Inflammation",
  "SIGLEC5"="Innate Sensing/Inflammasome",
  "OLFM2"="Neuronal/Sensory/Reproductive",
  "TMEM176B"="Innate Sensing/Inflammasome",
  "A2ML1"="Innate Humoral Defense",
  "WFDC5"="Innate Humoral Defense",
  "TIRAP"="Innate Sensing/Inflammasome",
  "ECM1"="ECM/Matrix/Proteolysis",
  "SIGLEC8"="Innate Sensing/Inflammasome",
  "SELPLG"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "TGM1"="Epithelial Barrier",
  "AGER"="Innate Sensing/Inflammasome",
  "FGL2"="Cytokine/Inflammation",
  "FCRL4"="Adaptive/Leukocyte",
  "GCA"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "KNG1"="Cytokine/Inflammation",
  "RN7SL1"="Other/Uncharacterised",
  "LSP1"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "INHBB"="Cytokine/Inflammation",
  "WARS"="Interferon/Antiviral",
  "TGM3"="Epithelial Barrier",
  "THBD"="Coagulation/Haemostasis",
  "FBXL16"="Apoptosis/Stress/Remodeling",
  "NMI"="Interferon/Antiviral",
  "SP140"="Innate Sensing/Inflammasome",
  "SPINK7"="Epithelial Barrier",
  "OAF"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "OSCAR"="Innate Sensing/Inflammasome",
  "PI3"="Innate Humoral Defense",
  "FUT7"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "PLAC8"="Innate Immunity (general)",
  "FBXO32-2"="Apoptosis/Stress/Remodeling",
  "ARG2"="Metabolism",
  "OLFM4"="Innate Humoral Defense",
  "SCIMP"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "ZFP36"="Cytokine/Inflammation",
  "LCN12"="Innate Humoral Defense",
  "RGS2-2"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "SNX20"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "SH3BP2"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "TGM2"="ECM/Matrix/Proteolysis",
  "SIRPB2"="Innate Sensing/Inflammasome",
  "THPO"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "SRGN"="Adaptive/Leukocyte",
  "F5"="Coagulation/Haemostasis",
  "BTN3A3"="Innate Sensing/Inflammasome",
  "EPO"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "FBLN1"="ECM/Matrix/Proteolysis",
  "FBXO39"="Apoptosis/Stress/Remodeling",
  "FCRL3"="Adaptive/Leukocyte",
  "CSTA"="Epithelial Barrier",
  "SLA"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "RETN"="Cytokine/Inflammation",
  "FBLN7"="ECM/Matrix/Proteolysis",
  "FBLL1"="ECM/Matrix/Proteolysis",
  "APLN"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "FBXO43"="Apoptosis/Stress/Remodeling",
  "PLA2G2A-2"="Innate Humoral Defense",
  "MS4A1"="Adaptive/Leukocyte",
  "A2M"="Innate Humoral Defense",
  "MS4A15"="Adaptive/Leukocyte",
  "IGCGAMMA"="Adaptive/Leukocyte",
  "FFAR2"="Innate Sensing/Inflammasome",
  "TFPI2"="ECM/Matrix/Proteolysis",
  "MS4A3"="Adaptive/Leukocyte",
  "GPR42"="Innate Sensing/Inflammasome",
  # genes shown in the figures
  "FRMPD2"="Epithelial Barrier",
  "GAS2L2"="Ciliary/Mucociliary",
  "PROK2"="Cytokine/Inflammation",
  "ERICH3"="Other/Uncharacterised",
  "WDR63"="Ciliary/Mucociliary",
  "ZMYND10"="Ciliary/Mucociliary",
  # reviewed figure genes (reasons below)
  "JCHAIN"="Adaptive/Leukocyte",
  "PTX3"="Innate Humoral Defense",
  "GIMAP7"="Adaptive/Leukocyte",
  "ZC3H12A"="Cytokine/Inflammation",
  "TNF"="Cytokine/Inflammation",
  "NLRP3"="Innate Sensing/Inflammasome",
  "MEFV"="Innate Sensing/Inflammasome",
  "FCAR"="Adaptive/Leukocyte",
  # loci without an NCBI/BioMart match (RefSeq "-N" suffix = additional locus sharing the symbol)
  "BOLA"="Adaptive/Leukocyte",
  "BOLA-2"="Adaptive/Leukocyte",
  "GIMAP7-2"="Adaptive/Leukocyte",
  "GIMAP7-3"="Adaptive/Leukocyte",
  "GIMAP7-4"="Adaptive/Leukocyte",
  "ITLN2-2"="Innate Humoral Defense",
  "APOL3"="Innate Humoral Defense",
  "APOL3-2"="Innate Humoral Defense",
  "PLAC8-2"="Innate Immunity (general)",
  "MT1E-2"="Apoptosis/Stress/Remodeling",
  "SERPINA5-2"="Coagulation/Haemostasis",
  "BLOC1S3-2"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "PTPN13-2"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "SSTR1-2"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "HIST1H2BI"="Transcription/RNA/Translation",
  "MNF1"="Metabolism",
  "MGC127055"="Metabolism",
  "VNN3"="Metabolism",
  "BEX2"="Apoptosis/Stress/Remodeling",
  "AGBL4"="Ciliary/Mucociliary",
  "BICD1"="Signaling/Trafficking/Cytoskeleton/Adhesion",
  "GUCY1A2"="Signaling/Trafficking/Cytoskeleton/Adhesion"
)


# reasons for manual calls made during the review of figure genes (shown in the evidence table)
manual_reasons <- c(
  "OLFM4"   = "olfactomedin 4, neutrophil granule protein",
  "ZMYND10" = "axonemal dynein-arm assembly factor; loss causes primary ciliary dyskinesia",
  "JCHAIN"  = "immunoglobulin J chain, required for polymeric IgA/IgM assembly",
  "PTX3"    = "pentraxin 3, soluble pattern-recognition molecule",
  "GIMAP7"  = "GTPase of immunity-associated protein family, expressed in lymphocytes",
  "ZC3H12A" = "Regnase-1, endoribonuclease degrading inflammatory cytokine mRNAs",
  "TNF"     = "tumour necrosis factor, pro-inflammatory cytokine",
  "NLRP3"   = "NLR family pyrin domain containing 3, inflammasome sensor",
  "MEFV"    = "pyrin, inflammasome sensor",
  "FCAR"    = "Fc alpha receptor (CD89), myeloid receptor for IgA",
  "UNC5C"   = "netrin-1 receptor; axon guidance and dependence-receptor apoptosis",
  "OLFM2"   = "olfactomedin 2, secreted glycoprotein of neuronal and retinal tissue",
  "ACE2"    = "angiotensin-converting enzyme 2, membrane carboxypeptidase in renin-angiotensin peptide signalling",
  "MME"     = "neprilysin, membrane metalloendopeptidase inactivating signalling peptides (substance P, bradykinin, natriuretic peptides)",
  "MMEL1"   = "membrane metalloendopeptidase-like 1, neprilysin-family peptidase",
  "KNG1"    = "kininogen 1, bradykinin precursor (kallikrein-kinin system)",
  "TGM2"    = "tissue transglutaminase; cross-links extracellular matrix proteins and binds fibronectin",
  "SRGN"    = "serglycin, proteoglycan of leukocyte secretory granules (mast cells, neutrophils, cytotoxic lymphocytes)",
  "ALB"     = "albumin; no specific host-defence role",
  "PLAC8"   = "placenta-associated 8 (onzin); antibacterial host defence in mouse (inferred); matched to PLAC8-2",
  "F2"      = "prothrombin, coagulation serine protease",
  "F5"      = "coagulation factor V, prothrombinase cofactor",
  "F13A1"   = "coagulation factor XIII A, transglutaminase cross-linking fibrin",
  "THBD"    = "thrombomodulin, thrombin cofactor for protein C activation",
  "F2R"     = "thrombin receptor (PAR1)",
  # loci without an NCBI/BioMart match
  "BOLA"       = "bovine MHC class I (BoLA), antigen presentation",
  "BOLA-2"     = "additional bovine MHC class I (BoLA) locus",
  "GIMAP7-2"   = "additional GIMAP7 locus (lymphocyte GTPase)",
  "GIMAP7-3"   = "additional GIMAP7 locus (lymphocyte GTPase)",
  "GIMAP7-4"   = "additional GIMAP7 locus (lymphocyte GTPase)",
  "ITLN2-2"    = "additional intelectin-2 locus (microbial glycan-binding lectin)",
  "APOL3"      = "apolipoprotein L3; antibacterial effector in human (inferred)",
  "APOL3-2"    = "additional apolipoprotein L3 locus",
  "PLAC8-2"    = "additional PLAC8 (onzin) locus; antibacterial host defence in mouse (inferred)",
  "MT1E-2"     = "additional metallothionein 1E locus, zinc binding and oxidative stress",
  "SERPINA5-2" = "additional protein C inhibitor (SERPINA5) locus; regulates the protein C anticoagulant pathway",
  "BLOC1S3-2"  = "additional BLOC-1 subunit 3 locus, endosomal trafficking",
  "PTPN13-2"   = "additional PTPN13 locus, non-receptor tyrosine phosphatase",
  "SSTR1-2"    = "additional somatostatin receptor 1 locus (GPCR)",
  "HIST1H2BI"  = "histone H2B, chromatin",
  "MNF1"       = "UQCC2 (formerly MNF1; C6orf125 homologue), mitochondrial complex III assembly",
  "MGC127055"  = "uncharacterised (LOC785824); orthologous to cytochrome P450 2D family (putative CYP2D43)",
  "VNN3"       = "vanin 3, pantetheinase",
  "BEX2"       = "brain-expressed X-linked 2, apoptosis and cell-cycle regulation",
  "AGBL4"      = "cytosolic carboxypeptidase 6 (CCP6), tubulin deglutamylation and cilium assembly",
  "BICD1"   = "BICD cargo adaptor 1, dynein-dynactin cargo adaptor in vesicle transport; GO 'RNA processing' is peripheral",
  "GUCY1A2" = "soluble guanylate cyclase alpha-2 subunit, nitric oxide-cGMP signalling"
)


# genes whose category, route and evidence are printed in the log for checking (character(0) = none)
check_genes <- c("JCHAIN", "PTX3", "GIMAP7", "ZC3H12A", "RNF40", "EIF3G", "VLDLR", "NWD1", "WAS", "FCAR",
                 "TMOD2", "TUBB1", "PCDH20", "ICAM4",
                 "UNC5C", "OLFM2", "ACE2", "MME", "MMEL1", "KNG1", "TGM2", "SRGN", "ALB", "PLAC8",
                 "F2", "F5", "F13A1", "THBD", "F2R", "SERPINA5-2")


# Route 2: gene-family rules on the NCBI description (lower case); first match wins
family_rules <- tibble::tribble(
  ~rule,                       ~pattern,                                                                                               ~category,
  "immunoglobulin chain",      "immunoglobulin (heavy|lambda|kappa|light|iota|j chain)|\\big (heavy|lambda|kappa|mu\\b)",             "Adaptive/Leukocyte",
  "granzyme / perforin",       "\\bgranzyme|\\bperforin",                                                                              "Adaptive/Leukocyte",
  "MHC / BoLA",                "major histocompatibility complex|\\bbola\\b",                                                          "Adaptive/Leukocyte",
  "T-cell receptor chain",     "t[- ]cell receptor (alpha|beta|gamma|delta)",                                                          "Adaptive/Leukocyte",
  "interferon",                "^interferon (alpha|beta|gamma|lambda|tau|omega|epsilon|kappa)",                                        "Interferon/Antiviral",
  "interferon-induced",        "interferon[- ]induced|interferon[- ]stimulated|2'-5'-oligoadenylate|guanylate[- ]binding protein|\\bmx dynamin|radical s-adenosyl methionine domain containing 2", "Interferon/Antiviral",
  "LILR receptor",             "leukocyte immunoglobulin-like receptor",                                                               "Cytokine/Inflammation",
  "toll-like receptor",        "toll-like receptor",                                                                                   "Innate Sensing/Inflammasome",
  "intelectin",                "\\bintelectin",                                                                                        "Innate Humoral Defense",
  "interleukin / chemokine",   "\\binterleukin|\\bchemokine|c-c motif|c-x-c motif",                                                    "Cytokine/Inflammation",
  "complement",                "\\bcomplement (c[0-9]|factor|component)",                                                              "Innate Humoral Defense",
  "antimicrobial / acute phase", "defensin|antimicrobial peptide|cathelicidin|lysozyme|lactotransferrin|bactericidal|lipocalin|haptoglobin|serum amyloid|ceruloplasmin|peptidoglycan recognition", "Innate Humoral Defense",
  "ciliary / axonemal",        "\\bcilia|\\bciliary|dynein axonemal|radial spoke|\\btektin|cilia and flagella associated|dynein regulatory complex|flagell|axonemal",            "Ciliary/Mucociliary",
  "keratin / mucin / junction", "\\bkeratin|\\bmucin\\b|mucin-|\\bclaudin|small proline[- ]rich|\\bdesmo(glein|collin|plakin)",       "Epithelial Barrier",
  "matrix / protease",         "matrix metallopeptidase|metalloproteinase|\\bcollagen\\b|collagen type|\\bserpin",                     "ECM/Matrix/Proteolysis",
  "solute carrier",            "solute carrier",                                                                                       "Ion/Solute Transport",
  "olfactory / taste receptor", "olfactory receptor|taste receptor",                                                                   "Neuronal/Sensory/Reproductive",
  "ribosomal / histone / ZNF", "ribosomal protein|\\bhistone h|zinc finger protein [0-9]",                                                   "Transcription/RNA/Translation"
)

# Route 3: GO BP terms defining each category (descendants included).
# Order = priority when a gene's GO annotation supports several categories.
category_go <- list(
  "Interferon/Antiviral"          = c("GO:0051607", "GO:0034340", "GO:0060337", "GO:0034341"),
  "Adaptive/Leukocyte"            = c("GO:0002250", "GO:0042110", "GO:0042113", "GO:0019882", "GO:0002449"),
  "Innate Sensing/Inflammasome"   = c("GO:0002221", "GO:0070269"),
  "Cytokine/Inflammation"         = c("GO:0006954", "GO:0019221", "GO:0006935", "GO:0001816"),
  "Innate Humoral Defense"        = c("GO:0006956", "GO:0019730", "GO:0042742", "GO:0006953"),
  "Autophagy/Lysosomal Defense"   = c("GO:0006914", "GO:0007040"),
  "Ciliary/Mucociliary"           = c("GO:0044782", "GO:0003341"),
  "Epithelial Barrier"            = c("GO:0008544", "GO:0045216", "GO:0002064"),
  "ECM/Matrix/Proteolysis"        = c("GO:0030198", "GO:0022617"),   # ECM organisation, ECM disassembly
  "Apoptosis/Stress/Remodeling"   = c("GO:0006915", "GO:0034976", "GO:0006986"),
  "Ion/Solute Transport"          = c("GO:0006811", "GO:0055085"),
  "Metabolism"                    = c("GO:0044281", "GO:0006629", "GO:0006091"),
  "Transcription/RNA/Translation" = c("GO:0006355", "GO:0006412", "GO:0006396"),
  "Neuronal/Sensory/Reproductive" = c("GO:0007399", "GO:0007606", "GO:0022414"),
  "Coagulation/Haemostasis"       = c("GO:0007596"),   # blood coagulation; just above Signaling, so it only takes genes that would otherwise fall to Signaling or the general categories
  "Signaling/Trafficking/Cytoskeleton/Adhesion" = c("GO:0007165", "GO:0016192", "GO:0007010", "GO:0007155"),
  # lowest priority: only genes with no more specific category term
  "Innate Immunity (general)"     = c("GO:0045087"),
  "Immune Response (general)"     = c("GO:0006955")
)


############################################################
# 1 — SETUP, LOG AND HELPERS
############################################################

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Bt.eg.db)
  library(org.Hs.eg.db)
  library(GO.db)
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
    if (is.data.frame(val)) {
      cat(nm, ":\n"); print(as.data.frame(val), row.names = FALSE)
    } else if (is.list(val)) {
      cat(nm, ":\n"); for (k in names(val)) cat(sprintf("  %-32s %s\n", k, paste(val[[k]], collapse = ", ")))
    } else if (!is.null(names(val))) {
      cat(nm, "(", length(val), "):\n"); cat(sprintf("  %-14s %s\n", names(val), val), sep = "")
    } else {
      cat(sprintf("%-22s = %s\n", nm, if (length(val) == 0) "(none)" else paste(val, collapse = ", ")))
    }
  }
}

cat("05a_gene_annotation.R\nRun:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

in_rds <- file.path(in_dir, "03a_objects.rds")
if (!file.exists(in_rds)) fail("Run 03a_host_DE.R first: ", in_rds, " not found.")
a <- readRDS(in_rds)
fdr_cutoff <- a$params$fdr_cutoff
lfc_cutoff <- a$params$lfc_cutoff
lfc_cutoff_high <- if (is.null(a$params$lfc_cutoff_high)) 2 else a$params$lfc_cutoff_high
cat("Input:", in_rds, "(created", format(a$created), ")\n")
cat("From 03a: fdr_cutoff =", fdr_cutoff, "| lfc_cutoff =", lfc_cutoff, "\n")

log_parameters(c("network_edge_files", "network_statuses", "annot_resp_top_n", "ensembl_host", "biomart_input", "query_ensembl",
                 "manual_categories", "manual_reasons", "family_rules", "category_go", "check_genes"))


############################################################
# 2 — GENES TO ANNOTATE
############################################################

section("GENES TO ANNOTATE")
DE_sig <- a$DE_all %>% filter(FDR < fdr_cutoff, abs(logFC) > lfc_cutoff)
resp_top <- a$DE_all %>%
  filter(comparison == "respiratory_vs_healthy", FDR < fdr_cutoff) %>%
  arrange(FDR) %>% slice_head(n = annot_resp_top_n)
de_genes    <- unique(str_remove(DE_sig$gene, "^gene-"))
resp_genes  <- unique(str_remove(resp_top$gene, "^gene-"))
high_genes  <- unique(str_remove(DE_sig$gene[abs(DE_sig$logFC) > lfc_cutoff_high], "^gene-"))
cat("DE genes with |log2FC| >", lfc_cutoff_high, ":", length(high_genes), "\n")

# host genes in the stable network edges (columns from/to; host genes start with "gene-")
net_genes <- bind_rows(lapply(names(network_edge_files), function(src) {
  bind_rows(lapply(network_statuses, function(st) {
    f <- str_replace(network_edge_files[[src]], fixed("{status}"), st)
    if (!file.exists(f)) { cat("WARNING: network edge file not found:", f, "\n"); return(NULL) }
    e <- read_csv(f, show_col_types = FALSE, col_types = cols(.default = col_character()))
    if (nrow(e) == 0) return(NULL)
    tibble(id = c(e$from, e$to)) %>% filter(str_detect(id, "^gene-")) %>%
      distinct() %>% transmute(gene = str_remove(id, "^gene-"), network = paste0(src, " ", st))
  }))
}))
net_membership <- if (nrow(net_genes) > 0) {
  net_genes %>% group_by(gene) %>% summarise(networks = paste(sort(unique(network)), collapse = "; "), .groups = "drop")
} else tibble(gene = character(), networks = character())
network_genes <- net_membership$gene

annot_genes <- union(union(de_genes, resp_genes), network_genes)
cat("DE genes:", length(de_genes),
    "| added from respiratory_vs_healthy top", annot_resp_top_n, ":", length(setdiff(resp_genes, de_genes)),
    "| network genes:", length(network_genes), "(not DE:", length(setdiff(network_genes, de_genes)), ")",
    "| total:", length(annot_genes), "\n")
if (nrow(net_genes) > 0) {
  cat("Network genes per network:\n")
  net_genes %>% count(network) %>% as.data.frame() %>% print(row.names = FALSE)
}


############################################################
# 3 — ENSEMBL BIOMART: HUMAN ORTHOLOGUES AND ENTREZ IDS (cached)
############################################################

section("ENSEMBL BIOMART")
attr_basic  <- c("external_gene_name", "description", "hsapiens_homolog_associated_gene_name")
attr_entrez <- c("external_gene_name", "entrezgene_id")

# fallback: the BioMart web service on the same host, in batches (base R only)
biomart_rest <- function(attributes, genes, host = ensembl_host, batch_size = 200) {
  run_batch <- function(gene_batch) {
    xml <- paste0(
      '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>',
      '<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" datasetConfigVersion="0.6">',
      '<Dataset name="btaurus_gene_ensembl" interface="default">',
      '<Filter name="external_gene_name" value="', paste(gene_batch, collapse = ","), '"/>',
      paste0('<Attribute name="', attributes, '"/>', collapse = ""),
      '</Dataset></Query>'
    )
    query_url <- paste0(host, "/biomart/martservice?query=", URLencode(xml, reserved = TRUE))
    txt <- NULL
    for (attempt in 1:3) {
      txt <- tryCatch({
        con <- suppressWarnings(url(query_url, open = "r"))
        on.exit(close(con), add = TRUE)
        paste(readLines(con, warn = FALSE), collapse = "\n")
      }, error = function(e) { message("Attempt ", attempt, " failed: ", conditionMessage(e)); NULL })
      if (!is.null(txt)) break
      Sys.sleep(5 * attempt)
    }
    # stop(), not fail(): fail() closes the log, and this error is caught below
    if (is.null(txt)) stop("BioMart query failed after 3 attempts", call. = FALSE)
    if (grepl("ERROR", substr(txt, 1, 300))) stop("BioMart error: ", substr(txt, 1, 300), call. = FALSE)
    if (!nzchar(trimws(txt))) return(as_tibble(setNames(rep(list(character()), length(attributes)), attributes)))
    read_tsv(I(txt), col_names = attributes, na = c("", "NA"), col_types = cols(.default = col_character()))
  }
  batches <- split(genes, ceiling(seq_along(genes) / batch_size))
  cat("Querying the BioMart web service in", length(batches), "batches...\n")
  bind_rows(lapply(batches, run_batch)) %>% distinct()
}

query_biomart <- function(genes) {
  tryCatch({
    cow_mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "btaurus_gene_ensembl",
                                 host = ensembl_host)
    list(basic  = biomaRt::getBM(attributes = attr_basic,  filters = "external_gene_name",
                                 values = genes, mart = cow_mart, useCache = FALSE),
         entrez = biomaRt::getBM(attributes = attr_entrez, filters = "external_gene_name",
                                 values = genes, mart = cow_mart, useCache = FALSE),
         method = "biomaRt")
  }, error = function(e) {
    message("biomaRt failed (", conditionMessage(e), "); trying the BioMart web service directly.")
    tryCatch(list(basic  = biomart_rest(attr_basic,  genes),
                  entrez = biomart_rest(attr_entrez, genes),
                  method = "BioMart web service"),
             error = function(e2) NULL)
  })
}

empty_bm <- list(basic  = as_tibble(setNames(rep(list(character()), length(attr_basic)),  attr_basic)),
                 entrez = as_tibble(setNames(rep(list(character()), length(attr_entrez)), attr_entrez)),
                 method = NA_character_, release = NA_character_, host = ensembl_host,
                 date = as.Date(NA), genes = character(0))
if (file.exists(biomart_input)) {
  bm <- readRDS(biomart_input)
  if (!is.list(bm) || !all(c("basic", "entrez", "genes") %in% names(bm)))
    fail(biomart_input, " is not a BioMart cache (expected a list with basic, entrez and genes). ",
         "It must be the BioMart cache file (biomart_btaurus_ensembl115.rds), not an annotation output.")
  cat("BioMart cache:", biomart_input, "| method:", bm$method, "| release:", bm$release,
      "| date:", format(bm$date), "\n")
} else {
  cat("WARNING: BioMart cache not found:", biomart_input, "— all genes are queried.\n")
  bm <- empty_bm
}
new_genes <- setdiff(annot_genes, bm$genes)
cat("Cache covers", length(intersect(annot_genes, bm$genes)), "of",
    length(annot_genes), "genes\n")

if (length(new_genes) > 0 && query_ensembl) {
  cat("Querying", length(new_genes), "genes from", ensembl_host, "\n")
  release <- tryCatch(biomaRt::listMarts(host = ensembl_host)$version[1], error = function(e) "unknown")
  res <- query_biomart(new_genes)
  if (is.null(res)) {
    cat("WARNING: Ensembl BioMart could not be reached;", length(new_genes),
        "genes have no BioMart orthologue/Entrez ID in this run. Re-run later to complete the cache.\n")
  } else {
    to_chr <- function(d) mutate(as_tibble(d), across(everything(), as.character))
    bm$basic   <- bind_rows(to_chr(bm$basic),  to_chr(res$basic))  %>% distinct()
    bm$entrez  <- bind_rows(to_chr(bm$entrez), to_chr(res$entrez)) %>% distinct()
    bm$method  <- paste(na.omit(unique(c(bm$method, res$method))), collapse = " + ")
    bm$release <- if (!is.na(bm$release) && bm$release != "unknown") bm$release else release
    bm$host    <- ensembl_host
    bm$date    <- Sys.Date()
    bm$genes   <- union(bm$genes, new_genes)
    # the deposited input is updated so that it covers every annotated gene;
    # a dated copy of the updated cache is kept in the output folder
    dir.create(dirname(biomart_input), showWarnings = FALSE, recursive = TRUE)
    saveRDS(bm, biomart_input)
    saveRDS(bm, biomart_cache)
    cat("Added", length(new_genes), "genes to the cache via", res$method, "\n")
    cat("Input cache updated:", biomart_input, "| copy:", biomart_cache, "\n")
  }
} else if (length(new_genes) > 0) {
  cat("query_ensembl = FALSE:", length(new_genes), "genes have no BioMart data in this run.\n")
}
genes_without_biomart <- setdiff(annot_genes, bm$genes)

bm_basic  <- as.data.frame(bm$basic)  %>% mutate(across(everything(), as.character))
bm_entrez <- as.data.frame(bm$entrez) %>% mutate(across(everything(), as.character))


############################################################
# 4 — NCBI DESCRIPTIONS, ENTREZ IDS, GO BP, ORTHOLOGUES
############################################################

section("NCBI DESCRIPTIONS, GO AND ORTHOLOGUES")
# AnnotationDbi::select() stops when none of the keys is in the database: return an empty table
# instead (only for that error; any other error stops the script)
safe_select <- function(db, keys, columns, keytype) {
  keys <- unique(na.omit(as.character(keys)))
  empty <- as_tibble(setNames(rep(list(character()), length(c(keytype, columns))), c(keytype, columns)))
  if (length(keys) == 0) return(empty)
  res <- tryCatch(suppressMessages(AnnotationDbi::select(db, keys = keys, columns = columns, keytype = keytype)),
                  error = function(e) if (grepl("valid keys", conditionMessage(e))) NULL else stop(e))
  if (is.null(res)) empty else as_tibble(res) %>% mutate(across(everything(), as.character))
}

# Entrez Gene ID, first match in this order:
#   Only unique matches are used at each step (a symbol matching more than one Entrez ID
#   falls through to the next step; listed in the log).
#   (a) RefSeq symbol is the current NCBI symbol         (org.Bt.eg.db SYMBOL)
#   (b) RefSeq symbol is a previous/alternative symbol   (org.Bt.eg.db ALIAS; only aliases of one gene)
#   (c) LOC genes: the number after "LOC" is the NCBI Gene ID (kept only if in org.Bt.eg.db)
#   (d) Entrez ID from BioMart
ambiguous_ids <- list()
unique_ids <- function(d, source) {   # keep genes with exactly one Entrez ID; record the others
  amb <- d %>% group_by(gene) %>% filter(n_distinct(ENTREZID) > 1) %>%
    summarise(entrez_ids = paste(sort(unique(ENTREZID)), collapse = ";"), .groups = "drop")
  if (nrow(amb) > 0) ambiguous_ids[[source]] <<- amb %>% mutate(step = source)
  d %>% group_by(gene) %>% filter(n_distinct(ENTREZID) == 1) %>%
    summarise(ENTREZID = dplyr::first(ENTREZID), .groups = "drop") %>%
    mutate(entrez_source = source)
}

id_symbol <- safe_select(org.Bt.eg.db, annot_genes, "ENTREZID", "SYMBOL") %>%
  filter(!is.na(ENTREZID)) %>%
  transmute(gene = SYMBOL, ENTREZID) %>%
  unique_ids("org.Bt.eg.db symbol")

left_genes <- setdiff(annot_genes, id_symbol$gene)
id_alias <- safe_select(org.Bt.eg.db, left_genes, "ENTREZID", "ALIAS") %>%
  filter(!is.na(ENTREZID)) %>%
  group_by(gene = ALIAS) %>% filter(n_distinct(ENTREZID) == 1) %>%
  summarise(ENTREZID = dplyr::first(ENTREZID), .groups = "drop") %>%
  mutate(entrez_source = "org.Bt.eg.db alias")

left_genes <- setdiff(left_genes, id_alias$gene)
loc_genes  <- left_genes[str_detect(left_genes, "^LOC[0-9]+$")]
id_loc <- safe_select(org.Bt.eg.db, str_remove(loc_genes, "^LOC"), "SYMBOL", "ENTREZID") %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(ENTREZID) %>%
  transmute(gene = paste0("LOC", ENTREZID), ENTREZID, entrez_source = "LOC number (NCBI Gene ID)")

left_genes <- setdiff(left_genes, id_loc$gene)
id_bm <- bm_entrez %>%
  filter(external_gene_name %in% left_genes, !is.na(entrezgene_id), entrezgene_id != "") %>%
  transmute(gene = external_gene_name, ENTREZID = as.character(entrezgene_id)) %>%
  unique_ids("BioMart")

ids <- bind_rows(id_symbol, id_alias, id_loc, id_bm)

# NCBI description and current NCBI symbol for every Entrez ID
bt_info <- safe_select(org.Bt.eg.db, ids$ENTREZID, c("SYMBOL", "GENENAME"), "ENTREZID") %>%
  group_by(ENTREZID) %>%
  summarise(ncbi_symbol = dplyr::first(na.omit(SYMBOL)), description = dplyr::first(na.omit(GENENAME)), .groups = "drop")

desc <- tibble(gene = annot_genes) %>%
  left_join(ids, by = "gene") %>%
  left_join(bt_info, by = "ENTREZID")

orth <- bm_basic %>%
  filter(external_gene_name %in% annot_genes, !is.na(hsapiens_homolog_associated_gene_name),
         hsapiens_homolog_associated_gene_name != "") %>%
  group_by(gene = external_gene_name) %>%
  summarise(human_orthologues = paste(sort(unique(hsapiens_homolog_associated_gene_name)), collapse = ";"),
            n_orthologues = n_distinct(hsapiens_homolog_associated_gene_name), .groups = "drop")

# GO BP: bovine (org.Bt.eg.db) PLUS the GO BP of a human gene (org.Hs.eg.db) for every gene with
#   (i)  exactly one human orthologue in BioMart, or, for genes without a BioMart orthologue,
#   (ii) a human gene with the identical symbol (current NCBI symbol, else the RefSeq symbol; not LOC;
#        only symbols of a single human gene)
go_bt <- safe_select(org.Bt.eg.db, desc$ENTREZID, c("GO", "ONTOLOGY"), "ENTREZID") %>%
  filter(ONTOLOGY == "BP") %>% distinct(ENTREZID, GO) %>%
  inner_join(desc %>% filter(!is.na(ENTREZID)) %>% dplyr::select(gene, ENTREZID),
             by = "ENTREZID", relationship = "many-to-many") %>%
  transmute(gene, GO, go_source = "bovine GO")

hs_orth <- orth %>% filter(n_orthologues == 1) %>%
  transmute(gene, hs_symbol = human_orthologues, hs_match = "BioMart orthologue")
hs_same <- desc %>% filter(!gene %in% orth$gene) %>%
  transmute(gene, hs_symbol = coalesce(ncbi_symbol, gene)) %>%
  filter(!str_detect(hs_symbol, "^LOC"))
hs_same_ids <- safe_select(org.Hs.eg.db, hs_same$hs_symbol, "ENTREZID", "SYMBOL") %>%
  filter(!is.na(ENTREZID)) %>% group_by(SYMBOL) %>% filter(n_distinct(ENTREZID) == 1) %>% ungroup()
hs_same <- hs_same %>% filter(hs_symbol %in% hs_same_ids$SYMBOL) %>% mutate(hs_match = "identical symbol")
hs_genes <- bind_rows(hs_orth, hs_same)

hs_ids <- safe_select(org.Hs.eg.db, hs_genes$hs_symbol, "ENTREZID", "SYMBOL") %>%
  filter(!is.na(ENTREZID)) %>% distinct(SYMBOL, ENTREZID)
go_hs <- safe_select(org.Hs.eg.db, hs_ids$ENTREZID, c("GO", "ONTOLOGY"), "ENTREZID") %>%
  filter(ONTOLOGY == "BP") %>% distinct(ENTREZID, GO) %>%
  inner_join(hs_ids, by = "ENTREZID", relationship = "many-to-many") %>%
  inner_join(hs_genes, by = c("SYMBOL" = "hs_symbol"), relationship = "many-to-many") %>%
  transmute(gene, GO, go_source = paste0("human ", SYMBOL, " GO (", hs_match, ")"))

# bovine rows first, so a term present in both is attributed to the bovine annotation
go_by_gene <- bind_rows(go_bt, go_hs) %>% distinct(gene, GO, .keep_all = TRUE)
go_source_gene <- go_by_gene %>% group_by(gene) %>%
  summarise(go_source = paste(unique(go_source), collapse = " + "), .groups = "drop")
has_go <- unique(go_by_gene$gene)

ann <- desc %>%
  left_join(orth, by = "gene") %>%
  left_join(go_source_gene, by = "gene") %>%
  left_join(hs_genes %>% dplyr::select(gene, hs_gene_for_go = hs_symbol, hs_match), by = "gene") %>%
  mutate(n_orthologues = coalesce(n_orthologues, 0L),
         gene_label = if_else(str_detect(gene, "^LOC") & n_orthologues == 1,
                              paste0(gene, " (", human_orthologues, ")"), gene),
         desc_lc = str_to_lower(coalesce(description, "")))

cat("Genes:", nrow(ann), "| LOC:", sum(str_detect(ann$gene, "^LOC")),
    "| NCBI description:", sum(!is.na(ann$description)),
    "| Entrez ID:", sum(!is.na(ann$ENTREZID)),
    "| GO BP:", length(has_go), "(bovine:", n_distinct(go_bt$gene), "| human:", n_distinct(go_hs$gene),
    "| human only:", length(setdiff(unique(go_hs$gene), go_bt$gene)), ")",
    "| human orthologue >= 1:", sum(ann$n_orthologues > 0),
    "(exactly one:", sum(ann$n_orthologues == 1), ")",
    "| no BioMart data this run:", length(genes_without_biomart), "\n")
cat("\nHuman gene used for GO BP:\n")
ann %>% count(hs_match) %>% as.data.frame() %>% print(row.names = FALSE)
cat("\nEntrez ID source:\n")
ann %>% count(entrez_source) %>% as.data.frame() %>% print(row.names = FALSE)
cat("\nAmbiguous Entrez matches not used (symbol -> more than one Entrez ID; next step tried):\n")
if (length(ambiguous_ids) == 0) cat("none\n") else {
  amb_tbl <- bind_rows(ambiguous_ids) %>%
    left_join(ann %>% dplyr::select(gene, final_ENTREZID = ENTREZID, final_source = entrez_source), by = "gene")
  print(as.data.frame(amb_tbl), row.names = FALSE)
}
cat("\nBioMart: genes queried with no BioMart record (no orthologue or Entrez ID in Ensembl 115):",
    sum(annot_genes %in% bm$genes & !annot_genes %in% c(bm_basic$external_gene_name, bm_entrez$external_gene_name)),
    "of", sum(annot_genes %in% bm$genes), "queried\n")
cat("\nGenes matched through a previous/alternative symbol (RefSeq symbol -> current NCBI symbol):\n")
ann %>% filter(entrez_source == "org.Bt.eg.db alias") %>% transmute(gene, ncbi_symbol, description) %>%
  as.data.frame() %>% print(row.names = FALSE)


############################################################
# 5 — FUNCTIONAL CATEGORIES
############################################################

section("GO BP TERMS DEFINING EACH CATEGORY (check names)")
go_def <- tibble(category = rep(names(category_go), lengths(category_go)), GOID = unlist(category_go))
go_def <- go_def %>%
  left_join(suppressMessages(AnnotationDbi::select(GO.db, keys = go_def$GOID, columns = c("TERM", "ONTOLOGY"),
                                                   keytype = "GOID")), by = "GOID")
print(as.data.frame(go_def), row.names = FALSE)
bad_go <- go_def$GOID[is.na(go_def$TERM) | go_def$ONTOLOGY != "BP"]
if (length(bad_go) > 0) fail("GO term not found or not BP: ", paste(bad_go, collapse = ", "))

bp_offspring <- as.list(GO.db::GOBPOFFSPRING)
category_go_all <- lapply(category_go, function(ids) unique(c(ids, unlist(bp_offspring[ids]))))

go_hits <- bind_rows(lapply(names(category_go_all), function(cat_name)
  go_by_gene %>% filter(GO %in% category_go_all[[cat_name]]) %>% mutate(category = cat_name))) %>%
  mutate(priority = match(category, names(category_go)))
if (nrow(go_hits) > 0) {
  go_hits <- go_hits %>%
    left_join(suppressMessages(AnnotationDbi::select(GO.db, keys = unique(go_hits$GO), columns = "TERM",
                                                     keytype = "GOID")), by = c("GO" = "GOID"))
}
go_best <- go_hits %>%
  group_by(gene) %>%
  arrange(priority, GO, .by_group = TRUE) %>%
  summarise(go_category = dplyr::first(category),
            go_evidence = paste0(dplyr::first(GO), " ", dplyr::first(TERM),
                                 if_else(dplyr::first(go_source) == "bovine GO", "", paste0(" [", dplyr::first(go_source), "]"))),
            go_all_categories = paste(unique(category[order(priority)]), collapse = "; "),
            .groups = "drop")

fam <- ann %>%
  transmute(gene,
            rule_idx = vapply(desc_lc, function(d) {
              hit <- which(str_detect(d, family_rules$pattern)); if (length(hit)) hit[1] else NA_integer_
            }, integer(1))) %>%
  mutate(family_rule = family_rules$rule[rule_idx], family_category = family_rules$category[rule_idx]) %>%
  dplyr::select(gene, family_rule, family_category)

ann <- ann %>%
  left_join(fam, by = "gene") %>%
  left_join(go_best, by = "gene") %>%
  mutate(
    manual_category = unname(manual_categories[gene]),
    category = case_when(!is.na(manual_category) ~ manual_category,
                         !is.na(family_category) ~ family_category,
                         !is.na(go_category)     ~ go_category,
                         TRUE                    ~ "Other/Uncharacterised"),
    route = case_when(!is.na(manual_category) ~ "1 manual",
                      !is.na(family_category) ~ "2 gene family",
                      !is.na(go_category)     ~ "3 GO BP",
                      gene %in% has_go        ~ "4 GO BP, no category term",
                      TRUE                    ~ "5 unannotated"),
    manual_reason = unname(manual_reasons[gene]),
    evidence = case_when(route == "1 manual"      ~ if_else(is.na(manual_reason), "manual curation",
                                                            paste0("manual curation: ", manual_reason)),
                         route == "2 gene family" ~ paste0("NCBI description \"", description, "\" (rule: ", family_rule, ")"),
                         route == "3 GO BP"       ~ go_evidence,
                         route == "4 GO BP, no category term" ~ paste0(go_source, ": no term within the category terms"),
                         TRUE                     ~ NA_character_),
    # manual/family category not among the categories supported by the gene's GO terms
    conflict = route %in% c("1 manual", "2 gene family") & !is.na(go_all_categories) &
      !str_detect(coalesce(go_all_categories, ""), fixed(category))
  )
if (anyDuplicated(ann$gene)) fail("Duplicated genes in the annotation table.")

section("GENES PER ROUTE AND CATEGORY")
ann %>% count(route) %>% as.data.frame() %>% print(row.names = FALSE)
id_sym <- ann %>% filter(route == "3 GO BP", str_detect(coalesce(evidence, ""), "identical symbol"))
cat("\nGenes whose category rests on GO of a human gene with the identical symbol (no BioMart orthologue):",
    nrow(id_sym), "\n")
if (nrow(id_sym) > 0) id_sym %>% transmute(gene, category, evidence) %>% as.data.frame() %>% print(row.names = FALSE)
cat("\n")
ann %>% count(category, route) %>% pivot_wider(names_from = route, values_from = n, values_fill = 0) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nManual assignments used (", sum(ann$route == "1 manual"), "genes; * = GO supports other categories only):\n")
ann %>% filter(route == "1 manual") %>%
  transmute(gene = paste0(gene, if_else(conflict, " *", "")), category, description,
            go_supported_categories = coalesce(go_all_categories, "-")) %>%
  arrange(category, gene) %>% as.data.frame() %>% print(row.names = FALSE)
cat("\nManual entries not among the annotated genes (not used):",
    paste(setdiff(names(manual_categories), ann$gene), collapse = ", "), "\n")

cat("\nGene-family assignments where GO supports other categories only (for review):\n")
ann %>% filter(route == "2 gene family", conflict) %>%
  transmute(gene, category, family_rule, description, go_supported_categories = go_all_categories) %>%
  as.data.frame() %>% print(row.names = FALSE)

if (length(check_genes) > 0) {
  section("CHECK GENES")
  tibble(gene = check_genes) %>%
    left_join(ann %>% transmute(gene, category, route, evidence, go_supported_categories = go_all_categories,
                                hs_gene_for_go, hs_match), by = "gene") %>%
    mutate(category = coalesce(category, "(not among the annotated genes)")) %>%
    as.data.frame() %>% print(row.names = FALSE)
}


############################################################
# 6 — OUTPUT FILES
############################################################

evidence_tbl <- ann %>%
  transmute(gene, gene_label, ncbi_symbol, description, ENTREZID, entrez_source, human_orthologues, n_orthologues,
            hs_gene_for_go, hs_match, go_source, category, route, evidence, go_supported_categories = go_all_categories,
            family_rule, manual_category, conflict, is_DE = gene %in% de_genes,
            is_high_FC = gene %in% high_genes) %>%
  left_join(net_membership, by = "gene") %>%
  arrange(route, category, gene)
write_csv(evidence_tbl, file.path(out_dir, "gene_annotation_evidence.csv"))

# one row per annotated gene; DE columns from its lowest-FDR DE test (NA if not DE)
de_lowest <- DE_sig %>%
  mutate(gene = str_remove(gene, "^gene-")) %>%
  group_by(gene) %>%
  slice_min(FDR, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(gene, logFC, FDR, taxon_lowest_FDR = taxon, comparison_lowest_FDR = comparison)
gene_annot_for_export <- ann %>%
  dplyr::select(gene, gene_label, category, ENTREZID) %>%
  left_join(de_lowest, by = "gene") %>%
  left_join(net_membership, by = "gene") %>%
  mutate(is_DE = gene %in% de_genes,
         is_high_FC = gene %in% high_genes) %>%
  dplyr::select(gene, gene_label, logFC, FDR, category, ENTREZID, taxon_lowest_FDR, comparison_lowest_FDR,
                is_DE, is_high_FC, networks) %>%
  arrange(FDR)
write_csv(gene_annot_for_export, file.path(out_dir, "gene_annotation_table.csv"))

a$DE_all %>%
  mutate(gene_clean = str_remove(gene, "^gene-")) %>%
  left_join(ann %>% dplyr::select(gene_clean = gene, ENTREZID), by = "gene_clean") %>%
  dplyr::select(-gene_clean, -taxon) %>%
  relocate(ENTREZID, .after = gene) %>%
  write_csv(file.path(out_dir, "host_DE_all_results_entrezid.csv"))

cat("\nGenes per functional category (DE genes | high-FC DE genes | network genes | all annotated):\n")
gene_annot_for_export %>%
  group_by(category) %>%
  summarise(DE = sum(is_DE), high_FC = sum(is_high_FC), network = sum(!is.na(networks)), all = n(), .groups = "drop") %>%
  arrange(desc(all)) %>% as.data.frame() %>% print(row.names = FALSE)

saveRDS(list(
  params = list(annot_resp_top_n = annot_resp_top_n, network_edge_files = network_edge_files,
                manual_categories = manual_categories,
                family_rules = family_rules, category_go = category_go),
  gene_annot = ann %>% dplyr::select(gene_clean = gene, gene_label, description, ENTREZID,
                                     human_orthologues, category, route, evidence),
  categories = c(names(category_go), "Other/Uncharacterised"),
  created = Sys.time()
), file.path(out_dir, "05a_annotation.rds"))
cat("\nSaved: gene_annotation_evidence.csv, gene_annotation_table.csv,",
    "host_DE_all_results_entrezid.csv, 05a_annotation.rds\n")


############################################################
# 7 — VERSIONS AND SESSION INFORMATION
############################################################

section("ANNOTATION SOURCES")
cat("Ensembl:", bm$release, "|", bm$host, "| via", bm$method, "| last queried", format(bm$date), "\n")
if (requireNamespace("biomaRt", quietly = TRUE)) cat("biomaRt:", as.character(packageVersion("biomaRt")), "\n")
cat("\norg.Bt.eg.db:\n")
print(AnnotationDbi::metadata(org.Bt.eg.db) %>% filter(str_detect(name, "SOURCEDATE|SOURCENAME|DBSCHEMA")))
cat("\norg.Hs.eg.db:\n")
print(AnnotationDbi::metadata(org.Hs.eg.db) %>% filter(str_detect(name, "SOURCEDATE|SOURCENAME|DBSCHEMA")))
cat("\nGO.db:\n")
print(AnnotationDbi::metadata(GO.db) %>% filter(str_detect(name, "SOURCEDATE|SOURCENAME")))

section("SESSION INFORMATION")
print(sessionInfo())
cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

sink(type = "message")
sink()
close(log_con)
