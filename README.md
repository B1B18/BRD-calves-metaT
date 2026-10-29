# Nasal metatranscriptomics of pre-weaned dairy calves with and without respiratory disease

This repository contains the R scripts used to analyse the data in the manuscript "An integrated metatranscriptomic study of viruses, bacteria and host transcriptional signatures in nasal swabs from pre-weaned dairy calves with and without respiratory disease" by Brito et al. 

The study profiles viruses, bacteria and the host transcriptome in nasal swab metatranscriptomes from pre-weaned dairy calves sampled on dairy farms in Australia. The scripts reproduce the statistical analyses, tables and figures reported in the manuscript. The usage of each script is described below and in the header of each script.

## Data

The raw sequencing reads (Illumina NovaSeq X, 150 bp paired-end) are deposited in the NCBI Sequence Read Archive under BioProject [PRJNA1191164](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1191164).

Reads were processed upstream of these scripts, as described in the Methods of the manuscript:

- Host read removal and rRNA depletion (Hostile, SortMeRNA)
- Viral abundance: de novo assembly and read mapping to a custom bovine viral reference database ([doi:10.5281/zenodo.21867952](https://doi.org/10.5281/zenodo.21867952)), summarised as normalised abundance per viral species
- Bacterial taxonomy: MetaPhlAn 4 with the GTDB taxonomy (r220), summarised as read counts per genus
- Host gene expression: nf-core/rnaseq (STAR–RSEM) against the bovine reference genome ARS-UCD1.2, summarised as gene-level counts

The processed tables produced by these steps, which are the inputs of the scripts, are in the `input_files/` folder. Its README describes each file and its columns:

- `metadata.csv`: sample metadata (clinical status, age, anonymised farm code, sequencing depth)
- `taxa_long.csv`: viral normalised abundance and bacterial read counts per sample and taxon
- `host_counts.csv.gz`: host gene-level read counts
- `all_merged_meta_gtdb_genus.txt`: MetaPhlAn 4 read counts for all bacterial genera (used by `01_permanova.R` only)
- `input_annotations.rds`: cache of the Ensembl BioMart query (release 115) used by `05a_gene_annotation.R`, so the published gene annotation can be reproduced without querying Ensembl

## Scripts

The scripts are in `scripts/` and are run in numerical order from the repository root. Each script writes its results, a log file (`log.txt`, with the numbers reported in the Results and the R session information) and its figures to a folder with the same name as the script.

1. `01_permanova.R` performs PERMANOVA (vegan `adonis2`, marginal terms, 9,999 permutations) of the viral, bacterial and host expression profiles against clinical status, age group (<25 / >=25 days) and farm, on one common set of calves. Viral and bacterial profiles use Euclidean distance on log1p-transformed abundance; host expression uses Euclidean distance on TMM-normalised log CPM. The script also runs a sensitivity analysis with permutations restricted within farms and tests homogeneity of multivariate dispersion (`betadisper`). Output: Table 2.
2. `02_hurdle_models.R` fits two-part (hurdle) mixed models per viral species and bacterial genus, with farm as a random intercept: a binomial model for detection and a linear model for log abundance in samples where the taxon was detected. Respiratory and unspecific calves are each compared with healthy calves, separately within each age group, and p-values are BH-adjusted. Output: Figure 3 and the model tables.
3. `03a_host_DE.R` performs host differential expression with limma-voom, with farm as a blocking factor (`duplicateCorrelation`) and adjusted for age. Two contrasts are fitted: respiratory vs healthy calves, and, for each taxon, taxon-positive vs taxon-negative calves adjusted for clinical status. Taxon detection is defined by the lower quartile of non-zero abundance (Q25 rule). Genes are filtered and TMM-normalised separately for each model. Output: DE result tables and QC plots; `03a_objects.rds` is passed to scripts 05a and 05b.
4. `04a_networks_DE_genes.R` builds host–taxon Spearman correlation networks separately for healthy and respiratory calves. Edges are kept if |rho| >= 0.4, BH-adjusted p <= 0.05 and the edge is stable in >= 80% of 1,000 bootstrap resamples of the calves. The script also summarises the overlap between correlated genes and the DE results from 03a. Output: network edge and node tables.
5. `05a_gene_annotation.R` annotates every host gene reported in the results (DE genes and network genes) with Entrez Gene IDs, NCBI descriptions, human orthologues (Ensembl BioMart) and GO biological process terms. Each gene is assigned one functional category by a fixed, documented order of rules (manual curation, gene-family rules, GO terms). No enrichment testing is performed. Output: the gene annotation table (supplementary table).
6. `05b_figures.R` draws Figures 6 and 7 (DE gene overlaps between taxa, heatmaps, category summaries and the correlation networks) and writes the supplementary table of all DE results. It reads the outputs of 03a, 04a and 05a and re-fits nothing.

Run order: `01` and `02` are independent; `03a` → `04a` → `05a` → `05b` must be run in this order.

For example, from the repository root:

```sh
Rscript scripts/03a_host_DE.R
```

## Software

The analyses were run in R 4.5.2. Packages (CRAN and Bioconductor):

- Data handling and plotting: tidyverse, tidytext, ggraph, igraph, pheatmap, ComplexHeatmap, circlize, eulerr, gridExtra
- Statistics: vegan, permute, lme4, lmerTest, broom.mixed, psych, foreach, doParallel, doRNG
- Expression analysis: edgeR, limma
- Annotation: AnnotationDbi, org.Bt.eg.db, org.Hs.eg.db, GO.db, biomaRt

The exact package versions used are listed in the session information at the end of each script's `log.txt`.

## Citation

If you use these scripts, please cite the manuscript:

Brito B, et al. An integrated metatranscriptomic study of viruses, bacteria and host transcriptional signatures in nasal swabs from pre-weaned dairy calves with and without respiratory disease. [JOURNAL, YEAR]. [DOI]

and this repository: [ZENODO DOI]

## License

The code in `scripts/` is released under the MIT License (see `LICENSE`). The data in `input_files/` are released under the Creative Commons Attribution 4.0 International licence (CC BY 4.0).
