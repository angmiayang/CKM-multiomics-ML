This repository contains a 10-step multi-omic machine learning pipeline designed to prioritize genes and therapeutic targets for Cardiovascular-Kidney-Metabolic (CKM) Syndrome. The project integrates high-power GWAS summary statistics with tissue-specific transcriptomic features from GTEx V10 to identify key drivers of CKM risk.

🧬 Pipeline Architecture

The workflow is divided into three functional phases: Data Harmonization, Machine Learning Prioritization, and Functional Validation.

Phase I: Data Processing & Harmonization

Genomic Harmonization (01): Standardizes GWAS summary statistics for six high-power pillars (Heart Failure, eGFR, eGFRcys, BUN, BMI, and WHR) using the MungeSumstats framework and GRCh37 reference genome.

GTEx V10 Feature Extraction (02): Extracts median TPM expression values for target tissues representing the three CKM pillars: Heart (Left Ventricle/Atrial Appendage), Kidney (Cortex), and Metabolic (Subcutaneous Adipose).

Matrix Construction (03): Maps SNPs to genes within a 50kb window and assembles a master multi-omic matrix, handling missing values through targeted imputation (P=1.0 for non-significance).

Phase II: Machine Learning & Filtration

Random Forest Prioritization (04): Trains a Random Forest model with 1,000 trees to calculate Gini Importance scores and a consolidated "CKM_Score" for all genes.

Protein-Coding Filtration (05): Filters results to retain only protein-coding genes, removing non-coding RNAs and pseudogenes to prepare a top 1,000 candidate list for STRING-db analysis.

Phase III: Functional & Therapeutic Insight

Batch Enrichr Filtering (06): A Python utility to filter pathway enrichment results (GO, KEGG, MSigDB) based on a Gene Count > 3 and Adjusted P-value < 0.05.

GO Enrichment Visualization (07): Generates faceted dot plots of top Gene Ontology terms across Biological Process, Molecular Function, and Cellular Component categories.

KEGG Pathway Network (08): Constructs a Jaccard-based gene overlap network for the top 20 KEGG pathways to visualize functional clusters.

Druggability Scoring (09): Utilizes DGIdb data to calculate a weighted clinical-evidence score for 58 hub genes, identifying the most viable targets for repurposing.

Therapeutic Visualization (10): Produces "Nature-style" visualizations of the drug repurposing landscape for the top 10 druggable CKM hub genes.

🛠️ Installation & Requirements

R Dependencies

The scripts require R (≥ 4.0) and the following libraries:

data.table, MungeSumstats, BiocParallel, biomaRt, randomForest, ggplot2, igraph, dplyr, patchwork

Python Dependencies

pandas, numpy, os
