# ==============================================================================
# PHASE II: MULTI-OMIC MATRIX CONSTRUCTION
# ==============================================================================
# Project Root: ~/CKM_GWAS/
# Input:        CKM_Munged_Results/ (6 High-Power Traits)
#               GTEx_v10/CKM_GTEx_V10_Features.csv
# Output:       Master_CKM_ML_Matrix.csv
# ==============================================================================

library(data.table)
library(biomaRt)

message("--- STARTING PHASE II: MATRIX CONSTRUCTION ---")

# 1. LOAD GTEx FEATURES (The CSV you just verified)
gtex_features <- fread("GTEx_v10/CKM_GTEx_V10_Features.csv")
message("Loaded ", nrow(gtex_features), " genes from GTEx V10.")

# 2. RETRIEVE GENE COORDINATES (GRCh37)
# Connect to Ensembl to get the genomic locations of all genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://grch37.ensembl.org")
genes_ref <- as.data.table(getBM(
  attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"),
  filters = "chromosome_name", values = c(1:22, "X", "Y"), mart = ensembl
))
setnames(genes_ref, c("Gene", "CHR", "Start", "End"))

# Define 50kb Window (to capture nearby regulatory SNPs)
genes_ref[, `:=`(W_Start = Start - 50000, W_End = End + 50000)]
genes_ref[, CHR := as.character(CHR)]

# 3. MAPPING FUNCTION (SNP to Gene)
map_trait <- function(trait_name, file_path) {
  message("Mapping Pillar: ", trait_name)
  snps <- fread(file_path, select = c("CHR", "BP", "P", "BETA"))
  snps[, CHR := as.character(CHR)]
  
  # Non-equi join to find SNPs inside each gene's 50kb window
  mapped <- snps[genes_ref, on = .(CHR, BP >= W_Start, BP <= W_End), nomatch = 0]
  
  # Pick the best P-value for each gene
  res <- mapped[, .(P = min(P, na.rm = TRUE), BETA = BETA[which.min(P)]), by = Gene]
  setnames(res, c("P", "BETA"), paste0(c("P_", "BETA_"), trait_name))
  return(res)
}

# 4. PROCESS THE 6 MUNGED PILLARS
munged_dir <- "CKM_Munged_Results"

heart_all_df   <- map_trait("heart",    file.path(munged_dir, "munged_heart_hf_all.tsv.gz"))
egfr_df        <- map_trait("egfr",     file.path(munged_dir, "munged_kidney_egfr.tsv.gz"))
egfrcys_df     <- map_trait("egfrcys",  file.path(munged_dir, "munged_kidney_egfrcys.tsv.gz"))
bun_df         <- map_trait("bun",      file.path(munged_dir, "munged_kidney_bun.tsv.gz"))
bmi_df         <- map_trait("bmi",      file.path(munged_dir, "munged_metab_bmi.tsv.gz"))
whr_df         <- map_trait("whr",      file.path(munged_dir, "munged_metab_whr.tsv.gz"))

# 5. ASSEMBLE MASTER MATRIX
message("Consolidating all features...")
master <- merge(gtex_features, heart_all_df, by = "Gene", all.x = TRUE)
master <- merge(master, egfr_df,      by = "Gene", all.x = TRUE)
master <- merge(master, egfrcys_df,   by = "Gene", all.x = TRUE)
master <- merge(master, bun_df,       by = "Gene", all.x = TRUE)
master <- merge(master, bmi_df,       by = "Gene", all.x = TRUE)
master <- merge(master, whr_df,       by = "Gene", all.x = TRUE)

# 6. RANDOM FOREST PREPARATION (Imputation)
# Machine Learning models cannot handle missing values (NAs)
# We fill missing P-values with 1.0 (no significance) and BETAs with 0
p_cols <- grep("^P_", names(master), value = TRUE)
for(col in p_cols) master[is.na(get(col)), (col) := 1.0]
master[is.na(master)] <- 0

# Save the final matrix
fwrite(master, "Master_CKM_ML_Matrix.csv")
message("--- SUCCESS: Master Matrix Constructed ---")