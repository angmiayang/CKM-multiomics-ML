# ==============================================================================
# SCRIPT: GTEx V10 FEATURE EXTRACTION
# ==============================================================================
# Project Root: ~/CKM_GWAS/
# Input:        ~/CKM_GWAS/GTEx_v10/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz
# Output:       ~/CKM_GWAS/GTEx_v10/CKM_GTEx_V10_Features.csv
# ==============================================================================

library(data.table)

# 1. PATH CONFIGURATION
# Pointing to your specific v10 path
gtex_input  <- "GTEx_v10/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"
gtex_output <- "GTEx_v10/CKM_GTEx_V10_Features.csv"

message("--- STARTING GTEx V10 EXTRACTION ---")

# 2. DEFINE TARGET TISSUES
# These columns correspond to the 3 pillars of CKM Syndrome
target_cols <- c(
    "Description",              # This contains the Gene Symbols (e.g., ACE2)
    "Heart_Left_Ventricle",     # Cardiovascular Pillar
    "Heart_Atrial_Appendage",   # Cardiovascular Pillar
    "Kidney_Cortex",            # Renal Pillar
    "Adipose_Subcutaneous"      # Metabolic Pillar
)

# 3. LOAD DATA
# We skip the first 2 metadata lines of the GCT file and use tab separation
message("Reading large GTEx v10 file (57,853 genes)...")
gtex_data <- fread(gtex_input, 
                  skip = 2, 
                  select = target_cols, 
                  sep = "\t", 
                  header = TRUE)

# 4. CLEANING & STANDARDIZATION
# Rename 'Description' to 'Gene' to match our GWAS mapping scripts
setnames(gtex_data, "Description", "Gene")

# Handle duplicate Gene symbols by taking the maximum expression value
# (Ensures 1 row per Gene for the Random Forest matrix)
gtex_final <- gtex_data[, lapply(.SD, max), by = Gene]

# 5. SAVE
fwrite(gtex_final, gtex_output)

message("--- SUCCESS ---")
message("Extracted features for: ", nrow(gtex_final), " genes.")
message("Saved to: ", gtex_output)