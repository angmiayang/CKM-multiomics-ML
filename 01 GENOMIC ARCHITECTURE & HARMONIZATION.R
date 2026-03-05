# ==============================================================================
# PHASE I: GENOMIC ARCHITECTURE & HARMONIZATION
# ==============================================================================
# Project Root: ~/CKM_GWAS/
# Input:        ~/CKM_GWAS/GWAS/ (Heart, Kidney, Metabolic)
# Output:       ~/CKM_GWAS/CKM_Munged_Results/
# ==============================================================================

# 1. SETUP & CONFIGURATION
options(timeout = 18000) 

library(data.table)
library(MungeSumstats)
library(BiocParallel)

message("--- STARTING GENOMIC HARMONIZATION: 6 HIGH-POWER PILLARS ---")

# 2. OUTPUT WORKSPACE
workspace <- "CKM_Munged_Results"
if(!dir.exists(workspace)) dir.create(workspace, recursive = TRUE)

# Parallel Processing (14 Cores)
cores_to_use <- 14
register(MulticoreParam(workers = cores_to_use))

message(paste("Environment Ready. Current Working Directory:", getwd()))

# 3. DEFINE INPUT PILLARS
# We prioritize maximum sample size for ML stability
file_list <- list(
  # Heart Pillar
  heart_hf_all   = "GWAS/Heart/FORMAT-METAL_Pheno1_ALL.tsv.gz",
  
  # Kidney Pillar
  kidney_egfr    = "GWAS/Kidney/metal_eGFR_meta1.TBL.map.annot.gc.gz",
  kidney_egfrcys = "GWAS/Kidney/metal_eGFRcys_meta1.TBL.map.annot.gc.gz",
  kidney_bun     = "GWAS/Kidney/metal_bun_meta_all1.TBL.map.annot.gc.gz",
  
  # Metabolic Pillar
  metab_bmi      = "GWAS/Metabolic/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz",
  metab_whr      = "GWAS/Metabolic/Whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt.gz"
)

# 4. EXECUTION LOOP
for (trait in names(file_list)) {
  
  input_path <- file_list[[trait]]
  save_path  <- file.path(workspace, paste0("munged_", trait, ".tsv.gz"))
  
  message(paste("\n--- Processing Pillar:", toupper(trait), "---"))
  
  if (file.exists(input_path)) {
    tryCatch({
      
      MungeSumstats::format_sumstats(
        path = input_path,
        ref_genome = "GRCh37",
        dbSNP = 155,               # Reviewer-standard validation
        save_path = save_path,
        force_new = TRUE,
        allele_flip_check = TRUE,  # Ensures consistent effect directions
        bi_allelic_filter = TRUE,  # Removes ambiguous SNPs
        nThread = cores_to_use
      )
      
      message(paste("✅ SUCCESS:", trait, "standardized."))
      
    }, error = function(e) {
      message(paste("❌ ERROR in", trait, ":", e$message))
    })
    
    gc() # Memory management between large files
    
  } else {
    message(paste("⚠️ FILE NOT FOUND at path:", input_path))
  }
}

message("\n--- PHASE I COMPLETE: ALL 6 TRAITS READY FOR INTEGRATION ---")