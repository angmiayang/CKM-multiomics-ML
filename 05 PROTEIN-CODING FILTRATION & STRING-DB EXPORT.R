# ==============================================================================
# SCRIPT: PROTEIN-CODING FILTRATION & STRING-DB EXPORT
# ==============================================================================

library(data.table)

# 1. LOAD ML RESULTS
results <- fread("CKM_Prioritized_Genes_Full.csv")

# 2. STRICT PROTEIN-CODING FILTRATION
# Removes non-coding RNAs (MIR, RNU, SNOR), pseudogenes (ending in P), 
# and antisense transcripts (-AS) which lack PPI data.
clean_results <- results[
    !grepl("^MIR|^RNU|^Y_RNA|^SNOR|^LINC|^LOC|^U[0-9]|^RN7SL|P[0-9]+$|-AS[0-9]$|ENSG", Gene)
]

# 3. SAVE THE CLEANED RANKINGS
fwrite(clean_results, "CKM_Prioritized_Genes_CLEAN.csv")

# 4. EXTRACT TOP 1000 PROTEINS FOR FUNCTIONAL VALIDATION
# We use 1000 as the standard cutoff for building a robust PPI network
top_1000_genes <- clean_results$Gene[1:1000]

# 5. WRITE STRING-DB INPUT LIST
# This creates a simple text file for easy copy-pasting into the STRING-db web interface
write.table(top_1000_genes, "STRING_Input_List.txt", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)

message("--- SUCCESS ---")
message("Protein-coding filtration complete.")
message("Current Top Gene: ", clean_results$Gene[1])
message("Total protein-coding genes identified: ", nrow(clean_results))
message("Top 1000 list saved for STRING-db analysis.")