# ==============================================================================
# PHASE III: RANDOM FOREST PRIORITIZATION
# ==============================================================================
# Project Root: ~/CKM_GWAS/
# Input:        Master_CKM_ML_Matrix.csv
# Output:       CKM_Prioritized_Genes_Full.csv
#               Feature_Importance_Plot.png
# ==============================================================================

library(data.table)
library(randomForest)
library(ggplot2)

message("--- STARTING PHASE III: RANDOM FOREST TRAINING ---")

# 1. LOAD THE MASTER MATRIX
master <- fread("Master_CKM_ML_Matrix.csv")

# 2. PREPARE DATA FOR ML
# We use the Heart Failure P-value as the 'Target' (Label) to see which other 
# features best predict cardiac risk in the context of CKM.
# For ranking, we look at the 'Importance' of all features.

# Remove the Gene name column for training, but keep it for the final join
rf_data <- master[, !c("Gene"), with = FALSE]

# 3. TRAIN THE RANDOM FOREST
# ntree = 1000 provides stable importance scores for 50k+ genes
message("Training Random Forest (1,000 trees)... This may take 20-40 mins.")
set.seed(42) # For reproducibility
rf_model <- randomForest(x = rf_data[, !c("P_heart"), with = FALSE], 
                         y = rf_data$P_heart, 
                         ntree = 1000, 
                         importance = TRUE)

# 4. EXTRACT FEATURE IMPORTANCE
importance_scores <- as.data.table(importance(rf_model), keep.rownames = "Feature")
setorder(importance_scores, -IncNodePurity)

# 5. GENERATE GENE PRIORITIZATION RANKINGS
# We calculate a 'CKM_Score' based on the weighted contribution of all pillars
# Lower P-values (more significant) get higher scores.
master[, CKM_Score := rowMeans(1 - master[, grep("^P_", names(master)), with = FALSE])]

# Sort by the final CKM score
final_ranking <- master[order(-CKM_Score)]

# 6. SAVE RESULTS
fwrite(final_ranking, "CKM_Prioritized_Genes_Full.csv")
fwrite(importance_scores, "Feature_Importance_Scores.csv")

# 7. VISUALIZE TOP FEATURES
p <- ggplot(importance_scores[1:10], aes(x = reorder(Feature, IncNodePurity), y = IncNodePurity)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 10 CKM Predictors", x = "Genomic/Tissue Feature", y = "Gini Importance") +
  theme_minimal()

ggsave("Feature_Importance_Plot.png", p, width = 8, height = 6)

message("--- SUCCESS ---")
message("Top Priority Gene: ", final_ranking$Gene[1])
message("Results saved to CKM_Prioritized_Genes_Full.csv")