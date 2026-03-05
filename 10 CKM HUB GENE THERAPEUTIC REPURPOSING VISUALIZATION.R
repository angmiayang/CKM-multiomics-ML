# ==============================================================================
# SCRIPT: CKM HUB GENE THERAPEUTIC REPURPOSING VISUALIZATION
# ==============================================================================

# 1. SETUP LIBRARIES
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("patchwork")) install.packages("patchwork")

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# 2. LOAD & PREPARE DATA
# ------------------------------------------------------------------------------
# Ensure the file is in your working directory
df <- read.csv("CKM_Top_10_Druggable_Genes.csv", stringsAsFactors = FALSE)

# Sort by Druggability Score and set Gene as factor for correct plotting order
df <- df %>%
  arrange(desc(Druggability_Score)) %>%
  mutate(Gene = factor(Gene, levels = rev(Gene)))

# Prepare long-format data for the drug comparison plot
# Order factor levels so "Total" appears before "Approved"
df_long <- df %>%
  select(Gene, Total_Drugs, Approved_Drugs) %>%
  pivot_longer(cols = c("Total_Drugs", "Approved_Drugs"), 
               names_to = "Drug_Status", 
               values_to = "Count") %>%
  mutate(Drug_Status = factor(Drug_Status, levels = c("Total_Drugs", "Approved_Drugs")))

# 3. DEFINE NATURE-INSPIRED THEME & COLORS
# ------------------------------------------------------------------------------
# Hex codes based on Nature Publishing Group (NPG) palettes
nature_dark <- "#3C5488" # Dark Blue
nature_sky  <- "#4DBBD5" # Sky Blue
nature_red  <- "#E64B35" # Coral Red

nature_theme <- theme_bw() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    text = element_text(family = "sans", size = 12),
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "italic", color = "black"), # Italics for Genes
    strip.background = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5)
  )

# 4. CREATE PLOT A: DRUGGABILITY SCORE
# ------------------------------------------------------------------------------
plot_a <- ggplot(df, aes(x = Druggability_Score, y = Gene)) +
  geom_col(fill = nature_dark, width = 0.7) +
  geom_text(aes(label = round(Druggability_Score, 1)), 
            hjust = 1.2, color = "white", size = 3, fontface = "bold") +
  labs(title = "Druggability Score",
       x = "Score",
       y = "Gene Symbol") +
  nature_theme

# 5. CREATE PLOT B: DRUG LANDSCAPE (TOTAL VS APPROVED)
# ------------------------------------------------------------------------------
plot_b <- ggplot(df_long, aes(x = Count, y = Gene, fill = Drug_Status)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("Total_Drugs" = nature_sky, "Approved_Drugs" = nature_red),
                    labels = c("Total Drugs", "Approved Drugs")) +
  labs(title = "Drug Repurposing Landscape",
       x = "Number of Drugs",
       y = NULL,
       fill = NULL) +
  nature_theme +
  theme(axis.text.y = element_blank(), # Remove redundant gene names
        axis.ticks.y = element_blank())

# 6. COMBINE AND SAVE
# ------------------------------------------------------------------------------
message("Generating CKM Hub Gene Plots...")

# Combine using patchwork
final_plot <- (plot_a | plot_b) + 
  plot_layout(widths = c(1, 1.2))

# Save high-resolution outputs
ggsave("CKM_Hub_Genes_NatureStyle.png", plot = final_plot, width = 12, height = 7, dpi = 300)
ggsave("CKM_Hub_Genes_NatureStyle.pdf", plot = final_plot, width = 12, height = 7)

# Print to viewer
print(final_plot)

message("✅ Visualizations saved: CKM_Hub_Genes_NatureStyle.png and .pdf")