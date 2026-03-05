# ==============================================================================
# SCRIPT: GENERATE GO ENRICHMENT DOT PLOT (SORTED BY P-VALUE)
# ==============================================================================
# Purpose: Create a faceted dot plot of Top 10 GO terms.
#          Selection AND Sorting are based on Adjusted P-value.
# Input:   Filtered_GO_*.csv files
# Output:  GO_Enrichment_Facet_Plot_Sorted.pdf
# ==============================================================================

# 1. SETUP LIBRARIES
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("data.table")) install.packages("data.table")
if (!require("dplyr")) install.packages("dplyr")
if (!require("stringr")) install.packages("stringr")

library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)

# 2. LOAD DATA
# ------------------------------------------------------------------------------
# Load the 3 filtered datasets
bp <- fread("Filtered_GO_Biological_Process_2025_table.csv")
bp$Category <- "(1) Biological Process"

mf <- fread("Filtered_GO_Molecular_Function_2025_table.csv")
mf$Category <- "(2) Molecular Function"

cc <- fread("Filtered_GO_Cellular_Component_2025_table.csv")
cc$Category <- "(3) Cellular Component"

# 3. SELECT TOP 10 TERMS (By Adjusted P-value)
# ------------------------------------------------------------------------------
process_data <- function(df, n = 10) {
  df %>%
    arrange(`Adjusted P-value`) %>%   # 1. Sort by significance (Smallest P first)
    head(n) %>%                       # 2. Grab top 10
    as.data.table()
}

combined_df <- rbind(process_data(bp), process_data(mf), process_data(cc))

# 4. PREPARE FOR PLOTTING
# ------------------------------------------------------------------------------
# Calculate -log10 P-value for coloring and sorting
combined_df[, LogP := -log10(`Adjusted P-value`)]

# Clean term names (Remove GO IDs)
combined_df[, Term_Clean := str_remove(Term, " \\(GO:.*\\)")]

# --- CRITICAL STEP: SORTING ---
# We set the Factor Levels based on LogP (ascending).
combined_df$Term_Clean <- factor(combined_df$Term_Clean, 
                                 levels = unique(combined_df$Term_Clean[order(combined_df$LogP)]))

# 5. GENERATE PLOT
# ------------------------------------------------------------------------------
message("Generating Plot...")

p <- ggplot(combined_df, aes(x = `Odds Ratio`, y = Term_Clean)) +
  
  # The Dots
  geom_point(aes(size = Gene_Count, color = LogP)) +
  
  # Facet Split
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  
  # Color Scale
  scale_color_gradient(low = "#A6D96A", high = "#D7191C", name = "-log10(Adj. P)") +
  
  # Labels
  labs(x = "Odds Ratio", 
       y = "Enriched GO Terms",   # Label for the y-axis
       size = "Gene Count",
       title = NULL) +           # No top title
  
  theme_bw() +
  theme(
    # --- Y-AXIS LABEL (PLAIN TEXT, CENTERED) ---
    axis.title.y = element_text(size = 12, face = "plain", color = "black", 
                                hjust = 0.5, vjust = 0.5),
    axis.text.y = element_text(size = 10, color = "black"),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank()
  )

# 6. SAVE
# ------------------------------------------------------------------------------
ggsave("GO_Enrichment_Facet_Plot_Sorted.pdf", plot = p, width = 8, height = 12)

message("✅ Plot saved to: GO_Enrichment_Facet_Plot_Sorted.pdf")