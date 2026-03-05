# ==============================================================================
# SCRIPT: KEGG PATHWAY NETWORK (LIGHTER GREEN THEME)
# ==============================================================================
# Purpose: Pairwise network of Top 20 KEGG pathways.
# Updates: 
#   1. Fixed "Vertex Name" error (Term_Clean is now Col 1).
#   2. Fixed "Label" error (Using V(g)$name).
#   3. CHANGED COLOR: Used a custom, lighter green palette.
# ==============================================================================

# 1. SETUP LIBRARIES
if (!require("igraph")) install.packages("igraph")
if (!require("data.table")) install.packages("data.table")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("stringr")) install.packages("stringr")
if (!require("scales")) install.packages("scales")

library(igraph)
library(data.table)
library(RColorBrewer)
library(stringr)
library(scales)

# 2. LOAD & PREPARE DATA
# ------------------------------------------------------------------------------
df <- fread("Filtered_KEGG_2026_table.csv")

# Sort by Adjusted P-value (Smallest first) and take Top 20
df <- df[order(`Adjusted P-value`)][1:min(20, nrow(df))]

# Clean Term Names: Title Case
df[, Term_Clean := str_to_title(Term)]

# --- CRITICAL FIX: Set Term_Clean as the First Column ---
setcolorder(df, c("Term_Clean", setdiff(names(df), "Term_Clean")))

# 3. CALCULATE GENE OVERLAP (JACCARD INDEX)
# ------------------------------------------------------------------------------
gene_list <- str_split(df$Genes, ";")
names(gene_list) <- df$Term_Clean

n <- nrow(df)
edge_list <- data.frame()

message("Calculating Pairwise Overlap...")

for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    genes_i <- gene_list[[i]]
    genes_j <- gene_list[[j]]
    
    intersection <- length(intersect(genes_i, genes_j))
    union_len    <- length(union(genes_i, genes_j))
    jaccard      <- intersection / union_len
    
    if (jaccard > 0.05) { 
      edge_list <- rbind(edge_list, data.frame(
        from = df$Term_Clean[i],
        to   = df$Term_Clean[j],
        weight = jaccard,
        overlap = intersection
      ))
    }
  }
}

# 4. BUILD THE GRAPH
# ------------------------------------------------------------------------------
g <- graph_from_data_frame(d = edge_list, vertices = df, directed = FALSE)

# 5. DEFINE VISUAL ATTRIBUTES
# ------------------------------------------------------------------------------

# A. NODE SIZE: Rescaled 10-25
V(g)$size <- rescale(V(g)$Gene_Count, to = c(10, 25))

# B. NODE COLOR (LIGHTER THEME)
# ------------------------------------------------------------------------------
logp <- -log10(V(g)$`Adjusted P-value`)

# Custom Palette: Starts Pale Green -> Ends Vivid Green (No Black-Green)
# Hex codes: #E5F5E0 (Light) -> #74C476 (Medium) -> #238B45 (Vibrant)
pal <- colorRampPalette(c("#E5F5E0", "#74C476", "#238B45"))(100)

V(g)$color <- pal[cut(logp, breaks = 100, labels = FALSE)]

# C. EDGE THICKNESS
E(g)$width <- E(g)$weight * 15  

# D. LABELS
V(g)$label <- V(g)$name 
V(g)$label.cex <- 0.8           
V(g)$label.color <- "black"
V(g)$label.family <- "sans"      
V(g)$label.font <- 2             

# 6. PLOT AND SAVE
# ------------------------------------------------------------------------------
message("Generating Network Plot...")

pdf("KEGG_Pathway_Network_Lighter.pdf", width = 14, height = 14)

set.seed(123) 
layout <- layout_with_fr(g)

plot(g, 
     layout = layout,
     vertex.frame.color = "#999999",  # Lighter Grey Border
     vertex.label.dist = 0,           
     edge.color = "#BEBEBE80",        
     main = "Top 20 Enriched KEGG Pathways (Gene Overlap Network)"
)

# Legend (Updated to match new colors)
legend("topright", 
       legend = c("Less Significant", "More Significant"), 
       fill = c(pal[1], pal[100]), 
       title = "Significance (-logP)",
       bty = "n",
       cex = 0.8)

dev.off()

message("✅ Network Diagram saved to: KEGG_Pathway_Network.pdf")