import pandas as pd
import numpy as np
import os

# 1. LOAD YOUR 58 HUB GENES
# ------------------------------------------------------------------------------
hub_file = '58_Hub_Genes.txt'
if os.path.exists(hub_file):
    with open(hub_file, 'r') as f:
        hub_genes = [line.strip() for line in f if line.strip()]
    print(f"✅ Loaded {len(hub_genes)} hub genes for filtering.")
else:
    print(f"❌ Error: {hub_file} not found.")
    exit()

# 2. LOAD DGIdb RAW DATA (TSV format)
# ------------------------------------------------------------------------------
print(">>> Loading interactions, categories, and drugs metadata...")
interactions = pd.read_csv('interactions.tsv', sep='\t', low_memory=False)
categories = pd.read_csv('categories.tsv', sep='\t', low_memory=False)

# 3. FILTERING FOR CKM HUB GENES
# ------------------------------------------------------------------------------
ckm_interactions = interactions[interactions['gene_name'].isin(hub_genes)].copy()
ckm_categories = categories[categories['name'].isin(hub_genes)].copy()

# 4. CALCULATE DRUGGABILITY SCORE
# ------------------------------------------------------------------------------
# We quantify 'druggability' using a weighted clinical-evidence formula:
# Score = (Approved_Drugs * 3) + (Total_Drugs * 1) + (Categories * 2) + (Avg_Evidence_Score * 5)
druggability_list = []

for gene in hub_genes:
    # Filter subsets
    gene_ints = ckm_interactions[ckm_interactions['gene_name'] == gene]
    gene_cats = ckm_categories[ckm_categories['name'] == gene]['name-2'].unique()
    
    # Calculate Metrics
    num_total_drugs = gene_ints['drug_name'].nunique()
    
    # Check approval status (handles boolean or string 'TRUE')
    if gene_ints.empty:
        num_approved_drugs = 0
    elif gene_ints['approved'].dtype == 'bool':
        num_approved_drugs = gene_ints[gene_ints['approved'] == True]['drug_name'].nunique()
    else:
        num_approved_drugs = gene_ints[gene_ints['approved'].astype(str).str.upper() == 'TRUE']['drug_name'].nunique()
        
    avg_score = gene_ints['interaction_score'].mean() if not gene_ints.empty else 0
    num_cats = len(gene_cats)
    
    # Calculate the weighted score
    score = (num_approved_drugs * 3) + (num_total_drugs * 1) + (num_cats * 2) + (avg_score * 5)
    
    druggability_list.append({
        'Gene': gene,
        'Approved_Drugs': num_approved_drugs,
        'Total_Drugs': num_total_drugs,
        'Num_Categories': num_cats,
        'Avg_Interaction_Score': round(avg_score, 3),
        'Druggability_Score': round(score, 2),
        'Categories': ", ".join(gene_cats) if num_cats > 0 else "None"
    })

# 5. SORT AND IDENTIFY TOP 10
# ------------------------------------------------------------------------------
druggability_df = pd.DataFrame(druggability_list).sort_values(by='Druggability_Score', ascending=False)
top_10 = druggability_df.head(10)

# 6. EXPORT RESULTS
# ------------------------------------------------------------------------------
druggability_df.to_csv('CKM_Gene_Druggability_Full_List.csv', index=False)
top_10.to_csv('CKM_Top_10_Druggable_Genes.csv', index=False)

# Filter interactions for only approved drugs to clean up your roadmap list
if ckm_interactions['approved'].dtype == 'bool':
    approved_interactions = ckm_interactions[ckm_interactions['approved'] == True]
else:
    approved_interactions = ckm_interactions[ckm_interactions['approved'].astype(str).str.upper() == 'TRUE']
approved_interactions.to_csv('CKM_Approved_Interactions_List.csv', index=False)

print("-" * 50)
print(f"✅ ANALYSIS COMPLETE")
print(f"Full results saved to: CKM_Gene_Druggability_Full_List.csv")
print(f"Top 10 summary saved to: CKM_Top_10_Druggable_Genes.csv")
print("-" * 50)

print("\n🏆 TOP 10 MOST DRUGGABLE CKM HUB GENES:")
print(top_10[['Gene', 'Approved_Drugs', 'Total_Drugs', 'Druggability_Score']].to_string(index=False))