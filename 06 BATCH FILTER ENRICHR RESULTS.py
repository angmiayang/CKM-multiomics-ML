import pandas as pd
import os

# ==============================================================================
# SCRIPT: BATCH FILTER ENRICHR RESULTS
# ==============================================================================
# Purpose: Filter multiple Enrichr files for Gene Count > 3 & Adj. P-value < 0.05
# Input:   List of .txt files
# Output:  Filtered_*.csv files
# ==============================================================================

# 1. DEFINE YOUR FILE LIST
# These are the exact files you provided
file_list = [
    "GO_Biological_Process_2025_table.txt",
    "GO_Cellular_Component_2025_table.txt",
    "GO_Molecular_Function_2025_table.txt",
    "KEGG_2026_table.txt",
    "MSigDB_Hallmark_2020_table.txt",
    "WikiPathways_2024_Human_table.txt"
]

print(f"--- STARTING BATCH FILTERING FOR {len(file_list)} FILES ---\n")

for input_file in file_list:
    
    # Check if file exists to prevent crashing
    if not os.path.exists(input_file):
        print(f"⚠️  Missing: {input_file} (Skipping)")
        continue

    # 2. LOAD DATA
    # Enrichr tables are tab-separated
    df = pd.read_csv(input_file, sep='\t')
    original_count = len(df)

    # 3. EXTRACT GENE COUNT
    # Split "13/156" -> take "13" -> convert to integer
    df['Gene_Count'] = df['Overlap'].apply(lambda x: int(str(x).split('/')[0]))

    # 4. FILTER DATA
    # Criteria: Genes > 3 AND Adjusted P-value < 0.05
    filtered_df = df[
        (df['Gene_Count'] > 3) & 
        (df['Adjusted P-value'] < 0.05)
    ].copy()

    # 5. SAVE RESULT
    # We add "Filtered_" to the start of the new filename
    output_file = "Filtered_" + input_file.replace(".txt", ".csv")
    filtered_df.to_csv(output_file, index=False)

    # 6. PRINT SUMMARY
    print(f"✅ {input_file}")
    print(f"   Original: {original_count} -> Retained: {len(filtered_df)}")
    print(f"   Saved to: {output_file}\n")

print("--- BATCH PROCESSING COMPLETE ---")