'''
Script to compile expression data from multiple plant sources and tissues into a single DataFrame.
This script is specifically designed to work with large expression matrices (CPM) containing
IM62 and IM767 samples. It filters for relevant columns (IM samples) and translates
IM62 gene names to the IM767 reference genome.

@author: Luis Javier Madrigal-Roca
@date: 2024-06-10 (Updated 2026-01-14 by Co-pilot)
'''

import pandas as pd
import os
import sys
import argparse
import re

def build_gene_name_dictionary_IM62(GO_based_txt):
    """
    Builds a dictionary mapping IM62 identifiers to IM62 v2 identifiers.
    Source: GO_based_txt file.
    Output: { 'Migut.01G000100': 'MiIM6v31000001m', ... }
    """
    gene_name = {}
    with open(GO_based_txt, 'r') as infile:
        for idx, line in enumerate(infile):
            if idx == 0:
                continue
            parts = line.strip().split('\t')
            if len(parts) > 2:
                # Remove version suffix from key if present (Migut.01G000100.v3.1 -> Migut.01G000100)
                key = parts[0].rsplit('.', 2)[0]
                # Value is the v2 name (MiIM6v31000001m)
                gene_name[key] = parts[1].split('.')[0]
    return gene_name

def build_gene_name_dictionary_from_refs(crossref_file):
    """
    Builds a dictionary mapping IM62 identifiers to IM767 identifiers.
    Source: Cross-reference CSV.
    Output: { 'Migut.01G000300': 'MgIM767.01G000300', ... }
    """
    gene_name_62_to_767 = {}
    ommited_count = 0
    with open(crossref_file, 'r') as infile:
        for idx, line in enumerate(infile):
            parts = line.strip().split(',')
            
            # Omit header and non_IM62 references
            if idx == 0 or parts[0] != "IM62":
                continue
            
            # Avoid complex mappings for now (containing |)
            if len(parts) == 7 and "|" not in parts[4]:
                gene_name_62_to_767[parts[1]] = parts[4]
            else:
                ommited_count += 1 
    
    print(f"Total mapped genes from IM62 to IM767: {len(gene_name_62_to_767)}")
    print(f"Omitted mappings (complex/multi): {ommited_count}")

    return gene_name_62_to_767

def main():
    parser = argparse.ArgumentParser(description="Compile and remap expression data from multiple sources (IM62/IM767).")
    parser.add_argument("--input_file", required=True, help="Input expression data file (CPM matrix).")
    parser.add_argument("--output_file", required=True, help="Path for the output compiled data.")
    parser.add_argument("--crossref_file", required=True, help="Cross-reference CSV file for gene mapping.")
    parser.add_argument("--go_based_txt", required=True, help="GO based annotation TXT file for gene mapping.")
    parser.add_argument("--column_pattern", default="^IM(62|767)", help="Regex pattern to select sample columns (default: starts with IM62 or IM767).")

    args = parser.parse_args()

    # --- 1. Build Mapping Dictionaries ---
    print("\n--- Building Gene Mapping Dictionaries ---")
    
    # Dict 1: IM62 (Migut...) -> IM62_v2 (MiIM6...)
    # NOTE: In bud_expression_compiler, the dict was built as {Migut... : MiIM6...}
    # But later iterated as `for im62, im62_v2 in gene_name_dict_stage1.items()`
    # We want to map FROM the file's ID. 
    # If the file has IM62 v2 IDs (MiIM6...), we need MiIM6 -> IM767.
    dict_stage1 = build_gene_name_dictionary_IM62(args.go_based_txt)
    
    # Dict 2: IM62 (Migut...) -> IM767 (MgIM767...)
    dict_stage2 = build_gene_name_dictionary_from_refs(args.crossref_file)

    # Combine: MiIM6... (file ID) -> [via Migut...] -> MgIM767... (target ID)
    final_gene_map = {}
    linked_count = 0
    
    for im62, im62_v2 in dict_stage1.items():
        if im62 in dict_stage2:
            target_id = dict_stage2[im62]
            final_gene_map[im62_v2] = target_id
            linked_count += 1
            
    print(f"Final mapping dictionary (IM62_v2 -> IM767): {len(final_gene_map)} entries")

    # --- 2. Read and Filter Expression Data ---
    print(f"\n--- Reading Expression Data: {args.input_file} ---")
    
    try:
        # Detected as Long Format based on user input (Sample, Gene, CPM)
        # No header
        df = pd.read_csv(args.input_file, sep='\t', header=None, names=['Sample', 'Gene', 'CPM'])
        print(f"Full dataset shape: {df.shape}")
        
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    # Filter Rows based on Sample Column
    print(f"Filtering rows where Sample matches pattern: '{args.column_pattern}'")
    
    # Filter rows
    df_filtered = df[df['Sample'].str.contains(args.column_pattern, regex=True, na=False)].copy()
    
    if df_filtered.empty:
        print("ERROR: No rows matched the pattern.")
        print("First 5 samples in file:", df['Sample'].head().tolist())
        sys.exit(1)
        
    unique_samples = df_filtered['Sample'].unique()
    print(f"Found {len(unique_samples)} matching samples.")
    print("Examples:", unique_samples[:5])

    # --- 3. Remap Gene Rows and Aggregate ---
    print("\n--- Remapping Gene Rows ---")
    
    # Map Genes
    df_filtered['Remapped_Gene'] = df_filtered['Gene'].map(final_gene_map)
    
    # Report mapping stats
    total_rows = len(df_filtered)
    mapped_rows = df_filtered['Remapped_Gene'].notna().sum()
    print(f"Rows with valid mapping: {mapped_rows} / {total_rows} ({mapped_rows/total_rows*100:.1f}%)")
    
    # Drop unmapped
    df_filtered = df_filtered.dropna(subset=['Remapped_Gene'])
    
    # Aggregate (Sum CPM for duplicate mappings within the same sample)
    # Group by [Sample, Remapped_Gene] -> Sum CPM
    print("Aggregating duplicates (summing CPM)...")
    df_grouped = df_filtered.groupby(['Remapped_Gene', 'Sample'])['CPM'].sum().reset_index()
    
    # Pivot to Wide Format (Genes x Samples)
    print("Pivoting to wide format (Matrix)...")
    df_matrix = df_grouped.pivot(index='Remapped_Gene', columns='Sample', values='CPM')
    
    # Fill NaN with 0 (missing gene in a sample = 0 expression)
    df_matrix = df_matrix.fillna(0)
    
    print(f"Final Matrix Shape: {df_matrix.shape}")
    
    # --- 4. Save Output ---
    print(f"\n--- Saving to {args.output_file} ---")
    df_matrix.to_csv(args.output_file, sep='\t')
    print("Done successfully.")

if __name__ == "__main__":
    main()
