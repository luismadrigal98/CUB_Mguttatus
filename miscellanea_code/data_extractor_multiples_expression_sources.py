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
    
    # Read CSV (assuming tab separated based on context, header=0, index_col=0 for Genes)
    try:
        # Check if file is tab or comma separated
        with open(args.input_file, 'r') as f:
            first_line = f.readline()
            sep = '\t' if '\t' in first_line else ','
            
        print(f"Detected separator: '{'tab' if sep=='\t' else 'comma'}'")
        df = pd.read_csv(args.input_file, sep=sep, header=0, index_col=0)
        print(f"Full dataset shape: {df.shape}")
        
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    # Filter Columns
    print(f"Filtering columns with pattern: '{args.column_pattern}'")
    im_cols = [c for c in df.columns if re.search(args.column_pattern, c)]
    
    if not im_cols:
        print("ERROR: No columns matched the pattern. Check your input file header or pattern.")
        print("First 10 columns in file:", df.columns[:10].tolist())
        sys.exit(1)
        
    print(f"Found {len(im_cols)} matching columns.")
    # Example columns
    print("Examples:", im_cols[:5])
    
    df_filtered = df[im_cols]

    # --- 3. Remap Gene Rows and Aggregate ---
    print("\n--- Remapping Gene Rows ---")
    
    # Create new column for mapped ID
    df_filtered = df_filtered.copy()
    df_filtered['Remapped_Gene'] = df_filtered.index.map(final_gene_map)
    
    # Report mapping stats
    total_rows = len(df_filtered)
    mapped_rows = df_filtered['Remapped_Gene'].notna().sum()
    print(f"Rows with valid mapping: {mapped_rows} / {total_rows} ({mapped_rows/total_rows*100:.1f}%)")
    
    # Drop unmapped
    df_filtered = df_filtered.dropna(subset=['Remapped_Gene'])
    
    # Check for duplicates before summing
    n_unique_targets = df_filtered['Remapped_Gene'].nunique()
    print(f"Unique target genes: {n_unique_targets}")
    if mapped_rows > n_unique_targets:
        print(f"Collapsing {mapped_rows - n_unique_targets} rows by summing (multiple input genes -> same target).")
    
    # Group by new ID and SUM
    # "This is correct for CPM because we're combining reads that map to the same locus"
    df_final = df_filtered.groupby('Remapped_Gene').sum()
    
    # --- 4. Save Output ---
    print(f"\n--- Saving to {args.output_file} ---")
    df_final.to_csv(args.output_file, sep='\t')
    print("Done successfully.")

if __name__ == "__main__":
    main()
