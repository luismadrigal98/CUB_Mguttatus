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
import random

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

def validate_mapping(df_original, df_final_matrix, gene_map, num_checks=5):
    """
    Validates the mapping and aggregation process by tracking random random points 
    back to the original data.
    """
    print("\n--- Performing Quality Control / Validation ---")
    
    # Invert gene map to find original genes for a remapped gen (IM767 -> [IM62, IM62...])
    reverse_map = {}
    for original, remapped in gene_map.items():
        if remapped not in reverse_map:
            reverse_map[remapped] = []
        reverse_map[remapped].append(original)
        
    # Get list of remapped genes and samples
    remapped_genes = df_final_matrix.index.tolist()
    samples = df_final_matrix.columns.tolist()
    
    checks_passed = 0
    
    print(f"Checking {num_checks} random data points...")

    for i in range(num_checks):
        # 1. Select random sample
        sample = random.choice(samples)
        
        # 2. Select random gene that has expression in this sample (prefer non-zero for meaningful check)
        expressed_genes = df_final_matrix.index[df_final_matrix[sample] > 0].tolist()
        
        if not expressed_genes:
            target_gene = random.choice(remapped_genes) # Fallback if sample is empty
        else:
            target_gene = random.choice(expressed_genes)
            
        observed_val = df_final_matrix.loc[target_gene, sample]
        
        # 3. Find contributing original genes (IM62 IDs)
        if target_gene in reverse_map:
            original_sources = reverse_map[target_gene]
        else:
            print(f"Check {i+1}: FAIL - Gene {target_gene} not found in reverse map.")
            continue
            
        # 4. Retrieve values from original dataframe
        # We assume df_original has columns ['Sample', 'Gene', 'CPM']
        # Filter for Sample AND (Gene is in original_sources)
        subset = df_original[
            (df_original['Sample'] == sample) & 
            (df_original['Gene'].isin(original_sources))
        ]
        
        expected_val = subset['CPM'].sum()
        
        # 5. Check equality (with floating point tolerance)
        match = abs(observed_val - expected_val) < 1e-6
        
        status = "PASS" if match else "FAIL"
        if match: checks_passed += 1
        
        print(f"\n[Check {i+1}]: {status}")
        print(f"  Sample: {sample}")
        print(f"  Target Gene (IM767): {target_gene}")
        print(f"  Final Value (Matrix): {observed_val:.6f}")
        print(f"  Contributing Original Genes (IM62): {len(subset)} found out of {len(original_sources)} potential mapping sources")
        if not subset.empty:
            for _, row in subset.iterrows():
                print(f"    - {row['Gene']}: {row['CPM']:.6f}")
        else:
            print("    (No expression found for contributing genes in source)")
        print(f"  Expected Sum: {expected_val:.6f}")
        print(f"  Difference: {abs(observed_val - expected_val):.2e}")

    print("-" * 40)
    print(f"Validation Summary: {checks_passed}/{num_checks} checks passed.")
    if checks_passed != num_checks:
        print("WARNING: Some validation checks failed! Please review.")
    else:
        print("SUCCESS: Data compiling looks consistent.")
    print("-" * 40)

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

    # Explicitly filter out any Empty Strings in Remapped_Gene (which can cause the "empty" row issue)
    df_filtered = df_filtered[df_filtered['Remapped_Gene'] != ""]
    
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
    
    # --- 4. Quality Control ---
    # Pass the unfiltered df (but we need to manually filter sample inside validation to match context)
    # The validation function handles filtering df by sample.
    validate_mapping(df, df_matrix, final_gene_map, num_checks=5)

    # --- 5. Save Output ---
    print(f"\n--- Saving to {args.output_file} ---")
    df_matrix.to_csv(args.output_file, sep='\t')
    print("Done successfully.")

if __name__ == "__main__":
    main()
