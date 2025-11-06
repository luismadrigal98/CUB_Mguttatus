#!/usr/bin/env python3
"""
Compile expression data from bud tissue samples.
Reads tab-delimited transcript.*.filelist.txt files containing gene counts
from 767 crosses with different lines, merges them, and normalizes to CPM.

Usage:
    python leaf_expression_compiler.py <input_dir> <output_file>
    
Where input_dir contains transcript.*.filelist.txt files
"""

import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path
import re

def find_transcript_files(input_dir):
    """
    Find all transcript.*.filelist.txt files in the input directory.
    Returns a list of file paths.
    """
    input_path = Path(input_dir)
    
    # Find all files matching pattern transcript.*.filelist.txt
    transcript_files = list(input_path.glob("transcript.*.filelist.txt"))
    
    return transcript_files

def read_expression_file(file_path):
    """
    Read tab-delimited expression file.
    First column is gene_id, remaining columns are sample counts.
    Adds file identifier to column names to distinguish plants from different crosses.
    """
    print(f"  Reading: {file_path.name}")
    
    # Extract the cross identifier from filename (e.g., "62" from "transcript.62.filelist.txt")
    match = re.match(r'transcript\.(\d+)\.filelist\.txt', file_path.name)
    if match:
        cross_id = match.group(1)
    else:
        cross_id = file_path.stem  # fallback to full filename
    
    # Read the file
    df = pd.read_csv(file_path, sep='\t')
    
    # Ensure gene_id column exists
    if 'gene_id' not in df.columns:
        print(f"    WARNING: No 'gene_id' column found in {file_path.name}")
        return None
    
    # Rename columns to include cross identifier (except gene_id)
    # This ensures plants from different crosses are kept separate
    new_columns = {'gene_id': 'gene_id'}
    for col in df.columns:
        if col != 'gene_id':
            new_columns[col] = f"cross{cross_id}_{col}"
    
    df.rename(columns=new_columns, inplace=True)
    
    print(f"    Genes: {len(df)}, Samples: {len(df.columns) - 1}")
    
    return df

def sum_allele_counts(df):
    """
    Sum allele-specific counts from the same plant.
    Within each cross, columns like 'cross62_767_s1_767-P1' and 'cross62_62_s1_767-P1' 
    should be summed into 'cross62_s1_767-P1' representing total expression for that plant.
    """
    print("\nSumming allele-specific counts per plant...")
    
    count_columns = [col for col in df.columns if col != 'gene_id']
    
    # Dictionary to store summed columns: {plant_id: [col1, col2, ...]}
    plant_groups = {}
    
    for col in count_columns:
        # Extract plant identifier from cross-prefixed columns
        # Pattern: crossXXX_LINENUM_sX_PLANT-ID
        # Examples: cross62_767_s1_767-P1 -> cross62_s1_767-P1
        #           cross62_62_s1_767-P1 -> cross62_s1_767-P1
        match = re.match(r'(cross\d+)_(\d+)_(s\d+_.+)', col)
        if match:
            cross_id = match.group(1)  # cross62
            plant_id = match.group(3)  # s1_767-P1
            full_plant_id = f"{cross_id}_{plant_id}"  # cross62_s1_767-P1
            
            if full_plant_id not in plant_groups:
                plant_groups[full_plant_id] = []
            plant_groups[full_plant_id].append(col)
        else:
            # If pattern doesn't match, keep column as is
            plant_groups[col] = [col]
    
    print(f"  Found {len(plant_groups)} unique plants across all crosses")
    print(f"  Original columns: {len(count_columns)}")
    
    # Create new dataframe with summed counts
    summed_data = {'gene_id': df['gene_id']}
    
    for plant_id, cols in plant_groups.items():
        if len(cols) > 1:
            # Sum allele counts for this plant
            # Use sum(axis=1, min_count=1) to return NA if all values are NA
            # This preserves NA for genes not present in this cross
            summed_data[plant_id] = df[cols].sum(axis=1, min_count=1)
        else:
            # Only one column for this plant (no allele splitting)
            summed_data[plant_id] = df[cols[0]]
    
    summed_df = pd.DataFrame(summed_data)
    
    print(f"  New columns (after summing alleles): {len(summed_df.columns) - 1}")
    
    return summed_df

def compile_expression_data(input_dir, output_file):
    """
    Main function to compile expression data from all transcript files.
    """
    print(f"Searching for transcript.*.filelist.txt files in: {input_dir}")
    
    # Find all transcript files
    transcript_files = find_transcript_files(input_dir)
    
    if not transcript_files:
        print("ERROR: No transcript.*.filelist.txt files found!")
        return
    
    print(f"\nFound {len(transcript_files)} transcript files")
    
    # Read all files
    all_data = []
    
    for file_path in sorted(transcript_files):
        df = read_expression_file(file_path)
        if df is not None:
            all_data.append(df)
    
    if not all_data:
        print("ERROR: No valid data files could be read!")
        return
    
    # Merge all dataframes
    print("\nMerging all files...")
    print("  Using outer join to include all genes from all files")
    print("  Note: Plants from different crosses are kept separate with cross identifiers")
    
    merged_df = all_data[0]
    for i, df in enumerate(all_data[1:], 1):
        print(f"  Merging file {i+1}/{len(all_data)}...")
        # Simple merge - no overlapping columns since we added cross identifiers
        merged_df = merged_df.merge(df, on='gene_id', how='outer')
    
    # Keep NA for genes not present in some crosses (don't fill with 0)
    print("\nKeeping NA for genes not present in some crosses...")
    print("  (This allows proper mean/median calculation ignoring missing genes)")
    
    # Get count columns (all except gene_id)
    count_columns = [col for col in merged_df.columns if col != 'gene_id']
    
    # Convert counts to float (they may have decimal values from salmon)
    # Keep NAs for truly missing values
    print("Converting counts to numeric values...")
    for col in count_columns:
        merged_df[col] = pd.to_numeric(merged_df[col], errors='coerce')
    
    # Sum allele-specific counts per plant
    merged_df = sum_allele_counts(merged_df)
    
    # Get count columns (all except gene_id) after summing
    count_columns = [col for col in merged_df.columns if col != 'gene_id']
    
    # Count NA values
    na_counts = merged_df[count_columns].isna().sum().sum()
    total_cells = len(merged_df) * len(count_columns)
    na_percent = (na_counts / total_cells) * 100
    
    print(f"\nFinal matrix dimensions: {merged_df.shape}")
    print(f"  Genes: {len(merged_df)}")
    print(f"  Plants (samples): {len(count_columns)}")
    print(f"  NA values: {na_counts:,} / {total_cells:,} ({na_percent:.2f}%)")
    print(f"  (NAs represent genes not present in certain cross annotations)")
    
    # Save raw counts (summed per plant)
    raw_output = output_file.replace('.csv', '_raw_counts.csv')
    print(f"\nSaving raw counts to: {raw_output}")
    merged_df.to_csv(raw_output, sep='\t', index=False)
    
    # Calculate CPM (Counts Per Million) per sample
    print("\nCalculating CPM normalization...")
    print("  (Keeping NA for genes not in cross annotations)")
    cpm_df = merged_df.copy()
    
    samples_normalized = 0
    samples_with_zero_counts = 0
    
    for col in count_columns:
        # Sum only non-NA values for this sample
        total_counts = merged_df[col].sum(skipna=True)
        if total_counts > 0:
            # Calculate CPM, keeping NA where data is missing
            cpm_df[col] = (merged_df[col] / total_counts) * 1e6
            samples_normalized += 1
        else:
            # If all values are NA or zero, keep as is
            cpm_df[col] = merged_df[col]
            samples_with_zero_counts += 1
    
    print(f"  Normalized: {samples_normalized} samples")
    if samples_with_zero_counts > 0:
        print(f"  WARNING: {samples_with_zero_counts} samples had zero total counts")
    
    # Save CPM normalized data
    cpm_output = output_file.replace('.csv', '_cpm.csv')
    print(f"Saving CPM normalized data to: {cpm_output}")
    cpm_df.to_csv(cpm_output, sep='\t', index=False)
    
    # Calculate summary statistics per gene
    print("\nCalculating summary statistics...")
    print("  (Using skipna=True to ignore NAs - genes not present in all crosses)")
    
    mean_cpm = cpm_df[count_columns].mean(axis=1, skipna=True)
    median_cpm = cpm_df[count_columns].median(axis=1, skipna=True)
    max_cpm = cpm_df[count_columns].max(axis=1, skipna=True)
    total_raw = merged_df[count_columns].sum(axis=1, skipna=True)
    
    # Count number of samples where gene is detected (non-NA and > 0)
    samples_detected = (merged_df[count_columns] > 0).sum(axis=1)
    # Count number of samples where gene annotation exists (non-NA)
    samples_present = merged_df[count_columns].notna().sum(axis=1)
    
    summary_df = pd.DataFrame({
        'gene_id': merged_df['gene_id'],
        'mean_CPM': mean_cpm,
        'median_CPM': median_cpm,
        'max_CPM': max_cpm,
        'total_raw_counts': total_raw,
        'samples_detected': samples_detected,
        'samples_present': samples_present  # New: how many crosses have this gene
    })
    
    # Sort by mean_CPM descending
    summary_df = summary_df.sort_values('mean_CPM', ascending=False)
    
    summary_output = output_file.replace('.csv', '_summary.csv')
    print(f"Saving summary statistics to: {summary_output}")
    summary_df.to_csv(summary_output, sep='\t', index=False)
    
    # Print some summary statistics
    print("\n=== SUMMARY ===")
    print(f"Total genes: {len(merged_df)}")
    print(f"Total samples: {len(count_columns)}")
    print(f"Genes present in all {len(count_columns)} samples: {(summary_df['samples_present'] == len(count_columns)).sum()}")
    print(f"Genes expressed (>0 in at least 1 sample): {(summary_df['samples_detected'] > 0).sum()}")
    print(f"Mean CPM range: {summary_df['mean_CPM'].min():.2f} - {summary_df['mean_CPM'].max():.2f}")
    print(f"\nTop 5 most highly expressed genes (by mean CPM):")
    print(summary_df[['gene_id', 'mean_CPM', 'samples_detected', 'samples_present']].head())
    
    print(f"\n=== OUTPUT FILES ===")
    print(f"  1. Raw counts: {raw_output}")
    print(f"  2. CPM normalized: {cpm_output}")
    print(f"  3. Summary (per gene): {summary_output}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python leaf_expression_compiler.py <input_dir> <output_file>")
        print("\nExample:")
        print("  python leaf_expression_compiler.py /path/to/cases ./results/leaf_expression.csv")
        print("\nWhere input_dir contains transcript.*.filelist.txt files")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(input_dir):
        print(f"ERROR: Input directory not found: {input_dir}")
        sys.exit(1)
    
    # Create output directory if needed
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    compile_expression_data(input_dir, output_file)
