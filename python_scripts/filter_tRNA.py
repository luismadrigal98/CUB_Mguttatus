"""
This script filters tRNA gene categorized as possible pseudogenes and optionally excludes scaffold entries.

@author: Luis Javier Madrigal Roca

@date: 2025-10-15

"""

import pandas as pd
import sys
import argparse

def filter_tRNA(input_file, output_file, exclude_scaffolds=False):
    # Read the input tRNA file, skipping the first 3 header lines
    # Use sep='\s+' instead of deprecated delim_whitespace
    df = pd.read_csv(input_file, skiprows=3, sep='\s+', header=None)
    
    # Handle variable number of columns by ensuring we have at least 10 columns
    # Some rows have 9 columns (no Notes), some have 10 (with Notes like 'pseudo')
    max_cols = df.shape[1]
    if max_cols < 10:
        # Add empty Notes column if it doesn't exist
        df[9] = ''
    
    # Filter out rows where the last column contains 'pseudo' (indicating pseudogenes)
    # Check both column 9 (Notes) and handle cases where 'pseudo' might be in different positions
    if max_cols >= 10:
        filtered_df = df[~df[9].astype(str).str.contains('pseudo', na=False)]
    else:
        filtered_df = df.copy()

    # Optionally filter out scaffold entries
    if exclude_scaffolds:
        # Filter out rows where sequence name contains 'scaffold' (case-insensitive)
        filtered_df = filtered_df[~filtered_df[0].astype(str).str.contains('scaffold', case=False, na=False)]

    filtered_df.columns = ['Sequence_name', 'tRNA_number', 'tRNA_start', 'tRNA_end', 'tRNA_type', 'Anticodon', 'Intron_start', 'Intron_end', 'Score', 'Notes']

    # Write the filtered DataFrame to the output file
    filtered_df.to_csv(output_file, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(
        description='Filter tRNA genes by removing pseudogenes and optionally scaffold entries.'
    )
    parser.add_argument('input_file', help='Input tRNA file path')
    parser.add_argument('output_file', help='Output filtered tRNA file path')
    parser.add_argument('--exclude-scaffolds', action='store_true', 
                       help='Exclude scaffold entries from the output')
    
    args = parser.parse_args()
    
    filter_tRNA(args.input_file, args.output_file, args.exclude_scaffolds)

if __name__ == "__main__":
    main()