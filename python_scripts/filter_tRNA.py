"""
This script filters tRNA gene categorized as possible pseudogenes and optionally excludes scaffold entries.

@author: Luis Javier Madrigal Roca

@date: 2025-10-15

"""

import pandas as pd
import sys
import argparse

def filter_tRNA(input_file, output_file, exclude_scaffolds=False):
    # Read the file line by line to handle variable number of columns
    rows = []
    with open(input_file, 'r') as f:
        # Skip the first 3 header lines
        for _ in range(3):
            next(f)
        
        # Process each data line
        for line in f:
            line = line.strip()
            if line:  # Skip empty lines
                parts = line.split()
                # Ensure we have at least 9 columns, pad with empty string if needed
                while len(parts) < 10:
                    parts.append('')
                # Take only the first 10 columns to handle any extra whitespace
                rows.append(parts[:10])
    
    # Create DataFrame with consistent structure
    df = pd.DataFrame(rows, columns=['Sequence_name', 'tRNA_number', 'tRNA_start', 'tRNA_end', 
                                   'tRNA_type', 'Anticodon', 'Intron_start', 'Intron_end', 'Score', 'Notes'])
    
    # Filter out rows where the Notes column contains 'pseudo' (indicating pseudogenes)
    filtered_df = df[~df['Notes'].str.contains('pseudo', na=False)]

    # Optionally filter out scaffold entries
    if exclude_scaffolds:
        # Filter out rows where sequence name contains 'scaffold' (case-insensitive)
        filtered_df = filtered_df[~filtered_df['Sequence_name'].str.contains('scaffold', case=False, na=False)]

    # Write the filtered DataFrame to the output file, excluding the Notes column, because at this point is empty
    filtered_df = filtered_df.drop(columns=['Notes'])
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