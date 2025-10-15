"""
This script filters tRNA gene categorized as possible pseudogenes.

@author: Luis Javier Madrigal Roca

@date: 2025-10-15

"""

import pandas as pd
import sys

def filter_tRNA(input_file, output_file):
    # Read the input tRNA file, skipping the first 3 header lines
    df = pd.read_csv(input_file, skiprows=3, delim_whitespace=True, header=None)

    # Filter out rows where the 4th column is 'pseudogene'
    filtered_df = df[df[9] != 'pseudogene']

    filtered_df.columns = ['Sequence_name', 'tRNA_number', 'tRNA_start', 'tRNA_end', 'tRNA_type', 'Anticodon', 'Intron_start', 'Intron_end', 'Score', 'Notes']

    # Write the filtered DataFrame to the output file
    filtered_df.to_csv(output_file, sep='\t', index=False)

def main():
    if len(sys.argv) != 3:
        print("Usage: python filter_tRNA.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    filter_tRNA(input_file, output_file)

if __name__ == "__main__":
    main()