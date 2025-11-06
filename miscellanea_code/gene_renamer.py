"""
This script is designed for changing the names of the genes in a txt file.

@Author: Luis Javier Madrigal Roca & John K. Kelly
@Date: 2024-09-10

"""

import pandas as pd
import os

# Expression data directory
input_file = '/mnt/1692B2EF92B2D28B/Ongoing_projects/Codon_bias_analysis/data/leaf_gene_expression_mean_cpm.txt'
output_file = '/mnt/1692B2EF92B2D28B/Ongoing_projects/Codon_bias_analysis/data/leaf_gene_expression_mean_cpm_renamed.txt'

# Read the Excel file
df = pd.read_excel('/mnt/1692B2EF92B2D28B/Ongoing_projects/Codon_bias_analysis/data/new_name_old_name.genes.xlsx')

# Create a dictionary mapping old names to new names
name_mapping = df.set_index('old_name')['new_name'].to_dict()

# List all the lines of the input file
with open(input_file, "r") as file:
    lines = file.readlines()

# Open the output file in write mode
with open(output_file, "w") as out_file:
    for line in lines:
        # Extract the old gene name from the line (assuming it's the first element in a tab-separated line)
        old_gene_name = line.split('\t')[0]
        # Check if the old gene name is in the mapping
        if old_gene_name in name_mapping:
            # Get the new gene name from the mapping
            new_gene_name = name_mapping[old_gene_name]
            # Replace the old gene name with the new gene name in the line
            line = line.replace(old_gene_name, new_gene_name)
        # Write the modified line to the output file
        out_file.write(line)