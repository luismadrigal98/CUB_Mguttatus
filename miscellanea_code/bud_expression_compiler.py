#!/usr/bin/env python3

"""
Compiler of expression data from bud tissue samples from Mimulus guttatus

- Read expression data files
- Remap gene names from IM62 to IM767 reference genome

@author: Luis Javier Madrigal Roca
@date: 2025-10-24

"""

import pandas as pd
import os
import sys
import re

# Appendixs2
# refGenome,refGenomeID,IM62,IM62_v2,IM767,SF,LVR
# LVR,Mitil.01G000100,Migut.01G000300,Migut.A00001,MgIM767.01G000300,Minas.01G000100,NA

def build_gene_name_dictionary_IM62(GO_based_txt):
    #G_name	mg_name	chrom	stpos	endpos	direction	GOTERMs_>									
    #Migut.01G000100.v3.1	MiIM6v31000001m.g	Chr_01	18808	21531	+
    
    gene_name = {}
    with open(GO_based_txt, 'r') as infile:
        for idx, line in enumerate(infile):
            if idx == 0:
                continue
            parts = line.strip().split('\t')
            if len(parts) > 2:
                # Remove version suffix from key if present
                key = parts[0].rsplit('.', 2)[0]
                gene_name[key] = parts[1].split('.')[0]
    return gene_name

def build_gene_name_dictionary_from_refs(crossref_file):
    #refGenome,refGenomeID,IM62,IM62_v2,IM767,SF,LVR
    #LVR,Mitil.01G000100,Migut.01G000300,Migut.A00001,MgIM767.01G000300,Minas.01G000100,NA
    #LVR,Mitil.01G000200,Migut.01G000400,Migut.A00002,MgIM767.O003200|MgIM767.01G000400,Minas.01G000200,NA

    gene_name_62_to_767 = {}
    ommited_count = 0
    with open(crossref_file, 'r') as infile:
        for idx, line in enumerate(infile):
            parts = line.strip().split(',')
            
            # Omit header and non_IM62 references
            if idx == 0 or parts[0] != "IM62":
                continue
            
            # Avoid complex mappings for now
            if len(parts) == 7 and "|" not in parts[4]:
                gene_name_62_to_767[parts[1]] = parts[4]
            else:
                ommited_count += 1 
    
    ## DEBUGGING:
    print(f"Total mapped genes: {len(gene_name_62_to_767)}")
    print(f"Omitted mappings: {ommited_count}")

    return gene_name_62_to_767

def remap_expression_file(input_expression_file, output_expression_file, gene_name_dict):
    '''
    This function remaps gene names in the expression data file and sums expression
    values for genes that map to the same IM767 gene.

    @param input_expression_file: Path to the input expression data file
    @param output_expression_file: Path to the output remapped expression data file
    @param gene_name_dict: Dictionary mapping gene names from 62 to 767 reference.

    '''
    
    #MiIM6v31036764m	23.4928028061896
    #MiIM6v31036765m	39.38709663540622
    #MiIM6v31036766m	38.33155733102882

    expression_data = pd.read_csv(input_expression_file, sep='\t', header=None, names=['Gene', 'Expression'])

    # Add a column for remapped gene names
    expression_data['Remapped_Gene'] = expression_data['Gene'].map(gene_name_dict)

    total_genes = len(expression_data)
    expression_data = expression_data.dropna(subset=['Remapped_Gene'])
    remapped_genes = len(expression_data)
    print(f"Total genes in original file: {total_genes}")
    print(f"Total genes after remapping: {remapped_genes}")

    # Check for duplicates BEFORE summing
    duplicate_check = expression_data.groupby('Remapped_Gene').size()
    duplicated_genes = duplicate_check[duplicate_check > 1]
    
    if len(duplicated_genes) > 0:
        print(f"\nFound {len(duplicated_genes)} genes with multiple expression entries")
        print("Examples (first 5):")
        for gene in duplicated_genes.head(5).index:
            subset = expression_data[expression_data['Remapped_Gene'] == gene]
            print(f"  {gene}: {len(subset)} entries with values {subset['Expression'].tolist()}")
    
    # Sum expression values for duplicate gene mappings
    # This is correct for CPM because we're combining reads that map to the same locus
    expression_data_summed = expression_data.groupby('Remapped_Gene', as_index=False)['Expression'].sum()
    
    final_genes = len(expression_data_summed)
    duplicates_collapsed = remapped_genes - final_genes
    
    if duplicates_collapsed > 0:
        print(f"Collapsed {duplicates_collapsed} duplicate mappings by summing expression values")
    
    print(f"Final unique genes: {final_genes}")
    
    # Rename column for clarity
    expression_data_summed.columns = ['Gene', 'Expression']
    
    # Final sanity check AFTER renaming
    if len(expression_data_summed) != len(expression_data_summed['Gene'].unique()):
        print("\n*** ERROR: Output still contains duplicate genes! ***")
    
    expression_data_summed.to_csv(output_expression_file, sep='\t', index=False)

def main():
    # Store arguments

    if len(sys.argv) != 5:
        print("Usage: python bud_expression_compiler.py <input_expression_file> <output_expression_file> <crossref_file> <GO_based_txt>")
        sys.exit(1)

    input_expression_file = sys.argv[1]
    output_expression_file = sys.argv[2]
    crossref_file = sys.argv[3]
    GO_based_txt = sys.argv[4]

    # Build gene name mapping dictionary stage 1 (IM62 v2.0 to IM62)
    gene_name_dict_stage1 = build_gene_name_dictionary_IM62(GO_based_txt)

    # Build gene name mapping dictionary stage 2 (IM62 to IM767)
    gene_name_dict_stage2 = build_gene_name_dictionary_from_refs(crossref_file)

    # Combine both dictionaries: IM62_v2.0 -> IM62 -> IM767
    # For each gene in expression file, map to IM62, then to IM767
    final_gene_name_dict = {}
    duplicates_found = {}
    
    for im62, im62_v2 in gene_name_dict_stage1.items():
        im767 = gene_name_dict_stage2.get(im62)
        if im767:
            # Check if this IM62_v2 already exists in the final dict
            if im62_v2 in final_gene_name_dict:
                if im62_v2 not in duplicates_found:
                    duplicates_found[im62_v2] = []
                duplicates_found[im62_v2].append((im62, im767))
            else:
                final_gene_name_dict[im62_v2] = im767
    
    # Report duplicates
    if duplicates_found:
        print(f"\nWARNING: Found {len(duplicates_found)} IM62_v2 genes mapping to multiple IM62/IM767 genes")
        print("Examples (first 5):")
        for i, (im62_v2, mappings) in enumerate(list(duplicates_found.items())[:5]):
            print(f"  {im62_v2} -> {mappings}")
    
    # Check for duplicate values (multiple IM62_v2 -> same IM767)
    reverse_dict = {}
    for im62_v2, im767 in final_gene_name_dict.items():
        if im767 not in reverse_dict:
            reverse_dict[im767] = []
        reverse_dict[im767].append(im62_v2)
    
    many_to_one = {k: v for k, v in reverse_dict.items() if len(v) > 1}
    if many_to_one:
        print(f"\nWARNING: Found {len(many_to_one)} IM767 genes mapped from multiple IM62_v2 genes")
        print("Examples (first 5):")
        for i, (im767, im62_v2_list) in enumerate(list(many_to_one.items())[:5]):
            print(f"  {im767} <- {im62_v2_list}")

    # Remap expression file
    remap_expression_file(input_expression_file, output_expression_file, final_gene_name_dict)

if __name__ == "__main__":
    main()