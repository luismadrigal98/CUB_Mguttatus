#!/usr/bin/env python3
"""
Quick diagnostic tool to check chromosome names in input files.

This helps identify mismatches before running the full pipeline.

Usage:
    python check_chromosome_names.py <gff3> <genome_fasta> <cds_fasta>

Author: Luis Javier Madrigal-Roca & John K. Kelly
"""

import sys

def check_genome_fasta(fasta_file):
    """Extract chromosome names from genome FASTA."""
    print(f"\n📄 Checking: {fasta_file}")
    chroms = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line.strip()[1:].split()[0]
                chroms.append(header)
    
    print(f"   Found {len(chroms)} sequences")
    print(f"   First 10 chromosome names:")
    for chrom in chroms[:10]:
        print(f"     {chrom}")
    if len(chroms) > 10:
        print(f"     ... and {len(chroms) - 10} more")
    
    return set(chroms)

def check_cds_fasta(fasta_file):
    """Extract gene IDs from CDS FASTA."""
    print(f"\n📄 Checking: {fasta_file}")
    genes = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line.strip()[1:].split()[0]
                # Extract simplified gene ID (remove transcript suffix)
                if '.' in header:
                    gene_id = '.'.join(header.split('.')[:2])  # MgIM767.01G000100
                    gene_id = gene_id.split('.')[1]  # 01G000100
                else:
                    gene_id = header
                genes.append(gene_id)
    
    print(f"   Found {len(genes)} CDS sequences")
    print(f"   First 10 gene IDs (after simplification):")
    for gene in genes[:10]:
        print(f"     {gene}")
    if len(genes) > 10:
        print(f"     ... and {len(genes) - 10} more")
    
    return set(genes)

def check_gff3(gff3_file):
    """Extract chromosome names and gene IDs from GFF3."""
    print(f"\n📄 Checking: {gff3_file}")
    chroms = set()
    genes_by_chrom = {}
    
    with open(gff3_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            
            chr_name = cols[0]
            feature = cols[2]
            attrs = cols[8]
            
            chroms.add(chr_name)
            
            if feature == 'CDS':
                if chr_name not in genes_by_chrom:
                    genes_by_chrom[chr_name] = set()
                
                # Extract gene ID
                attr_dict = {}
                for attr in attrs.split(';'):
                    if '=' in attr:
                        key, val = attr.split('=', 1)
                        attr_dict[key] = val
                
                if 'Parent' in attr_dict:
                    parent = attr_dict['Parent']
                    gene_id = parent.split('.')[1] if '.' in parent else parent
                    genes_by_chrom[chr_name].add(gene_id)
    
    print(f"   Found {len(chroms)} chromosomes")
    print(f"   First 10 chromosome names:")
    for chrom in sorted(chroms)[:10]:
        n_genes = len(genes_by_chrom.get(chrom, set()))
        print(f"     {chrom} ({n_genes} genes with CDS)")
    if len(chroms) > 10:
        print(f"     ... and {len(chroms) - 10} more")
    
    return chroms, genes_by_chrom

def main():
    if len(sys.argv) != 4:
        print("Usage: python check_chromosome_names.py <gff3> <genome_fasta> <cds_fasta>")
        print("\nExample:")
        print("  python check_chromosome_names.py \\")
        print("    data/Mguttatusvar_IM767_887_v2.1.gene.gff3 \\")
        print("    data/Mguttatusvar_IM767_887_v2.0.fa \\")
        print("    data/Mguttatusvar_IM767_887_v2.1.cds_primaryTranscriptOnly.fa")
        sys.exit(1)
    
    gff3_file = sys.argv[1]
    genome_fasta = sys.argv[2]
    cds_fasta = sys.argv[3]
    
    print("=" * 70)
    print("Chromosome Name Diagnostic Tool")
    print("=" * 70)
    
    # Check all files
    genome_chroms = check_genome_fasta(genome_fasta)
    cds_genes = check_cds_fasta(cds_fasta)
    gff3_chroms, gff3_genes_by_chrom = check_gff3(gff3_file)
    
    # Compare
    print("\n" + "=" * 70)
    print("COMPARISON")
    print("=" * 70)
    
    # Chromosomes in genome but not GFF3
    missing_in_gff3 = genome_chroms - gff3_chroms
    if missing_in_gff3:
        print(f"\n⚠️  Chromosomes in genome FASTA but NOT in GFF3:")
        for chrom in sorted(missing_in_gff3)[:5]:
            print(f"   {chrom}")
        if len(missing_in_gff3) > 5:
            print(f"   ... and {len(missing_in_gff3) - 5} more")
    
    # Chromosomes in GFF3 but not genome
    missing_in_genome = gff3_chroms - genome_chroms
    if missing_in_genome:
        print(f"\n⚠️  Chromosomes in GFF3 but NOT in genome FASTA:")
        for chrom in sorted(missing_in_genome)[:5]:
            print(f"   {chrom}")
        if len(missing_in_genome) > 5:
            print(f"   ... and {len(missing_in_genome) - 5} more")
    
    # Common chromosomes
    common_chroms = genome_chroms & gff3_chroms
    if common_chroms:
        print(f"\n✓ {len(common_chroms)} chromosomes present in BOTH genome and GFF3:")
        for chrom in sorted(common_chroms)[:10]:
            n_genes = len(gff3_genes_by_chrom.get(chrom, set()))
            print(f"   {chrom} ({n_genes} genes)")
        if len(common_chroms) > 10:
            print(f"   ... and {len(common_chroms) - 10} more")
    
    # Gene comparison (sample)
    if common_chroms:
        sample_chrom = sorted(common_chroms)[0]
        sample_genes = list(gff3_genes_by_chrom.get(sample_chrom, set()))[:5]
        
        print(f"\n✓ Checking if sample genes from {sample_chrom} exist in CDS FASTA:")
        for gene in sample_genes:
            if gene in cds_genes:
                print(f"   {gene} ✓")
            else:
                print(f"   {gene} ✗ (NOT FOUND)")
    
    print("\n" + "=" * 70)
    print("RECOMMENDATIONS")
    print("=" * 70)
    
    if not missing_in_gff3 and not missing_in_genome and common_chroms:
        print("✓ All files look compatible!")
        print("\nYou can process these chromosomes:")
        for chrom in sorted(common_chroms)[:14]:
            print(f"  {chrom}")
        
        print("\nExample command:")
        sample_chrom = sorted(common_chroms)[0]
        print(f"  python describe_gene_positions_by_degeneracy.py {sample_chrom} \\")
        print(f"    {gff3_file} \\")
        print(f"    {genome_fasta} \\")
        print(f"    {cds_fasta}")
    else:
        print("⚠️  Found some mismatches!")
        print("   Check that all files are from the same genome version.")
        print("   Chromosome names must match exactly (case-sensitive).")
    
    print()

if __name__ == "__main__":
    main()
