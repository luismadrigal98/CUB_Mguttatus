#!/usr/bin/env python3
"""
Reformat codon frequency output to be PREFERRED-codon-centric.

Takes the reference-based codon frequency file and reformats it so that:
1. Preferred_Codon column shows the optimal codon for that amino acid
2. Codon_Variants always includes the preferred codon (even if count=0)
3. Frequencies shows preferred codon frequency first
4. Easy to calculate selection metrics (preferred vs non-preferred)

This makes it trivial to:
- Calculate frequency of preferred codon usage per gene
- Measure selection strength (how often preferred is used)
- Compare polymorphism at preferred vs non-preferred sites
- Test for selection on codon usage bias

Usage:
    python reformat_codon_freq_by_preference.py \\
        <ref_based_codon_freq.txt> \\
        <preferred_codons.txt> \\
        <output_file.txt>

Input format (reference-based):
    Gene    Codon_Pos    AA    Ref_Codon    Is_Preferred    Genomic_Positions    Codon_Variants    Frequencies
    01G000100    2    I    ATT    True    14530,14531,14532    ATT:186;ATC:1    ATT:0.995;ATC:0.005

Output format (preferred-centric):
    Gene    Codon_Pos    AA    Preferred_Codon    Ref_Codon    Genomic_Positions    Codon_Variants    Frequencies    Preferred_Freq    Non_Preferred_Freq
    01G000100    2    I    ATT    ATT    14530,14531,14532    ATT:186;ATC:1    ATT:0.995;ATC:0.005    0.995    0.005

Author: GitHub Copilot
Date: 2025-11-06
"""

import sys
import csv
from pathlib import Path
from collections import defaultdict


def load_preferred_codons(preferred_file):
    """
    Load preferred codons from file.
    
    Returns:
        dict: {amino_acid: preferred_codon}
    """
    # Genetic code for amino acid lookup
    genetic_code = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    
    aa_to_preferred = {}
    
    with open(preferred_file, 'r') as f:
        for line in f:
            codon = line.strip().upper()
            if codon in genetic_code:
                aa = genetic_code[codon]
                aa_to_preferred[aa] = codon
    
    return aa_to_preferred


def parse_codon_variants(variant_str):
    """
    Parse codon variants string.
    
    Args:
        variant_str: "ATT:186;ATC:1" or "ATT:186"
    
    Returns:
        dict: {codon: count}
    """
    if not variant_str or variant_str == "":
        return {}
    
    variants = {}
    for part in variant_str.split(';'):
        codon, count = part.split(':')
        variants[codon] = int(count)
    
    return variants


def format_variants(variant_dict, preferred_codon):
    """
    Format variants with preferred codon first.
    
    Args:
        variant_dict: {codon: count}
        preferred_codon: The preferred codon for this AA
    
    Returns:
        tuple: (variants_str, frequencies_str)
    """
    # Ensure preferred codon is included (even if count=0)
    if preferred_codon not in variant_dict:
        variant_dict[preferred_codon] = 0
    
    total = sum(variant_dict.values())
    
    # Sort: preferred first, then by count (descending)
    sorted_codons = sorted(
        variant_dict.keys(),
        key=lambda c: (c != preferred_codon, -variant_dict[c])
    )
    
    variants_parts = []
    freq_parts = []
    
    for codon in sorted_codons:
        count = variant_dict[codon]
        freq = count / total if total > 0 else 0.0
        variants_parts.append(f"{codon}:{count}")
        freq_parts.append(f"{codon}:{freq:.3f}")
    
    return ';'.join(variants_parts), ';'.join(freq_parts)


def calculate_preferred_freq(variant_dict, preferred_codon):
    """
    Calculate frequency of preferred vs non-preferred codons.
    
    Returns:
        tuple: (preferred_freq, non_preferred_freq)
    """
    total = sum(variant_dict.values())
    if total == 0:
        return 0.0, 0.0
    
    preferred_count = variant_dict.get(preferred_codon, 0)
    non_preferred_count = total - preferred_count
    
    return preferred_count / total, non_preferred_count / total


def reformat_to_preferred_centric(input_file, preferred_codons_file, output_file):
    """
    Reformat codon frequency file to be preferred-codon-centric.
    """
    # Load preferred codons
    print("Loading preferred codons...")
    aa_to_preferred = load_preferred_codons(preferred_codons_file)
    print(f"  Loaded {len(aa_to_preferred)} preferred codons")
    
    # Process input file
    print(f"Processing {input_file}...")
    
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')
        
        # Define output columns
        output_fields = [
            'Gene',
            'Codon_Pos',
            'AA',
            'Preferred_Codon',
            'Ref_Codon',
            'Genomic_Positions',
            'Codon_Variants',
            'Frequencies',
            'Preferred_Freq',
            'Non_Preferred_Freq'
        ]
        
        writer = csv.DictWriter(outfile, fieldnames=output_fields, delimiter='\t')
        writer.writeheader()
        
        rows_processed = 0
        
        for row in reader:
            aa = row['AA']
            
            # Get preferred codon for this amino acid
            if aa not in aa_to_preferred:
                # No preferred codon defined (e.g., Met, Trp, stop codons)
                # Use reference codon as "preferred"
                preferred_codon = row['Ref_Codon']
            else:
                preferred_codon = aa_to_preferred[aa]
            
            # Parse variants
            variant_dict = parse_codon_variants(row['Codon_Variants'])
            
            # Reformat with preferred codon first
            variants_str, freq_str = format_variants(variant_dict, preferred_codon)
            
            # Calculate preferred vs non-preferred frequencies
            pref_freq, non_pref_freq = calculate_preferred_freq(variant_dict, preferred_codon)
            
            # Write output row
            writer.writerow({
                'Gene': row['Gene'],
                'Codon_Pos': row['Codon_Pos'],
                'AA': aa,
                'Preferred_Codon': preferred_codon,
                'Ref_Codon': row['Ref_Codon'],
                'Genomic_Positions': row['Genomic_Positions'],
                'Codon_Variants': variants_str,
                'Frequencies': freq_str,
                'Preferred_Freq': f"{pref_freq:.3f}",
                'Non_Preferred_Freq': f"{non_pref_freq:.3f}"
            })
            
            rows_processed += 1
    
    print(f"✓ Processed {rows_processed} codon positions")
    print(f"✓ Output written to: {output_file}")


def main():
    if len(sys.argv) != 4:
        print("Usage: python reformat_codon_freq_by_preference.py <input_codon_freq.txt> <preferred_codons.txt> <output.txt>")
        print("\nExample:")
        print("  python reformat_codon_freq_by_preference.py \\")
        print("      Chr_01.codon_frequencies.txt \\")
        print("      preferred_codons.txt \\")
        print("      Chr_01.codon_frequencies_preferred_centric.txt")
        sys.exit(1)
    
    input_file = sys.argv[1]
    preferred_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Check files exist
    if not Path(input_file).exists():
        print(f"ERROR: Input file not found: {input_file}")
        sys.exit(1)
    
    if not Path(preferred_file).exists():
        print(f"ERROR: Preferred codons file not found: {preferred_file}")
        sys.exit(1)
    
    reformat_to_preferred_centric(input_file, preferred_file, output_file)
    
    print("\nDone!")


if __name__ == "__main__":
    main()
