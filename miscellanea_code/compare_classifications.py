#!/usr/bin/env python3
"""
Compare classification systems between proc2.py and TRUE degeneracy approach.

This script compares the site counts and polymorphism metrics between the two systems
to understand how sites are redistributed.

Usage:
    python compare_classifications.py <old_file> <new_file> <gene_id>

Example:
    python compare_classifications.py Chr_01.bygene Chr_01.bygene.pi.txt 01G000100
"""

import sys
import pandas as pd

def parse_old_format(file_path, gene_id):
    """Parse proc2.py output format."""
    with open(file_path, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 14:
                continue
            
            # Chr_01	01G000100	Sites_1	Poly_1	Pi_1	Sites_2	Poly_2	Pi_2	Sites_3_Not_4f	Poly_3_Not_4f	Pi_3_Not_4f	Sites_3_4f	Poly_3_4f	Pi_3_4f
            chrom, gene = cols[0], cols[1]
            
            if gene == gene_id:
                return {
                    'Chr': chrom,
                    'Gene': gene,
                    'Sites_1': int(cols[2]),
                    'Poly_1': int(cols[3]),
                    'Pi_1': float(cols[4]),
                    'Sites_2': int(cols[5]),
                    'Poly_2': int(cols[6]),
                    'Pi_2': float(cols[7]),
                    'Sites_3_Not_4f': int(cols[8]),
                    'Poly_3_Not_4f': int(cols[9]),
                    'Pi_3_Not_4f': float(cols[10]),
                    'Sites_3_fourfold': int(cols[11]),
                    'Poly_3_fourfold': int(cols[12]),
                    'Pi_3_fourfold': float(cols[13])
                }
    return None

def parse_new_format(file_path, gene_id):
    """Parse TRUE degeneracy output format."""
    with open(file_path, 'r') as f:
        header = f.readline()
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 18:
                continue
            
            chrom, gene = cols[0], cols[1]
            
            if gene == gene_id:
                return {
                    'Chr': chrom,
                    'Gene': gene,
                    'Sites_0fold': int(cols[2]),
                    'Poly_0fold': int(cols[3]),
                    'Pi_sum_0fold': float(cols[4]),
                    'Pi_mean_0fold': float(cols[5]),
                    'Sites_2fold': int(cols[6]),
                    'Poly_2fold': int(cols[7]),
                    'Pi_sum_2fold': float(cols[8]),
                    'Pi_mean_2fold': float(cols[9]),
                    'Sites_3fold': int(cols[10]),
                    'Poly_3fold': int(cols[11]),
                    'Pi_sum_3fold': float(cols[12]),
                    'Pi_mean_3fold': float(cols[13]),
                    'Sites_4fold': int(cols[14]),
                    'Poly_4fold': int(cols[15]),
                    'Pi_sum_4fold': float(cols[16]),
                    'Pi_mean_4fold': float(cols[17])
                }
    return None

def compare_gene(old_data, new_data):
    """Compare metrics between old and new classification."""
    print("=" * 80)
    print(f"Comparison for Gene: {old_data['Gene']}")
    print("=" * 80)
    
    print("\nOLD CLASSIFICATION (proc2.py - simplified):")
    print("-" * 80)
    print(f"  1-fold (non-synonymous):")
    print(f"    Sites: {old_data['Sites_1']:6d}  Polymorphic: {old_data['Poly_1']:4d}  π: {old_data['Pi_1']:10.6f}")
    print(f"  2-fold:")
    print(f"    Sites: {old_data['Sites_2']:6d}  Polymorphic: {old_data['Poly_2']:4d}  π: {old_data['Pi_2']:10.6f}")
    print(f"  3-fold (NOT 4-fold):")
    print(f"    Sites: {old_data['Sites_3_Not_4f']:6d}  Polymorphic: {old_data['Poly_3_Not_4f']:4d}  π: {old_data['Pi_3_Not_4f']:10.6f}")
    print(f"  3-fold (4-fold sites):")
    print(f"    Sites: {old_data['Sites_3_fourfold']:6d}  Polymorphic: {old_data['Poly_3_fourfold']:4d}  π: {old_data['Pi_3_fourfold']:10.6f}")
    
    old_total = old_data['Sites_1'] + old_data['Sites_2'] + old_data['Sites_3_Not_4f'] + old_data['Sites_3_fourfold']
    old_total_poly = old_data['Poly_1'] + old_data['Poly_2'] + old_data['Poly_3_Not_4f'] + old_data['Poly_3_fourfold']
    old_total_pi = old_data['Pi_1'] + old_data['Pi_2'] + old_data['Pi_3_Not_4f'] + old_data['Pi_3_fourfold']
    
    print(f"\n  TOTAL:")
    print(f"    Sites: {old_total:6d}  Polymorphic: {old_total_poly:4d}  π_sum: {old_total_pi:10.6f}")
    
    print("\n" + "=" * 80)
    print("\nNEW CLASSIFICATION (TRUE degeneracy - tests all 4 nucleotides):")
    print("-" * 80)
    print(f"  0-fold (non-degenerate):")
    print(f"    Sites: {new_data['Sites_0fold']:6d}  Polymorphic: {new_data['Poly_0fold']:4d}  π_sum: {new_data['Pi_sum_0fold']:10.6f}")
    print(f"  2-fold:")
    print(f"    Sites: {new_data['Sites_2fold']:6d}  Polymorphic: {new_data['Poly_2fold']:4d}  π_sum: {new_data['Pi_sum_2fold']:10.6f}")
    print(f"  3-fold:")
    print(f"    Sites: {new_data['Sites_3fold']:6d}  Polymorphic: {new_data['Poly_3fold']:4d}  π_sum: {new_data['Pi_sum_3fold']:10.6f}")
    print(f"  4-fold:")
    print(f"    Sites: {new_data['Sites_4fold']:6d}  Polymorphic: {new_data['Poly_4fold']:4d}  π_sum: {new_data['Pi_sum_4fold']:10.6f}")
    
    new_total = new_data['Sites_0fold'] + new_data['Sites_2fold'] + new_data['Sites_3fold'] + new_data['Sites_4fold']
    new_total_poly = new_data['Poly_0fold'] + new_data['Poly_2fold'] + new_data['Poly_3fold'] + new_data['Poly_4fold']
    new_total_pi = new_data['Pi_sum_0fold'] + new_data['Pi_sum_2fold'] + new_data['Pi_sum_3fold'] + new_data['Pi_sum_4fold']
    
    print(f"\n  TOTAL:")
    print(f"    Sites: {new_total:6d}  Polymorphic: {new_total_poly:4d}  π_sum: {new_total_pi:10.6f}")
    
    print("\n" + "=" * 80)
    print("\nCOMPARISON:")
    print("-" * 80)
    
    # Total sites comparison
    print(f"\nTotal sites:")
    print(f"  Old: {old_total:6d}")
    print(f"  New: {new_total:6d}")
    print(f"  Difference: {new_total - old_total:6d} ({(new_total - old_total) / old_total * 100:+.1f}%)")
    
    # Check if totals match
    if old_total == new_total:
        print("  ✓ Total sites MATCH - same positions analyzed")
    else:
        print("  ⚠ Total sites DIFFER - different positions analyzed")
    
    # Polymorphic sites comparison
    print(f"\nTotal polymorphic sites:")
    print(f"  Old: {old_total_poly:6d}")
    print(f"  New: {new_total_poly:6d}")
    print(f"  Difference: {new_total_poly - old_total_poly:6d}")
    
    # Pi sum comparison
    print(f"\nTotal π_sum:")
    print(f"  Old: {old_total_pi:10.6f}")
    print(f"  New: {new_total_pi:10.6f}")
    print(f"  Difference: {new_total_pi - old_total_pi:10.6f} ({(new_total_pi - old_total_pi) / old_total_pi * 100:+.2f}%)")
    
    print("\n" + "=" * 80)
    print("\nSITE REDISTRIBUTION ANALYSIS:")
    print("-" * 80)
    
    # Map corresponding categories
    print(f"\nNon-synonymous sites:")
    print(f"  Old '1-fold':     {old_data['Sites_1']:6d} sites")
    print(f"  New '0-fold':     {new_data['Sites_0fold']:6d} sites")
    print(f"  Redistribution:   {new_data['Sites_0fold'] - old_data['Sites_1']:+6d} sites")
    
    print(f"\n2-fold sites:")
    print(f"  Old '2-fold':     {old_data['Sites_2']:6d} sites")
    print(f"  New '2-fold':     {new_data['Sites_2fold']:6d} sites")
    print(f"  Redistribution:   {new_data['Sites_2fold'] - old_data['Sites_2']:+6d} sites")
    
    print(f"\n3-fold sites:")
    print(f"  Old '3_Not_4f':   {old_data['Sites_3_Not_4f']:6d} sites")
    print(f"  New '3-fold':     {new_data['Sites_3fold']:6d} sites")
    print(f"  Redistribution:   {new_data['Sites_3fold'] - old_data['Sites_3_Not_4f']:+6d} sites")
    
    print(f"\n4-fold sites:")
    print(f"  Old '3_fourfold': {old_data['Sites_3_fourfold']:6d} sites")
    print(f"  New '4-fold':     {new_data['Sites_4fold']:6d} sites")
    print(f"  Redistribution:   {new_data['Sites_4fold'] - old_data['Sites_3_fourfold']:+6d} sites")
    
    print("\n" + "=" * 80)
    print("\nKEY INSIGHTS:")
    print("-" * 80)
    
    # Determine redistribution pattern
    if new_data['Sites_0fold'] < old_data['Sites_1']:
        diff = old_data['Sites_1'] - new_data['Sites_0fold']
        print(f"• {diff} sites moved FROM 0-fold TO higher degeneracy classes")
        print(f"  → TRUE degeneracy reveals these sites have synonymous options")
    
    if new_data['Sites_4fold'] > old_data['Sites_3_fourfold']:
        diff = new_data['Sites_4fold'] - old_data['Sites_3_fourfold']
        print(f"• {diff} sites moved TO 4-fold class")
        print(f"  → More sites correctly identified as fully degenerate")
    
    if abs(new_total_pi - old_total_pi) / old_total_pi < 0.01:
        print(f"• Total π is NEARLY IDENTICAL (< 1% difference)")
        print(f"  → Reclassification doesn't change overall diversity estimate")
    
    print("\n" + "=" * 80)

def main():
    if len(sys.argv) != 4:
        print("Usage: python compare_classifications.py <old_bygene_file> <new_pi_file> <gene_id>")
        print("Example: python compare_classifications.py Chr_01.bygene Chr_01.bygene.pi.txt 01G000100")
        sys.exit(1)
    
    old_file = sys.argv[1]
    new_file = sys.argv[2]
    gene_id = sys.argv[3]
    
    print(f"\nLoading data for gene {gene_id}...")
    
    old_data = parse_old_format(old_file, gene_id)
    if old_data is None:
        print(f"ERROR: Gene {gene_id} not found in {old_file}")
        sys.exit(1)
    
    new_data = parse_new_format(new_file, gene_id)
    if new_data is None:
        print(f"ERROR: Gene {gene_id} not found in {new_file}")
        sys.exit(1)
    
    compare_gene(old_data, new_data)
    
    print("\nDone!")

if __name__ == "__main__":
    main()
