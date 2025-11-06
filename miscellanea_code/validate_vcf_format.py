#!/usr/bin/env python3
"""
Quick VCF format validator.

Checks:
- Number of columns
- Sample count
- Basic format validation

Usage:
    python validate_vcf_format.py <vcf_file> <sample_mapping_file>

Example:
    python validate_vcf_format.py Included_snp.sites.txt Lines_in_sequence.txt
"""

import sys

def validate_vcf(vcf_file, sample_file):
    print("=" * 70)
    print("VCF Format Validator")
    print("=" * 70)
    
    # Load sample mapping
    print(f"\n📄 Reading sample mapping: {sample_file}")
    sample_map = {}
    with open(sample_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                col_num = int(parts[0])
                sample_id = parts[1]
                sample_map[col_num] = sample_id
    
    print(f"   Found {len(sample_map)} samples")
    print(f"   Column range: {min(sample_map.keys())} to {max(sample_map.keys())}")
    
    # Check VCF
    print(f"\n📄 Reading VCF: {vcf_file}")
    print("   Checking first 10 lines...")
    
    with open(vcf_file, 'r') as f:
        line_count = 0
        for line in f:
            line_count += 1
            
            # Skip headers
            if line.startswith('#'):
                print(f"   Line {line_count}: Header (skipped)")
                continue
            
            # Parse columns
            cols = line.strip().split('\t')
            n_cols = len(cols)
            
            if line_count <= 10:
                # Show first few lines
                chrom = cols[0] if len(cols) > 0 else "?"
                pos = cols[1] if len(cols) > 1 else "?"
                ref = cols[3] if len(cols) > 3 else "?"
                alt = cols[4] if len(cols) > 4 else "?"
                
                print(f"   Line {line_count}: {n_cols} columns | {chrom}:{pos} {ref}→{alt}")
            
            if line_count == 10:
                break
    
    # Analysis
    print(f"\n" + "=" * 70)
    print("ANALYSIS")
    print("=" * 70)
    
    n_standard_cols = 9  # CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT
    n_genotype_cols = n_cols - n_standard_cols
    
    print(f"\n✓ Total columns: {n_cols}")
    print(f"  - Standard VCF columns: {n_standard_cols}")
    print(f"  - Genotype columns: {n_genotype_cols}")
    
    # Compare with sample map
    if n_genotype_cols == len(sample_map):
        print(f"\n✅ PERFECT MATCH!")
        print(f"   VCF has {n_genotype_cols} genotype columns")
        print(f"   Sample file has {len(sample_map)} samples")
        print(f"   All samples accounted for!")
    elif n_genotype_cols > len(sample_map):
        print(f"\n⚠️  WARNING: VCF has MORE columns than samples in mapping file")
        print(f"   VCF genotype columns: {n_genotype_cols}")
        print(f"   Samples in mapping: {len(sample_map)}")
        print(f"   Missing {n_genotype_cols - len(sample_map)} sample IDs")
    else:
        print(f"\n⚠️  WARNING: Sample file has MORE entries than VCF columns")
        print(f"   VCF genotype columns: {n_genotype_cols}")
        print(f"   Samples in mapping: {len(sample_map)}")
        print(f"   Extra {len(sample_map) - n_genotype_cols} sample IDs in file")
    
    # Check column numbering
    expected_first_col = n_standard_cols + 1  # Column 10 for first sample
    actual_first_col = min(sample_map.keys())
    
    print(f"\n✓ Column numbering:")
    print(f"  - First genotype column should be: {expected_first_col}")
    print(f"  - Sample file starts at column: {actual_first_col}")
    
    if expected_first_col == actual_first_col:
        print(f"  ✅ Column numbering matches!")
    else:
        print(f"  ⚠️  Column numbering mismatch!")
        print(f"     VCF columns 1-9 are standard fields")
        print(f"     Genotypes start at column {expected_first_col}")
        print(f"     But your sample file starts at column {actual_first_col}")
    
    # Script compatibility check
    print(f"\n" + "=" * 70)
    print("SCRIPT COMPATIBILITY")
    print("=" * 70)
    
    print(f"\nThe analysis scripts assume:")
    print(f"  - VCF has {n_genotype_cols} samples")
    print(f"  - Genotypes in columns {n_standard_cols + 1} to {n_cols}")
    print(f"  - Watterson's theta calculation uses n={n_genotype_cols}")
    
    if n_genotype_cols == 187:
        print(f"\n✅ Perfect! Your VCF has exactly 187 samples as expected.")
        print(f"   Scripts will work without modification.")
    else:
        print(f"\n⚠️  Your VCF has {n_genotype_cols} samples, not 187.")
        print(f"   You should update the hardcoded value in:")
        print(f"   - calculate_pi.py (line ~263)")
        print(f"   - calculate_pi_optimized.py (line ~249)")
        print(f"   - calculate_freq_preferred_codon.py (line ~327)")
        print(f"\n   Change: n_samples = 187")
        print(f"   To:     n_samples = {n_genotype_cols}")
    
    print()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python validate_vcf_format.py <vcf_file> <sample_mapping_file>")
        print("\nExample:")
        print("  python validate_vcf_format.py Included_snp.sites.txt Lines_in_sequence.txt")
        sys.exit(1)
    
    vcf_file = sys.argv[1]
    sample_file = sys.argv[2]
    
    validate_vcf(vcf_file, sample_file)
