#!/usr/bin/env python3
"""
Build Site Frequency Spectrum (SFS) for preferred codons by amino acid family.

This script reads the output from reformat_codon_freq_by_preference.py which already
contains the frequency of the preferred codon at each position across all samples.

For each polymorphic synonymous site, we extract the preferred codon frequency
and build an SFS grouped by amino acid and degeneracy family.

Prediction under selection for preferred codons:
  - Mean preferred codon frequency should be > 0.5
  - SFS should be skewed toward high frequencies

Under neutrality:
  - Mean preferred codon frequency ≈ 0.5
  - SFS should be symmetric

Input:
    - all_chromosomes.codon_frequencies_preferred.txt (from reformat_codon_freq_by_preference.py)

Output:
    - SFS summary statistics per amino acid
    - Raw preferred frequencies for custom analysis in R
    - Binned SFS for visualization

Author: Luis Javier Madrigal-Roca & GitHub Copilot
Date: 2025-12-05
"""

import sys
import csv
from collections import defaultdict

# Amino acid families by degeneracy
AA_FAMILIES = {
    '2-fold': ['F', 'Y', 'H', 'Q', 'N', 'K', 'D', 'E', 'C'],
    '3-fold': ['I'],
    '4-fold': ['A', 'G', 'P', 'T', 'V'],
    '6-fold': ['L', 'S', 'R'],
    '1-fold': ['M', 'W']  # No synonymous changes possible
}

# Reverse mapping: AA -> family
AA_TO_FAMILY = {}
for family, aas in AA_FAMILIES.items():
    for aa in aas:
        AA_TO_FAMILY[aa] = family


def process_codon_frequencies(input_file):
    """
    Read preferred codon frequencies and build SFS per amino acid.
    
    Returns:
        sfs_data: dict of amino_acid -> {
            'preferred_freqs': [list of preferred codon frequencies at polymorphic sites],
            'n_invariant_pref': count of sites fixed for preferred,
            'n_invariant_nonpref': count of sites fixed for non-preferred,
            'n_total': total sites
        }
    """
    sfs_data = defaultdict(lambda: {
        'preferred_freqs': [],
        'n_invariant_pref': 0,
        'n_invariant_nonpref': 0,
        'n_total': 0
    })
    
    print(f"Processing {input_file}...")
    
    line_count = 0
    polymorphic_count = 0
    
    with open(input_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        for row in reader:
            line_count += 1
            
            if line_count % 1000000 == 0:
                print(f"  Processed {line_count:,} codon positions...", file=sys.stderr)
            
            aa = row['AA']
            
            # Skip stop codons and single-codon amino acids
            if aa in ['*', 'M', 'W']:
                continue
            
            sfs_data[aa]['n_total'] += 1
            
            # Parse preferred frequency
            try:
                pref_freq = float(row['Preferred_Freq'])
            except (ValueError, KeyError):
                continue
            
            # Categorize by polymorphism status
            # Using small epsilon for floating point comparison
            if pref_freq > 0.999:
                # Fixed for preferred codon
                sfs_data[aa]['n_invariant_pref'] += 1
            elif pref_freq < 0.001:
                # Fixed for non-preferred codon
                sfs_data[aa]['n_invariant_nonpref'] += 1
            else:
                # Polymorphic - add to SFS
                sfs_data[aa]['preferred_freqs'].append(pref_freq)
                polymorphic_count += 1
    
    print(f"  Total codon positions: {line_count:,}")
    print(f"  Polymorphic positions: {polymorphic_count:,}")
    
    return dict(sfs_data)


def calculate_sfs_bins(freq_list, n_bins=10):
    """
    Bin frequencies into SFS categories.
    
    Returns:
        counts: list of counts per bin
        bin_edges: list of bin boundaries
    """
    bin_edges = [i / n_bins for i in range(n_bins + 1)]
    counts = [0] * n_bins
    
    for freq in freq_list:
        bin_idx = min(int(freq * n_bins), n_bins - 1)
        counts[bin_idx] += 1
    
    return counts, bin_edges


def calculate_stats(freq_list):
    """Calculate summary statistics for a frequency list."""
    if not freq_list:
        return {
            'n': 0,
            'mean': None,
            'median': None,
            'std': None,
            'skewness': None
        }
    
    n = len(freq_list)
    mean = sum(freq_list) / n
    
    # Sort for median
    sorted_freqs = sorted(freq_list)
    if n % 2 == 0:
        median = (sorted_freqs[n//2 - 1] + sorted_freqs[n//2]) / 2
    else:
        median = sorted_freqs[n//2]
    
    # Standard deviation
    variance = sum((x - mean) ** 2 for x in freq_list) / n
    std = variance ** 0.5
    
    # Skewness (Fisher's)
    if std > 0:
        skewness = sum(((x - mean) / std) ** 3 for x in freq_list) / n
    else:
        skewness = 0
    
    return {
        'n': n,
        'mean': mean,
        'median': median,
        'std': std,
        'skewness': skewness
    }


def write_output(sfs_data, output_prefix):
    """Write SFS results to output files."""
    
    # 1. Per-amino-acid summary
    summary_file = f"{output_prefix}.aa_summary.txt"
    print(f"Writing amino acid summary to {summary_file}...")
    
    with open(summary_file, 'w') as out:
        out.write("Amino_Acid\tFamily\t")
        out.write("N_polymorphic\tN_fixed_pref\tN_fixed_nonpref\tN_total\t")
        out.write("Mean_pref_freq\tMedian_pref_freq\tStd_pref_freq\tSkewness\t")
        out.write("Prop_fixed_pref\n")
        
        for aa in sorted(sfs_data.keys()):
            data = sfs_data[aa]
            family = AA_TO_FAMILY.get(aa, 'unknown')
            stats = calculate_stats(data['preferred_freqs'])
            
            n_poly = stats['n']
            n_fixed_pref = data['n_invariant_pref']
            n_fixed_nonpref = data['n_invariant_nonpref']
            n_total = data['n_total']
            
            # Proportion of fixed sites that are preferred
            n_fixed_total = n_fixed_pref + n_fixed_nonpref
            prop_fixed_pref = n_fixed_pref / n_fixed_total if n_fixed_total > 0 else 0
            
            out.write(f"{aa}\t{family}\t")
            out.write(f"{n_poly}\t{n_fixed_pref}\t{n_fixed_nonpref}\t{n_total}\t")
            
            if stats['mean'] is not None:
                out.write(f"{stats['mean']:.4f}\t{stats['median']:.4f}\t{stats['std']:.4f}\t{stats['skewness']:.4f}\t")
            else:
                out.write("NA\tNA\tNA\tNA\t")
            
            out.write(f"{prop_fixed_pref:.4f}\n")
    
    # 2. Binned SFS per amino acid
    sfs_file = f"{output_prefix}.sfs_binned.txt"
    print(f"Writing binned SFS to {sfs_file}...")
    
    n_bins = 10
    
    with open(sfs_file, 'w') as out:
        # Header with bin ranges
        out.write("Amino_Acid\tFamily\t")
        bin_headers = [f"Bin_{i+1}_({i/n_bins:.1f}-{(i+1)/n_bins:.1f})" for i in range(n_bins)]
        out.write("\t".join(bin_headers) + "\n")
        
        for aa in sorted(sfs_data.keys()):
            data = sfs_data[aa]
            family = AA_TO_FAMILY.get(aa, 'unknown')
            
            counts, _ = calculate_sfs_bins(data['preferred_freqs'], n_bins)
            
            out.write(f"{aa}\t{family}\t")
            out.write("\t".join(map(str, counts)) + "\n")
    
    # 3. Raw frequencies for R analysis
    raw_file = f"{output_prefix}.sfs_raw.txt"
    print(f"Writing raw frequencies to {raw_file}...")
    
    with open(raw_file, 'w') as out:
        out.write("Amino_Acid\tFamily\tPreferred_Freq\n")
        
        for aa in sorted(sfs_data.keys()):
            data = sfs_data[aa]
            family = AA_TO_FAMILY.get(aa, 'unknown')
            
            for freq in data['preferred_freqs']:
                out.write(f"{aa}\t{family}\t{freq:.4f}\n")
    
    # 4. Summary by degeneracy family
    family_file = f"{output_prefix}.family_summary.txt"
    print(f"Writing family summary to {family_file}...")
    
    # Aggregate by family
    family_data = defaultdict(lambda: {
        'preferred_freqs': [],
        'n_invariant_pref': 0,
        'n_invariant_nonpref': 0,
        'n_total': 0
    })
    
    for aa, data in sfs_data.items():
        family = AA_TO_FAMILY.get(aa, 'unknown')
        family_data[family]['preferred_freqs'].extend(data['preferred_freqs'])
        family_data[family]['n_invariant_pref'] += data['n_invariant_pref']
        family_data[family]['n_invariant_nonpref'] += data['n_invariant_nonpref']
        family_data[family]['n_total'] += data['n_total']
    
    with open(family_file, 'w') as out:
        out.write("Family\t")
        out.write("N_polymorphic\tN_fixed_pref\tN_fixed_nonpref\tN_total\t")
        out.write("Mean_pref_freq\tMedian_pref_freq\tStd_pref_freq\tSkewness\t")
        out.write("Prop_fixed_pref\n")
        
        for family in ['2-fold', '3-fold', '4-fold', '6-fold']:
            data = family_data[family]
            stats = calculate_stats(data['preferred_freqs'])
            
            n_poly = stats['n']
            n_fixed_pref = data['n_invariant_pref']
            n_fixed_nonpref = data['n_invariant_nonpref']
            n_total = data['n_total']
            
            n_fixed_total = n_fixed_pref + n_fixed_nonpref
            prop_fixed_pref = n_fixed_pref / n_fixed_total if n_fixed_total > 0 else 0
            
            out.write(f"{family}\t")
            out.write(f"{n_poly}\t{n_fixed_pref}\t{n_fixed_nonpref}\t{n_total}\t")
            
            if stats['mean'] is not None:
                out.write(f"{stats['mean']:.4f}\t{stats['median']:.4f}\t{stats['std']:.4f}\t{stats['skewness']:.4f}\t")
            else:
                out.write("NA\tNA\tNA\tNA\t")
            
            out.write(f"{prop_fixed_pref:.4f}\n")


def print_summary(sfs_data):
    """Print summary statistics to console."""
    
    print("\n" + "=" * 70)
    print("PREFERRED CODON FREQUENCY - SFS SUMMARY")
    print("=" * 70)
    
    # Aggregate all frequencies
    all_freqs = []
    total_fixed_pref = 0
    total_fixed_nonpref = 0
    
    for aa, data in sfs_data.items():
        all_freqs.extend(data['preferred_freqs'])
        total_fixed_pref += data['n_invariant_pref']
        total_fixed_nonpref += data['n_invariant_nonpref']
    
    stats = calculate_stats(all_freqs)
    
    print(f"\nPolymorphic sites: {stats['n']:,}")
    print(f"Fixed for preferred: {total_fixed_pref:,}")
    print(f"Fixed for non-preferred: {total_fixed_nonpref:,}")
    
    if stats['mean'] is not None:
        print(f"\nPreferred codon frequency at polymorphic sites:")
        print(f"  Mean:     {stats['mean']:.4f}")
        print(f"  Median:   {stats['median']:.4f}")
        print(f"  Std Dev:  {stats['std']:.4f}")
        print(f"  Skewness: {stats['skewness']:.4f}")
        
        print(f"\n--- Interpretation ---")
        
        # Test against neutral expectation (0.5)
        if stats['mean'] > 0.5:
            diff = stats['mean'] - 0.5
            print(f"Mean preferred frequency is ABOVE 0.5 (+{diff:.4f})")
            print("→ Consistent with SELECTION FOR preferred codons")
        elif stats['mean'] < 0.5:
            diff = 0.5 - stats['mean']
            print(f"Mean preferred frequency is BELOW 0.5 (-{diff:.4f})")
            print("→ Unexpected: Selection AGAINST preferred codons?")
        else:
            print("Mean preferred frequency equals 0.5")
            print("→ Consistent with NEUTRAL evolution")
        
        if stats['skewness'] > 0:
            print(f"\nPositive skewness ({stats['skewness']:.4f})")
            print("→ Distribution shifted toward high frequencies (selection signature)")
        elif stats['skewness'] < 0:
            print(f"\nNegative skewness ({stats['skewness']:.4f})")
            print("→ Distribution shifted toward low frequencies")
        
        # Proportion of fixed sites that are preferred
        total_fixed = total_fixed_pref + total_fixed_nonpref
        if total_fixed > 0:
            prop_pref = total_fixed_pref / total_fixed
            print(f"\nProportion of fixed sites with preferred codon: {prop_pref:.4f}")
            if prop_pref > 0.5:
                print("→ More sites fixed for preferred than non-preferred (selection signature)")
    
    print("\n" + "=" * 70)


def main():
    if len(sys.argv) != 3:
        print("Usage: python build_sfs_from_preferred_freqs.py <codon_frequencies_preferred.txt> <output_prefix>")
        print()
        print("Example:")
        print("  python build_sfs_from_preferred_freqs.py all_chromosomes.codon_frequencies_preferred.txt sfs_results")
        print()
        print("Input:")
        print("  Output from reformat_codon_freq_by_preference.py with columns:")
        print("  Gene, Codon_Pos, AA, Preferred_Codon, Ref_Codon, Genomic_Positions,")
        print("  Codon_Variants, Frequencies, Preferred_Freq, Non_Preferred_Freq")
        print()
        print("Output files:")
        print("  <prefix>.aa_summary.txt     - Per-amino-acid summary statistics")
        print("  <prefix>.sfs_binned.txt     - Binned SFS per amino acid")
        print("  <prefix>.sfs_raw.txt        - Raw frequencies for R analysis")
        print("  <prefix>.family_summary.txt - Summary by degeneracy family")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_prefix = sys.argv[2]
    
    print("=== Building SFS from Preferred Codon Frequencies ===")
    print()
    
    sfs_data = process_codon_frequencies(input_file)
    
    write_output(sfs_data, output_prefix)
    
    print_summary(sfs_data)
    
    print("\nDone!")


if __name__ == "__main__":
    main()
