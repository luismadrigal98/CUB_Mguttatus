#!/usr/bin/env python3
"""
Analyze preferred codon frequency distribution at polymorphic sites.

This script provides DESCRIPTIVE statistics on the frequency of preferred codons
at polymorphic synonymous sites, grouped by amino acid and degeneracy family.

IMPORTANT LIMITATIONS:
- Without outgroup sequences, cannot distinguish selection from drift
- Neutral expectation depends on degeneracy (NOT 0.5 for all amino acids)
- Cannot infer direction of selection without ancestral state polarization
- Results are descriptive and require additional tests for inference

For rigorous selection inference, you need:
1. Outgroup sequences to polarize mutations
2. McDonald-Kreitman test (polymorphism vs. divergence)
3. Correlation with gene expression levels
4. Comparison to neutral synonymous sites

What this analysis DOES provide:
- Current state of preferred codon usage at variable sites
- Comparison across amino acid families
- Baseline for between-population or temporal comparisons
- Identification of highly variable vs. conserved sites

Input:
    - Codon frequency table with preferred codon annotations
      (from reformat_codon_freq_by_preference.py)

Output:
    - Summary statistics per amino acid and degeneracy family
    - Frequency distributions for visualization
    - Raw data for additional statistical tests

Author: Luis Javier Madrigal-Roca & John K. Kelly
Date: 2025-12-05
"""

import sys
import os
import csv
import time
from collections import defaultdict

# =============================================================================
# Configuration Constants
# =============================================================================

# Thresholds for classifying sites as fixed vs polymorphic
FIXED_THRESHOLD_HIGH = 0.999  # Sites with freq > this are "fixed for preferred"
FIXED_THRESHOLD_LOW = 0.001   # Sites with freq < this are "fixed for non-preferred"

# Number of bins for frequency histogram
N_BINS = 10

# Amino acid families by degeneracy
# M (Met) and W (Trp) are excluded - single codon, no synonymous variation possible
# Stop codons (*) are also excluded
AA_FAMILIES = {
    '2-fold': ['F', 'Y', 'H', 'Q', 'N', 'K', 'D', 'E', 'C'],
    '3-fold': ['I'],
    '4-fold': ['A', 'G', 'P', 'T', 'V'],
    '6-fold': ['L', 'S', 'R']
}

# Amino acids to skip (no synonymous variation or not applicable)
SKIP_AMINO_ACIDS = {'*', 'M', 'W'}

# Number of synonymous codons per amino acid (for neutral expectation calculation)
CODON_COUNTS = {
    'F': 2, 'Y': 2, 'C': 2, 'H': 2, 'Q': 2,
    'N': 2, 'K': 2, 'D': 2, 'E': 2,  # 2-fold
    'I': 3,  # 3-fold
    'A': 4, 'G': 4, 'P': 4, 'T': 4, 'V': 4,  # 4-fold
    'L': 6, 'S': 6, 'R': 6  # 6-fold
}

# Reverse mapping: AA -> family
AA_TO_FAMILY = {}
for family, aas in AA_FAMILIES.items():
    for aa in aas:
        AA_TO_FAMILY[aa] = family


def get_neutral_expectation(aa):
    """
    Calculate expected preferred codon frequency under neutrality
    (equal mutation rates, no selection).
    
    For amino acid with N synonymous codons, if 1 is preferred:
    Expected frequency = 1/N
    
    Note: This assumes:
    - Equal mutation rates between all codons
    - No selection
    - No mutation bias
    - Equilibrium conditions
    
    Args:
        aa: Single-letter amino acid code
        
    Returns:
        float: Expected frequency under neutrality, or None if not applicable
    """
    n_codons = CODON_COUNTS.get(aa)
    if n_codons and n_codons > 1:
        return 1.0 / n_codons
    return None


def process_codon_frequencies(input_file):
    """
    Read preferred codon frequencies and summarize distribution per amino acid.
    
    Args:
        input_file: Path to codon frequency file with Preferred_Freq column
    
    Returns:
        preferred_freq_data: dict of amino_acid -> {
            'preferred_freqs': [list of preferred codon frequencies at polymorphic sites],
            'n_invariant_pref': count of sites fixed for preferred,
            'n_invariant_nonpref': count of sites fixed for non-preferred,
            'n_total': total sites
        }
        parse_errors: int, count of lines with parsing errors
    """
    preferred_freq_data = defaultdict(lambda: {
        'preferred_freqs': [],
        'n_invariant_pref': 0,
        'n_invariant_nonpref': 0,
        'n_total': 0
    })
    
    print(f"Processing {input_file}...", file=sys.stderr)
    
    line_count = 0
    polymorphic_count = 0
    parse_errors = 0
    out_of_range_count = 0
    skipped_aa_count = 0
    
    with open(input_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        # Validate required columns exist
        required_cols = ['AA', 'Preferred_Freq']
        if reader.fieldnames is None:
            print("ERROR: Input file appears to be empty or has no header", file=sys.stderr)
            sys.exit(1)
        
        missing_cols = [col for col in required_cols if col not in reader.fieldnames]
        if missing_cols:
            print(f"ERROR: Missing required columns: {missing_cols}", file=sys.stderr)
            print(f"Found columns: {reader.fieldnames}", file=sys.stderr)
            sys.exit(1)
        
        for row in reader:
            line_count += 1
            
            if line_count % 1000000 == 0:
                print(f"  Processed {line_count:,} codon positions...", file=sys.stderr)
            
            aa = row['AA']
            
            # Skip amino acids without synonymous variation
            if aa in SKIP_AMINO_ACIDS or aa not in AA_TO_FAMILY:
                skipped_aa_count += 1
                continue
            
            preferred_freq_data[aa]['n_total'] += 1
            
            # Parse preferred frequency with error handling
            try:
                pref_freq = float(row['Preferred_Freq'])
            except (ValueError, KeyError, TypeError):
                parse_errors += 1
                if parse_errors <= 5:
                    print(f"  Warning: Could not parse frequency at line {line_count + 1}: "
                          f"{row.get('Preferred_Freq', 'MISSING')}", file=sys.stderr)
                continue
            
            # Validate frequency is in valid range [0, 1]
            if pref_freq < 0 or pref_freq > 1:
                out_of_range_count += 1
                if out_of_range_count <= 5:
                    print(f"  Warning: Frequency out of range [0,1] at line {line_count + 1}: {pref_freq}",
                          file=sys.stderr)
                continue
            
            # Categorize by polymorphism status
            if pref_freq > FIXED_THRESHOLD_HIGH:
                # Fixed for preferred codon
                preferred_freq_data[aa]['n_invariant_pref'] += 1
            elif pref_freq < FIXED_THRESHOLD_LOW:
                # Fixed for non-preferred codon
                preferred_freq_data[aa]['n_invariant_nonpref'] += 1
            else:
                # Polymorphic - add to frequency distribution
                preferred_freq_data[aa]['preferred_freqs'].append(pref_freq)
                polymorphic_count += 1
    
    # Report parsing summary
    print(f"  Total codon positions: {line_count:,}", file=sys.stderr)
    print(f"  Polymorphic positions: {polymorphic_count:,}", file=sys.stderr)
    print(f"  Skipped (M/W/*): {skipped_aa_count:,}", file=sys.stderr)
    
    if parse_errors > 0:
        print(f"  WARNING: {parse_errors:,} lines had parsing errors", file=sys.stderr)
    if out_of_range_count > 0:
        print(f"  WARNING: {out_of_range_count:,} lines had frequencies outside [0,1]", file=sys.stderr)
    
    return dict(preferred_freq_data), parse_errors


def bin_preferred_frequencies(freq_list, n_bins=N_BINS):
    """
    Bin frequencies into histogram categories.
    
    Args:
        freq_list: list of frequencies in [0, 1]
        n_bins: number of bins (default from N_BINS constant)
    
    Returns:
        counts: list of counts per bin
        bin_edges: list of bin boundaries
    """
    bin_edges = [i / n_bins for i in range(n_bins + 1)]
    counts = [0] * n_bins
    
    for freq in freq_list:
        # Handle edge case: freq == 1.0 should go to last bin
        bin_idx = min(int(freq * n_bins), n_bins - 1)
        # Ensure bin_idx is not negative (shouldn't happen with valid freqs)
        bin_idx = max(0, bin_idx)
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


def write_output(preferred_freq_data, output_prefix):
    """Write preferred codon frequency analysis results to output files."""
    
    # ==========================================================================
    # 1. Per-amino-acid summary with caveats header
    # ==========================================================================
    summary_file = f"{output_prefix}.aa_summary.txt"
    print(f"Writing amino acid summary to {summary_file}...", file=sys.stderr)
    
    with open(summary_file, 'w') as out:
        # Add header with caveats
        out.write("# Preferred Codon Frequency Analysis\n")
        out.write("# \n")
        out.write("# IMPORTANT CAVEATS:\n")
        out.write("#   - These are DESCRIPTIVE statistics only\n")
        out.write("#   - Without outgroup: cannot distinguish selection from drift\n")
        out.write("#   - Neutral expectation is NOT 0.5 for all amino acids\n")
        out.write("#   - For 4-fold sites with 1 preferred codon: neutral expectation = 0.25\n")
        out.write("#   - Mean frequency depends on: selection, mutation bias, demography, drift\n")
        out.write("# \n")
        out.write("# To test for selection, compare to:\n")
        out.write("#   1. Neutral synonymous sites genome-wide\n")
        out.write("#   2. Gene expression levels (high expression → stronger CUB)\n")
        out.write("#   3. Between-population comparisons\n")
        out.write("# \n")
        
        # Column headers
        out.write("Amino_Acid\tFamily\t")
        out.write("N_polymorphic\tN_fixed_pref\tN_fixed_nonpref\tN_total\t")
        out.write("Mean_pref_freq\tMedian_pref_freq\tStd_pref_freq\tSkewness\t")
        out.write("Neutral_expectation\tDiff_from_neutral\t")
        out.write("Prop_fixed_pref\n")
        
        for aa in sorted(preferred_freq_data.keys()):
            data = preferred_freq_data[aa]
            family = AA_TO_FAMILY.get(aa, 'unknown')
            stats = calculate_stats(data['preferred_freqs'])
            
            n_poly = stats['n']
            n_fixed_pref = data['n_invariant_pref']
            n_fixed_nonpref = data['n_invariant_nonpref']
            n_total = data['n_total']
            
            # Proportion of fixed sites that are preferred
            n_fixed_total = n_fixed_pref + n_fixed_nonpref
            prop_fixed_pref = n_fixed_pref / n_fixed_total if n_fixed_total > 0 else 0
            
            # Neutral expectation
            neutral_exp = get_neutral_expectation(aa)
            
            out.write(f"{aa}\t{family}\t")
            out.write(f"{n_poly}\t{n_fixed_pref}\t{n_fixed_nonpref}\t{n_total}\t")
            
            if stats['mean'] is not None:
                out.write(f"{stats['mean']:.4f}\t{stats['median']:.4f}\t")
                out.write(f"{stats['std']:.4f}\t{stats['skewness']:.4f}\t")
                
                # Neutral expectation and difference
                if neutral_exp is not None:
                    diff_from_neutral = stats['mean'] - neutral_exp
                    out.write(f"{neutral_exp:.4f}\t{diff_from_neutral:+.4f}\t")
                else:
                    out.write("NA\tNA\t")
            else:
                out.write("NA\tNA\tNA\tNA\tNA\tNA\t")
            
            out.write(f"{prop_fixed_pref:.4f}\n")
    
    # ==========================================================================
    # 2. Binned frequency distribution per amino acid
    # ==========================================================================
    freq_dist_file = f"{output_prefix}.freq_distribution.txt"
    print(f"Writing frequency distribution to {freq_dist_file}...", file=sys.stderr)
    
    n_bins = N_BINS
    
    with open(freq_dist_file, 'w') as out:
        # Header with caveats
        out.write("# Preferred Codon Frequency Distribution (Binned)\n")
        out.write("# These are frequencies of preferred codons at POLYMORPHIC sites only\n")
        out.write("# \n")
        
        # Column headers with bin ranges
        out.write("Amino_Acid\tFamily\t")
        bin_headers = [f"Bin_{i+1}_({i/n_bins:.1f}-{(i+1)/n_bins:.1f})" for i in range(n_bins)]
        out.write("\t".join(bin_headers) + "\n")
        
        for aa in sorted(preferred_freq_data.keys()):
            data = preferred_freq_data[aa]
            family = AA_TO_FAMILY.get(aa, 'unknown')
            
            counts, _ = bin_preferred_frequencies(data['preferred_freqs'], n_bins)
            
            out.write(f"{aa}\t{family}\t")
            out.write("\t".join(map(str, counts)) + "\n")
    
    # ==========================================================================
    # 3. Raw frequencies for R analysis
    # ==========================================================================
    raw_file = f"{output_prefix}.raw_frequencies.txt"
    print(f"Writing raw frequencies to {raw_file}...", file=sys.stderr)
    
    with open(raw_file, 'w') as out:
        out.write("# Raw preferred codon frequencies at polymorphic sites\n")
        out.write("# For custom statistical analysis in R\n")
        out.write("# \n")
        out.write("Amino_Acid\tFamily\tPreferred_Freq\n")
        
        for aa in sorted(preferred_freq_data.keys()):
            data = preferred_freq_data[aa]
            family = AA_TO_FAMILY.get(aa, 'unknown')
            
            for freq in data['preferred_freqs']:
                out.write(f"{aa}\t{family}\t{freq:.4f}\n")
    
    # ==========================================================================
    # 4. Summary by degeneracy family
    # ==========================================================================
    family_file = f"{output_prefix}.family_summary.txt"
    print(f"Writing family summary to {family_file}...", file=sys.stderr)
    
    # Aggregate by family
    family_data = defaultdict(lambda: {
        'preferred_freqs': [],
        'n_invariant_pref': 0,
        'n_invariant_nonpref': 0,
        'n_total': 0
    })
    
    for aa, data in preferred_freq_data.items():
        family = AA_TO_FAMILY.get(aa, 'unknown')
        family_data[family]['preferred_freqs'].extend(data['preferred_freqs'])
        family_data[family]['n_invariant_pref'] += data['n_invariant_pref']
        family_data[family]['n_invariant_nonpref'] += data['n_invariant_nonpref']
        family_data[family]['n_total'] += data['n_total']
    
    with open(family_file, 'w') as out:
        # Header with caveats
        out.write("# Preferred Codon Frequency Summary by Degeneracy Family\n")
        out.write("# \n")
        out.write("# Neutral expectation varies by degeneracy:\n")
        out.write("#   2-fold: 0.50 (1 preferred out of 2 codons)\n")
        out.write("#   3-fold: 0.33 (1 preferred out of 3 codons)\n")
        out.write("#   4-fold: 0.25 (1 preferred out of 4 codons)\n")
        out.write("#   6-fold: 0.17 (1 preferred out of 6 codons)\n")
        out.write("# \n")
        
        out.write("Family\t")
        out.write("N_polymorphic\tN_fixed_pref\tN_fixed_nonpref\tN_total\t")
        out.write("Mean_pref_freq\tMedian_pref_freq\tStd_pref_freq\tSkewness\t")
        out.write("Neutral_expectation\tDiff_from_neutral\t")
        out.write("Prop_fixed_pref\n")
        
        # Neutral expectations by family
        family_neutral_exp = {
            '2-fold': 0.50,
            '3-fold': 0.33,
            '4-fold': 0.25,
            '6-fold': 0.17
        }
        
        for family in ['2-fold', '3-fold', '4-fold', '6-fold']:
            data = family_data[family]
            stats = calculate_stats(data['preferred_freqs'])
            
            n_poly = stats['n']
            n_fixed_pref = data['n_invariant_pref']
            n_fixed_nonpref = data['n_invariant_nonpref']
            n_total = data['n_total']
            
            n_fixed_total = n_fixed_pref + n_fixed_nonpref
            prop_fixed_pref = n_fixed_pref / n_fixed_total if n_fixed_total > 0 else 0
            
            neutral_exp = family_neutral_exp[family]
            
            out.write(f"{family}\t")
            out.write(f"{n_poly}\t{n_fixed_pref}\t{n_fixed_nonpref}\t{n_total}\t")
            
            if stats['mean'] is not None:
                diff_from_neutral = stats['mean'] - neutral_exp
                out.write(f"{stats['mean']:.4f}\t{stats['median']:.4f}\t")
                out.write(f"{stats['std']:.4f}\t{stats['skewness']:.4f}\t")
                out.write(f"{neutral_exp:.2f}\t{diff_from_neutral:+.4f}\t")
            else:
                out.write("NA\tNA\tNA\tNA\tNA\tNA\t")
            
            out.write(f"{prop_fixed_pref:.4f}\n")


def print_summary(preferred_freq_data):
    """Print summary statistics with appropriate caveats."""
    
    print("\n" + "=" * 70)
    print("PREFERRED CODON FREQUENCY DISTRIBUTION - DESCRIPTIVE ANALYSIS")
    print("=" * 70)
    
    # Aggregate all frequencies
    all_freqs = []
    total_fixed_pref = 0
    total_fixed_nonpref = 0
    
    for aa, data in preferred_freq_data.items():
        all_freqs.extend(data['preferred_freqs'])
        total_fixed_pref += data['n_invariant_pref']
        total_fixed_nonpref += data['n_invariant_nonpref']
    
    stats = calculate_stats(all_freqs)
    
    print(f"\n--- Summary Statistics (All Amino Acids Combined) ---")
    print(f"Polymorphic sites analyzed: {stats['n']:,}")
    
    if stats['mean'] is not None:
        print(f"Mean preferred codon frequency: {stats['mean']:.4f}")
        print(f"Median: {stats['median']:.4f}")
        print(f"Std Dev: {stats['std']:.4f}")
        print(f"Skewness: {stats['skewness']:.4f}")
        
        print(f"\n--- What These Numbers Mean ---")
        print(f"At polymorphic synonymous sites:")
        print(f"  • Preferred codons have mean frequency of {stats['mean']:.2%}")
        print(f"  • Half of sites have preferred frequency above {stats['median']:.2%}")
        
        if stats['mean'] > 0.5:
            print(f"\nPreferred codons are MORE COMMON than non-preferred at polymorphic sites.")
            print(f"This could indicate:")
            print(f"  • Selection favoring preferred codons")
            print(f"  • Mutation bias toward preferred codons")
            print(f"  • Ancestral state is often the preferred codon")
            print(f"  • Demographic history (bottleneck, population structure)")
        else:
            print(f"\nPreferred codons are LESS COMMON than non-preferred at polymorphic sites.")
            print(f"This could indicate:")
            print(f"  • Recent mutations creating non-preferred variants")
            print(f"  • Weak or absent selection on codon usage")
            print(f"  • Ancestral state is often non-preferred")
        
        # Skewness interpretation
        print(f"\n--- Skewness Interpretation ---")
        if stats['skewness'] > 0.5:
            print(f"Positive skewness ({stats['skewness']:.4f})")
            print("→ Distribution is right-skewed (long tail toward high frequencies)")
            print("  Most sites have lower preferred codon frequency,")
            print("  but some sites have very high preferred codon frequency")
        elif stats['skewness'] < -0.5:
            print(f"Negative skewness ({stats['skewness']:.4f})")
            print("→ Distribution is left-skewed (long tail toward low frequencies)")
            print("  Most sites have higher preferred codon frequency,")
            print("  but some sites have very low preferred codon frequency")
        else:
            print(f"Near-symmetric distribution (skewness: {stats['skewness']:.4f})")
            print("→ Distribution is approximately symmetric around the mean")
        
        print("\nNote: Skewness alone cannot distinguish selection from demography.")
    
    print(f"\n--- Fixed Sites ---")
    print(f"Sites fixed for preferred: {total_fixed_pref:,}")
    print(f"Sites fixed for non-preferred: {total_fixed_nonpref:,}")
    
    total_fixed = total_fixed_pref + total_fixed_nonpref
    if total_fixed > 0:
        prop_pref = total_fixed_pref / total_fixed
        print(f"Proportion fixed as preferred: {prop_pref:.2%}")
        
        if prop_pref > 0.5:
            print(f"More sites are fixed for preferred than non-preferred codons.")
    
    print(f"\n--- IMPORTANT LIMITATIONS ---")
    print(f"✗ These statistics are DESCRIPTIVE only")
    print(f"✗ Cannot distinguish selection from neutral processes without:")
    print(f"    • Outgroup sequences (to polarize ancestral/derived)")
    print(f"    • Neutral baseline (intergenic regions or introns)")
    print(f"    • Expression data (to test CUB-expression correlation)")
    print(f"    • Population structure information")
    print(f"✗ Neutral expectation is NOT 0.5 for most amino acids")
    print(f"    • For 2-fold degenerate sites: expected = 0.50")
    print(f"    • For 3-fold degenerate sites: expected = 0.33")
    print(f"    • For 4-fold degenerate sites: expected = 0.25")
    print(f"    • For 6-fold degenerate sites: expected = 0.17")
    
    print(f"\n--- What's Needed for Proper Inference ---")
    print(f"To test for selection on codon usage, you would need:")
    print(f"  1. Outgroup sequences to polarize ancestral vs. derived states")
    print(f"  2. Comparison to neutral synonymous sites genome-wide")
    print(f"  3. McDonald-Kreitman test (polymorphism vs. divergence)")
    print(f"  4. Correlation with gene expression levels")
    
    print("\n" + "=" * 70)


def main():
    if len(sys.argv) != 3:
        print("Usage: python analyze_preferred_codon_frequencies.py <codon_frequencies_preferred.txt> <output_prefix>")
        print()
        print("Example:")
        print("  python analyze_preferred_codon_frequencies.py all_chromosomes.codon_frequencies_preferred.txt codon_analysis")
        print()
        print("Input:")
        print("  Output from reformat_codon_freq_by_preference.py with columns:")
        print("  Gene, Codon_Pos, AA, Preferred_Codon, Ref_Codon, Genomic_Positions,")
        print("  Codon_Variants, Frequencies, Preferred_Freq, Non_Preferred_Freq")
        print()
        print("Output files:")
        print("  <prefix>.aa_summary.txt       - Per-amino-acid summary statistics")
        print("  <prefix>.freq_distribution.txt - Binned frequency distribution")
        print("  <prefix>.raw_frequencies.txt  - Raw frequencies for R analysis")
        print("  <prefix>.family_summary.txt   - Summary by degeneracy family")
        print()
        print("IMPORTANT: This analysis provides DESCRIPTIVE statistics only.")
        print("See output files for caveats and limitations.")
        sys.exit(1)
    
    start_time = time.time()
    
    input_file = sys.argv[1]
    output_prefix = sys.argv[2]
    
    # Validate input file exists
    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}", file=sys.stderr)
        sys.exit(1)
    
    # Check file is not empty
    if os.path.getsize(input_file) == 0:
        print(f"ERROR: Input file is empty: {input_file}", file=sys.stderr)
        sys.exit(1)
    
    print("=== Preferred Codon Frequency Analysis ===")
    print()
    print(f"Input file: {input_file}")
    print(f"Output prefix: {output_prefix}")
    print()
    
    preferred_freq_data, parse_errors = process_codon_frequencies(input_file)
    
    if not preferred_freq_data:
        print("ERROR: No data processed. Check input file format.", file=sys.stderr)
        sys.exit(1)
    
    write_output(preferred_freq_data, output_prefix)
    
    print_summary(preferred_freq_data)
    
    elapsed_time = time.time() - start_time
    print(f"\nRuntime: {elapsed_time:.1f} seconds")
    
    if parse_errors > 0:
        print(f"WARNING: {parse_errors} lines had parsing errors", file=sys.stderr)
    
    print("\nDone!")


if __name__ == "__main__":
    main()
