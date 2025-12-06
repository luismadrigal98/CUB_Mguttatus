#!/usr/bin/env python3
"""
Analyze allele frequency distributions at synonymous sites,
stratified by preferred/non-preferred codon status.

This script provides DESCRIPTIVE statistics on allele frequencies:
1. Frequency spectrum of changes creating preferred codons (nonpref→pref)
2. Frequency spectrum of changes creating non-preferred codons (pref→nonpref)

IMPORTANT LIMITATIONS:
- Without outgroup sequences, cannot distinguish selection from drift
- Frequency differences could reflect mutation bias, demography, or drift
- Cannot infer direction of selection without ancestral state polarization
- Results are descriptive and require additional tests for inference

Polarization:
  - Uses preferred/non-preferred status as a categorical axis (NOT ancestral/derived)
  - REF codon status defines the current reference state
  - ALT codon status defines the alternative allele state
  - Tracks frequency of the ALT allele at each polymorphic site

Output:
  - Allele frequency distributions binned by category
  - Separate distributions for each amino acid
  - Descriptive statistics (mean frequency, counts)

Input:
    - VCF file with variant sites
    - <chrom>.genic_bases.annotated.txt (with Strand column)
    - preferred_codons.txt

Author: Luis Javier Madrigal-Roca & GitHub Copilot
Date: 2025-12-05
"""

import sys
import math
from collections import defaultdict, Counter

# Genetic code
GENETIC_CODE = {
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

# Amino acid families by degeneracy
# M (Met) and W (Trp) are excluded - single codon, no synonymous variation possible
# Stop codons (*) are also excluded
AA_FAMILIES = {
    '2-fold': ['F', 'Y', 'H', 'Q', 'N', 'K', 'D', 'E', 'C'],
    '3-fold': ['I'],
    '4-fold': ['A', 'G', 'P', 'T', 'V'],
    '6-fold': ['L', 'S', 'R']
}

# Amino acids to skip (no synonymous variation)
SKIP_AMINO_ACIDS = {'*', 'M', 'W'}

# Reverse mapping: AA -> family
AA_TO_FAMILY = {}
for family, aas in AA_FAMILIES.items():
    for aa in aas:
        AA_TO_FAMILY[aa] = family


def complement(base):
    """Return complement of a DNA base."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return comp.get(base.upper(), 'N')


def load_preferred_codons(preferred_file):
    """Load list of preferred codons from file."""
    preferred = set()
    
    with open(preferred_file, 'r') as f:
        for line in f:
            codon = line.strip().upper()
            if len(codon) == 3 and codon in GENETIC_CODE:
                preferred.add(codon)
    
    return preferred


def load_annotated_sites(annotation_file, chrom):
    """
    Load position annotations from describe_gene_positions_by_degeneracy.py output.
    
    Returns:
        sites: dict of position -> {
            'gene': str,
            'base': str,
            'codon_pos': int,
            'degeneracy': str,
            'ref_codon': str,
            'amino_acid': str,
            'strand': str
        }
    """
    sites = {}
    
    with open(annotation_file, 'r') as f:
        header = f.readline().strip().split('\t')
        
        # Check if strand column exists
        has_strand = 'Strand' in header
        if not has_strand:
            print("  WARNING: Annotation file lacks Strand column.", file=sys.stderr)
            print("           Re-run describe_gene_positions_by_degeneracy.py", file=sys.stderr)
        
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 8:
                continue
            
            chr_name = cols[0]
            gene_id = cols[1]
            pos = cols[2]
            base = cols[3]
            codon_pos = cols[4]
            degeneracy = cols[5]
            ref_codon = cols[6]
            amino_acid = cols[7]
            strand = cols[8] if len(cols) > 8 else '+'
            
            if chr_name == chrom:
                sites[int(pos)] = {
                    'gene': gene_id,
                    'base': base,
                    'codon_pos': int(codon_pos),
                    'degeneracy': degeneracy,
                    'ref_codon': ref_codon,
                    'amino_acid': amino_acid,
                    'strand': strand
                }
    
    return sites


def change_base_in_codon(codon, position, new_base):
    """Change base at position (0, 1, or 2) in codon."""
    codon_list = list(codon)
    codon_list[position] = new_base
    return ''.join(codon_list)


def get_allele_frequency(genotypes, min_depth_ratio=5):
    """
    Calculate allele frequency from VCF genotypes.
    
    Uses strict depth-based calling (matching calculate_pi.py logic):
    - 0/0 with ref_depth > 5*alt_depth → count as REF homozygote
    - 1/1 with alt_depth > 5*ref_depth → count as ALT homozygote
    
    Returns:
        (alt_freq, n_chromosomes, ref_count, alt_count)
        or (None, 0, 0, 0) if not enough data
    """
    ref_hom = 0
    alt_hom = 0
    
    for gt, ref_depth, alt_depth in genotypes:
        if gt == "0/0" and ref_depth > min_depth_ratio * alt_depth:
            ref_hom += 1
        elif gt == "1/1" and alt_depth > min_depth_ratio * ref_depth:
            alt_hom += 1
    
    total = ref_hom + alt_hom
    
    if total < 2:
        return None, 0, 0, 0
    
    # For diploids, each individual = 2 chromosomes
    # But we're counting homozygotes, so:
    n_ref = ref_hom
    n_alt = alt_hom
    
    alt_freq = float(n_alt) / float(total)
    
    return alt_freq, total, n_ref, n_alt


def parse_vcf_line(line):
    """Parse a single VCF line."""
    cols = line.strip().split('\t')
    
    if len(cols) < 10:
        return None
    
    chrom = cols[0]
    pos = int(cols[1])
    ref = cols[3]
    alt = cols[4]
    
    genotypes = []
    for j in range(9, len(cols)):
        gt_field = cols[j]
        parts = gt_field.split(':')
        
        if len(parts) < 3:
            genotypes.append(("./.", 0, 0))
            continue
        
        gt = parts[0]
        
        try:
            ad_field = parts[2]
            ad_parts = ad_field.split(',')
            ref_count = int(ad_parts[0])
            alt_count = int(ad_parts[1]) if len(ad_parts) > 1 else 0
        except (ValueError, IndexError):
            ref_count, alt_count = 0, 0
        
        genotypes.append((gt, ref_count, alt_count))
    
    return chrom, pos, ref, alt, genotypes


def categorize_and_get_frequency(site_info, ref_base, alt_base, genotypes, preferred_codons):
    """
    Categorize a variant and get its allele frequency.
    
    Returns:
        (category, amino_acid, alt_freq, n_samples)
        
        category: 'to_preferred', 'to_nonpreferred', 'neutral_syn', 'non_synonymous', or None
    """
    ref_codon = site_info['ref_codon']
    codon_pos = site_info['codon_pos'] - 1  # Convert to 0-indexed
    amino_acid = site_info['amino_acid']
    strand = site_info.get('strand', '+')
    
    # Skip stop codons and single-codon amino acids
    if amino_acid in ['*', 'M', 'W']:
        return None, amino_acid, None, 0
    
    # Convert alt_base from reference strand to gene orientation
    if strand == '-':
        alt_base_gene = complement(alt_base)
    else:
        alt_base_gene = alt_base
    
    # Create alternate codon
    alt_codon = change_base_in_codon(ref_codon, codon_pos, alt_base_gene)
    
    # Check validity
    if alt_codon not in GENETIC_CODE:
        return None, amino_acid, None, 0
    
    alt_aa = GENETIC_CODE[alt_codon]
    
    # Skip non-synonymous
    if alt_aa != amino_acid:
        return 'non_synonymous', amino_acid, None, 0
    
    # Get allele frequency
    alt_freq, n_samples, n_ref, n_alt = get_allele_frequency(genotypes)
    
    if alt_freq is None:
        return None, amino_acid, None, 0
    
    # Categorize by preferred status
    ref_is_pref = ref_codon in preferred_codons
    alt_is_pref = alt_codon in preferred_codons
    
    if ref_is_pref and not alt_is_pref:
        # REF is preferred, ALT is non-preferred
        # ALT frequency = frequency of non-preferred allele
        # This is "away from preferred" or "to_nonpreferred"
        category = 'to_nonpreferred'
    elif not ref_is_pref and alt_is_pref:
        # REF is non-preferred, ALT is preferred
        # ALT frequency = frequency of preferred allele
        # This is "toward preferred" or "to_preferred"
        category = 'to_preferred'
    else:
        # Both same status (both preferred or both non-preferred)
        category = 'neutral_syn'
    
    return category, amino_acid, alt_freq, n_samples


def process_vcf(vcf_file, sites, preferred_codons, chrom, n_bins=20):
    """
    Process VCF file and build SFS per amino acid.
    
    Returns:
        sfs_data: dict of amino_acid -> {
            'to_preferred': [freq_list],
            'to_nonpreferred': [freq_list],
            'neutral_syn': [freq_list],
            'invariant_pref': int,
            'invariant_nonpref': int
        }
    """
    # Initialize SFS data structure
    sfs_data = defaultdict(lambda: {
        'to_preferred': [],
        'to_nonpreferred': [],
        'neutral_syn': [],
        'invariant_pref': 0,
        'invariant_nonpref': 0
    })
    
    # Count invariant sites per amino acid
    print("Counting invariant sites by amino acid...")
    for pos, site_info in sites.items():
        aa = site_info['amino_acid']
        if aa in ['*', 'M', 'W']:
            continue
        
        ref_codon = site_info['ref_codon']
        if ref_codon in preferred_codons:
            sfs_data[aa]['invariant_pref'] += 1
        else:
            sfs_data[aa]['invariant_nonpref'] += 1
    
    # Process polymorphic sites
    print(f"Processing VCF for {chrom}...")
    line_count = 0
    poly_count = 0
    
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            line_count += 1
            if line_count % 500000 == 0:
                print(f"  Processed {line_count:,} VCF lines, {poly_count:,} synonymous polymorphisms...", 
                      file=sys.stderr)
            
            parsed = parse_vcf_line(line)
            if not parsed:
                continue
            
            vcf_chrom, pos, ref, alt, genotypes = parsed
            
            if vcf_chrom != chrom:
                continue
            
            if alt == '.':
                continue  # Invariant
            
            if pos not in sites:
                continue  # Not in CDS
            
            site_info = sites[pos]
            
            category, aa, alt_freq, n_samples = categorize_and_get_frequency(
                site_info, ref, alt, genotypes, preferred_codons
            )
            
            if category is None or alt_freq is None:
                continue
            
            if category == 'non_synonymous':
                continue
            
            poly_count += 1
            
            # Subtract from invariant count (this site is polymorphic)
            ref_codon = site_info['ref_codon']
            if ref_codon in preferred_codons:
                sfs_data[aa]['invariant_pref'] -= 1
            else:
                sfs_data[aa]['invariant_nonpref'] -= 1
            
            # Add to appropriate frequency list
            sfs_data[aa][category].append(alt_freq)
    
    print(f"  Total VCF lines: {line_count:,}")
    print(f"  Synonymous polymorphisms: {poly_count:,}")
    
    return dict(sfs_data)


def calculate_sfs_bins(freq_list, n_bins=10):
    """
    Bin frequencies into SFS categories.
    
    For unfolded SFS with preferred/non-preferred polarization:
    - Bins from 0 to 1 in n_bins equal intervals
    
    Returns:
        counts: list of counts per bin
        bin_edges: list of bin boundaries
    """
    bin_edges = [i / n_bins for i in range(n_bins + 1)]
    counts = [0] * n_bins
    
    for freq in freq_list:
        # Find bin (handle edge case of freq=1.0)
        bin_idx = min(int(freq * n_bins), n_bins - 1)
        counts[bin_idx] += 1
    
    return counts, bin_edges


def calculate_mean_frequency(freq_list):
    """Calculate mean allele frequency."""
    if not freq_list:
        return 0.0
    return sum(freq_list) / len(freq_list)


def write_output(sfs_data, chrom, output_prefix):
    """Write allele frequency distribution results to output files."""
    
    # 1. Per-amino-acid summary
    summary_file = f"{output_prefix}.aa_summary.txt"
    print(f"Writing amino acid summary to {summary_file}...")
    
    with open(summary_file, 'w') as out:
        # Add caveats header
        out.write("# Allele Frequency Distribution by Amino Acid\n")
        out.write("# \n")
        out.write("# IMPORTANT CAVEATS:\n")
        out.write("#   - These are DESCRIPTIVE statistics only\n")
        out.write("#   - Without outgroup: cannot distinguish selection from drift\n")
        out.write("#   - Frequency differences could reflect mutation bias, demography, or drift\n")
        out.write("# \n")
        out.write("Chr\tAmino_Acid\tFamily\t")
        out.write("N_to_preferred\tN_to_nonpreferred\tN_neutral_syn\t")
        out.write("Mean_freq_to_pref\tMean_freq_to_nonpref\tMean_freq_neutral\t")
        out.write("Invariant_pref\tInvariant_nonpref\n")
        
        for aa in sorted(sfs_data.keys()):
            data = sfs_data[aa]
            family = AA_TO_FAMILY.get(aa, 'unknown')
            
            n_to_pref = len(data['to_preferred'])
            n_to_nonpref = len(data['to_nonpreferred'])
            n_neutral = len(data['neutral_syn'])
            
            mean_to_pref = calculate_mean_frequency(data['to_preferred'])
            mean_to_nonpref = calculate_mean_frequency(data['to_nonpreferred'])
            mean_neutral = calculate_mean_frequency(data['neutral_syn'])
            
            out.write(f"{chrom}\t{aa}\t{family}\t")
            out.write(f"{n_to_pref}\t{n_to_nonpref}\t{n_neutral}\t")
            out.write(f"{mean_to_pref:.4f}\t{mean_to_nonpref:.4f}\t{mean_neutral:.4f}\t")
            out.write(f"{data['invariant_pref']}\t{data['invariant_nonpref']}\n")
    
    # 2. Detailed SFS per amino acid (binned)
    sfs_file = f"{output_prefix}.sfs_binned.txt"
    print(f"Writing binned SFS to {sfs_file}...")
    
    n_bins = 10
    
    with open(sfs_file, 'w') as out:
        # Header
        out.write("Chr\tAmino_Acid\tFamily\tCategory\t")
        bin_headers = [f"Bin_{i+1}" for i in range(n_bins)]
        out.write("\t".join(bin_headers) + "\n")
        
        for aa in sorted(sfs_data.keys()):
            data = sfs_data[aa]
            family = AA_TO_FAMILY.get(aa, 'unknown')
            
            for category in ['to_preferred', 'to_nonpreferred', 'neutral_syn']:
                counts, _ = calculate_sfs_bins(data[category], n_bins)
                
                out.write(f"{chrom}\t{aa}\t{family}\t{category}\t")
                out.write("\t".join(map(str, counts)) + "\n")
    
    # 3. Raw frequencies (for custom analysis)
    raw_file = f"{output_prefix}.sfs_raw.txt"
    print(f"Writing raw frequencies to {raw_file}...")
    
    with open(raw_file, 'w') as out:
        out.write("Chr\tAmino_Acid\tFamily\tCategory\tAlt_Frequency\n")
        
        for aa in sorted(sfs_data.keys()):
            data = sfs_data[aa]
            family = AA_TO_FAMILY.get(aa, 'unknown')
            
            for category in ['to_preferred', 'to_nonpreferred', 'neutral_syn']:
                for freq in data[category]:
                    out.write(f"{chrom}\t{aa}\t{family}\t{category}\t{freq:.4f}\n")
    
    # 4. Summary by family
    family_file = f"{output_prefix}.family_summary.txt"
    print(f"Writing family summary to {family_file}...")
    
    # Aggregate by family
    family_data = defaultdict(lambda: {
        'to_preferred': [],
        'to_nonpreferred': [],
        'neutral_syn': [],
        'invariant_pref': 0,
        'invariant_nonpref': 0
    })
    
    for aa, data in sfs_data.items():
        family = AA_TO_FAMILY.get(aa, 'unknown')
        family_data[family]['to_preferred'].extend(data['to_preferred'])
        family_data[family]['to_nonpreferred'].extend(data['to_nonpreferred'])
        family_data[family]['neutral_syn'].extend(data['neutral_syn'])
        family_data[family]['invariant_pref'] += data['invariant_pref']
        family_data[family]['invariant_nonpref'] += data['invariant_nonpref']
    
    with open(family_file, 'w') as out:
        out.write("Chr\tFamily\t")
        out.write("N_to_preferred\tN_to_nonpreferred\tN_neutral_syn\t")
        out.write("Mean_freq_to_pref\tMean_freq_to_nonpref\tMean_freq_neutral\t")
        out.write("Invariant_pref\tInvariant_nonpref\n")
        
        for family in ['2-fold', '3-fold', '4-fold', '6-fold']:
            data = family_data[family]
            
            n_to_pref = len(data['to_preferred'])
            n_to_nonpref = len(data['to_nonpreferred'])
            n_neutral = len(data['neutral_syn'])
            
            mean_to_pref = calculate_mean_frequency(data['to_preferred'])
            mean_to_nonpref = calculate_mean_frequency(data['to_nonpreferred'])
            mean_neutral = calculate_mean_frequency(data['neutral_syn'])
            
            out.write(f"{chrom}\t{family}\t")
            out.write(f"{n_to_pref}\t{n_to_nonpref}\t{n_neutral}\t")
            out.write(f"{mean_to_pref:.4f}\t{mean_to_nonpref:.4f}\t{mean_neutral:.4f}\t")
            out.write(f"{data['invariant_pref']}\t{data['invariant_nonpref']}\n")


def print_summary(sfs_data):
    """Print descriptive summary statistics with appropriate caveats."""
    
    print("\n" + "=" * 70)
    print("ALLELE FREQUENCY DISTRIBUTION - DESCRIPTIVE SUMMARY")
    print("=" * 70)
    
    # Aggregate totals
    total_to_pref = 0
    total_to_nonpref = 0
    total_neutral = 0
    all_to_pref_freqs = []
    all_to_nonpref_freqs = []
    
    for aa, data in sfs_data.items():
        total_to_pref += len(data['to_preferred'])
        total_to_nonpref += len(data['to_nonpreferred'])
        total_neutral += len(data['neutral_syn'])
        all_to_pref_freqs.extend(data['to_preferred'])
        all_to_nonpref_freqs.extend(data['to_nonpreferred'])
    
    print(f"\nTotal synonymous polymorphisms:")
    print(f"  Creating preferred codon (nonpref→pref):     {total_to_pref:,}")
    print(f"  Creating non-preferred codon (pref→nonpref): {total_to_nonpref:,}")
    print(f"  Same preference status:                      {total_neutral:,}")
    
    if all_to_pref_freqs and all_to_nonpref_freqs:
        mean_to_pref = sum(all_to_pref_freqs) / len(all_to_pref_freqs)
        mean_to_nonpref = sum(all_to_nonpref_freqs) / len(all_to_nonpref_freqs)
        
        print(f"\nMean ALT allele frequency (DESCRIPTIVE):")
        print(f"  Creating preferred:     {mean_to_pref:.4f}")
        print(f"  Creating non-preferred: {mean_to_nonpref:.4f}")
        
        diff = abs(mean_to_pref - mean_to_nonpref)
        print(f"\n--- What This Means ---")
        if mean_to_pref > mean_to_nonpref:
            print(f"  ALT alleles creating preferred codons are at higher frequency (+{diff:.4f})")
        elif mean_to_nonpref > mean_to_pref:
            print(f"  ALT alleles creating non-preferred codons are at higher frequency (+{diff:.4f})")
        else:
            print(f"  Frequencies are approximately equal")
        
        print(f"\n--- IMPORTANT CAVEATS ---")
        print(f"  ✗ This is DESCRIPTIVE data only")
        print(f"  ✗ Cannot infer selection without:")
        print(f"      • Outgroup sequences (to polarize ancestral/derived)")
        print(f"      • Neutral baseline comparison")
        print(f"      • Demographic model")
        print(f"  ✗ Differences could reflect:")
        print(f"      • Selection (for or against preferred codons)")
        print(f"      • Mutation bias (toward/away from preferred)")
        print(f"      • Demographic history")
        print(f"      • Random drift")
    
    print("\n" + "=" * 70)


def main():
    if len(sys.argv) != 5:
        print("Usage: python calculate_sfs_by_amino_acid.py <chromosome> <vcf_file> <annotation_file> <preferred_codons>")
        print()
        print("Example:")
        print("  python calculate_sfs_by_amino_acid.py Chr_01 variants.vcf Chr_01.genic_bases.annotated.txt preferred_codons.txt")
        print()
        print("Output files:")
        print("  <chrom>.sfs.aa_summary.txt     - Per-amino-acid summary statistics")
        print("  <chrom>.sfs.sfs_binned.txt     - Binned SFS per amino acid per category")
        print("  <chrom>.sfs.sfs_raw.txt        - Raw frequencies for custom analysis")
        print("  <chrom>.sfs.family_summary.txt - Summary by degeneracy family")
        sys.exit(1)
    
    chrom = sys.argv[1]
    vcf_file = sys.argv[2]
    annotation_file = sys.argv[3]
    preferred_file = sys.argv[4]
    
    print(f"=== SFS Analysis for {chrom} ===")
    print()
    
    print("Loading preferred codons...")
    preferred_codons = load_preferred_codons(preferred_file)
    print(f"  Loaded {len(preferred_codons)} preferred codons")
    
    print("Loading annotated sites...")
    sites = load_annotated_sites(annotation_file, chrom)
    print(f"  Loaded {len(sites):,} annotated positions")
    
    print("Processing VCF and building SFS...")
    sfs_data = process_vcf(vcf_file, sites, preferred_codons, chrom)
    
    output_prefix = f"{chrom}.sfs"
    write_output(sfs_data, chrom, output_prefix)
    
    print_summary(sfs_data)
    
    print("\nDone!")


if __name__ == "__main__":
    main()
