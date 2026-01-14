#!/usr/bin/env python3
"""
Calculate nucleotide diversity (π) across different genomic compartments.

Compartments:
    1. Intergenic windows (50kb) - neutral baseline
    2. Introns (50bp trimmed from splice sites) - neutral within genes
    3. First exons - 4-fold degenerate sites only
    4. Non-first exons - 4-fold degenerate sites only

Input:
    - VCF file with variant AND invariant sites
    - GFF3 annotation file
    - Degeneracy annotation file (from describe_gene_positions_by_degeneracy.py)

Output:
    - pi_by_compartment.txt: Summary statistics per compartment
    - pi_by_window.txt: Per-window/region statistics

Author: Luis Javier Madrigal-Roca & GitHub Copilot
Date: 2026-01-14
"""

import sys
import os
import re
import argparse
import gzip
from collections import defaultdict
import math

# ====================== CONSTANTS ======================

INTRON_TRIM_BP = 50       # bp to trim from intron boundaries (splice site removal)
INTERGENIC_WINDOW = 50000 # 50kb windows for intergenic regions
MIN_SAMPLES = 10          # Minimum sample size for π calculation
MIN_DEPTH_RATIO = 5       # Depth ratio filter for homozygous calls

# Nucleotide categories for polarized analysis
NUC_CATEGORIES = ['all', 'C', 'G', 'AT']  # AT = A or T (weak nucleotides)


# ====================== NUCLEOTIDE POLARIZATION ======================

def get_complement(nuc):
    """Return complement of a nucleotide."""
    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', '.': '.'}
    return mapping.get(nuc.upper(), 'N')


def classify_nucleotide(nuc):
    """
    Classify nucleotide for polarized π analysis.
    Returns: 'C', 'G', 'AT', or None for invalid
    """
    nuc = nuc.upper()
    if nuc == 'C':
        return 'C'
    elif nuc == 'G':
        return 'G'
    elif nuc in ['A', 'T']:
        return 'AT'
    return None


def get_site_nucleotide_category(ref, alt, strand, is_invariant):
    """
    Determine which nucleotide category a site belongs to.
    For polarized π, we want to know if the site has C, G, or AT alleles.
    
    For invariant sites: category based on ref allele
    For polymorphic sites: we track BOTH alleles (site contributes to both categories)
    
    Returns: list of categories this site belongs to ['C'], ['G'], ['AT'], ['C', 'AT'], etc.
    """
    # Apply strand correction
    if strand == '-':
        ref_corrected = get_complement(ref)
        alt_corrected = get_complement(alt) if alt and alt not in ['.', '*', '<NON_REF>'] else None
    else:
        ref_corrected = ref.upper()
        alt_corrected = alt.upper() if alt and alt not in ['.', '*', '<NON_REF>'] else None
    
    categories = set()
    
    # Add ref category
    ref_cat = classify_nucleotide(ref_corrected)
    if ref_cat:
        categories.add(ref_cat)
    
    # Add alt category (if polymorphic)
    if not is_invariant and alt_corrected:
        alt_cat = classify_nucleotide(alt_corrected)
        if alt_cat:
            categories.add(alt_cat)
    
    return list(categories)


# ====================== GFF PARSING ======================

def parse_gff3(gff_file, target_chrom=None):
    """
    Parse GFF3 to extract:
    - Gene boundaries (for defining intergenic)
    - Exon coordinates with exon number (for first vs non-first)
    - Calculate introns from exon gaps
    
    Args:
        gff_file: Path to GFF3 file
        target_chrom: If set, only parse this chromosome (for parallel processing)
    
    Returns:
        genes: {gene_id: {'chrom': str, 'start': int, 'end': int, 'strand': str}}
        exons: {chrom: [(start, end, gene_id, exon_number, strand), ...]}
        introns: {chrom: [(start, end, gene_id, strand), ...]}  # After trimming
        gene_intervals: {chrom: [(start, end), ...]}  # For intergenic calculation
    """
    chrom_filter = f" (filtering for {target_chrom})" if target_chrom else ""
    print(f"Parsing GFF3: {gff_file}{chrom_filter}")
    
    genes = {}
    mrna_exons = defaultdict(list)  # mrna_id -> [(start, end)]
    mrna_to_gene = {}  # mrna_id -> gene_id
    gene_info = {}  # gene_id -> {chrom, strand}
    
    # First pass: collect all features
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            chrom, source, feature, start, end, score, strand, phase, attributes = parts
            
            # Skip scaffolds
            if "Chr" not in chrom:
                continue
            
            # Filter by target chromosome if specified
            if target_chrom and chrom != target_chrom:
                continue
            
            start, end = int(start) - 1, int(end)  # Convert to 0-based
            
            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, val = attr.split('=', 1)
                    attr_dict[key] = val
            
            if feature == 'gene':
                gene_id = attr_dict.get('ID', '')
                if gene_id:
                    genes[gene_id] = {'chrom': chrom, 'start': start, 'end': end, 'strand': strand}
                    gene_info[gene_id] = {'chrom': chrom, 'strand': strand}
            
            elif feature == 'mRNA':
                mrna_id = attr_dict.get('ID', '')
                parent = attr_dict.get('Parent', '')
                if mrna_id and parent:
                    mrna_to_gene[mrna_id] = parent
            
            elif feature in ['CDS', 'exon']:
                parent = attr_dict.get('Parent', '')
                if parent:
                    mrna_exons[parent].append((start, end))
    
    print(f"  Found {len(genes)} genes, {len(mrna_exons)} mRNAs with CDS/exon features")
    
    # Second pass: calculate exons with numbering and introns
    exons = defaultdict(list)
    introns = defaultdict(list)
    gene_intervals = defaultdict(list)
    
    # Build gene intervals for intergenic calculation
    for gene_id, info in genes.items():
        gene_intervals[info['chrom']].append((info['start'], info['end']))
    
    # Sort gene intervals per chromosome
    for chrom in gene_intervals:
        gene_intervals[chrom].sort()
    
    # Process each mRNA
    for mrna_id, exon_list in mrna_exons.items():
        gene_id = mrna_to_gene.get(mrna_id)
        if not gene_id or gene_id not in gene_info:
            continue
        
        chrom = gene_info[gene_id]['chrom']
        strand = gene_info[gene_id]['strand']
        
        # Sort exons by position
        sorted_exons = sorted(exon_list, key=lambda x: x[0])
        
        # Assign exon numbers (respecting strand)
        if strand == '-':
            # For minus strand, first exon is the last in genomic coordinates
            for i, (ex_start, ex_end) in enumerate(reversed(sorted_exons)):
                exon_num = i + 1
                exons[chrom].append((ex_start, ex_end, gene_id, exon_num, strand))
        else:
            # For plus strand, first exon is the first in genomic coordinates
            for i, (ex_start, ex_end) in enumerate(sorted_exons):
                exon_num = i + 1
                exons[chrom].append((ex_start, ex_end, gene_id, exon_num, strand))
        
        # Calculate introns (gaps between exons, with trimming)
        if len(sorted_exons) >= 2:
            for i in range(len(sorted_exons) - 1):
                intron_start = sorted_exons[i][1] + INTRON_TRIM_BP
                intron_end = sorted_exons[i + 1][0] - INTRON_TRIM_BP
                
                if intron_end > intron_start:
                    introns[chrom].append((intron_start, intron_end, gene_id, strand))
    
    # Sort for binary search
    for chrom in exons:
        exons[chrom].sort()
    for chrom in introns:
        introns[chrom].sort()
    
    # Count statistics
    total_exons = sum(len(v) for v in exons.values())
    total_introns = sum(len(v) for v in introns.values())
    first_exons = sum(1 for chrom in exons for ex in exons[chrom] if ex[3] == 1)
    
    print(f"  Exons: {total_exons} total, {first_exons} first exons")
    print(f"  Introns: {total_introns} (after {INTRON_TRIM_BP}bp trimming)")
    
    return genes, exons, introns, gene_intervals


def calculate_intergenic_windows(gene_intervals, chrom_sizes=None):
    """
    Calculate intergenic regions and divide into windows.
    
    Returns:
        {chrom: [(window_start, window_end), ...]}
    """
    intergenic_windows = defaultdict(list)
    
    for chrom, intervals in gene_intervals.items():
        if not intervals:
            continue
        
        # Sort intervals
        sorted_intervals = sorted(intervals)
        
        # Calculate intergenic gaps
        intergenic_regions = []
        
        # Gap before first gene (if we know chrom size)
        # For now, we start from position 0 or after the first gene
        prev_end = 0
        
        for gene_start, gene_end in sorted_intervals:
            if gene_start > prev_end:
                intergenic_regions.append((prev_end, gene_start))
            prev_end = max(prev_end, gene_end)
        
        # Split large intergenic regions into windows
        for ig_start, ig_end in intergenic_regions:
            ig_length = ig_end - ig_start
            
            if ig_length < INTERGENIC_WINDOW:
                # Small gap - keep as single window
                intergenic_windows[chrom].append((ig_start, ig_end))
            else:
                # Divide into windows
                for win_start in range(ig_start, ig_end, INTERGENIC_WINDOW):
                    win_end = min(win_start + INTERGENIC_WINDOW, ig_end)
                    intergenic_windows[chrom].append((win_start, win_end))
    
    total_windows = sum(len(v) for v in intergenic_windows.values())
    print(f"  Intergenic: {total_windows} windows of up to {INTERGENIC_WINDOW // 1000}kb")
    
    return intergenic_windows


# ====================== DEGENERACY ANNOTATION ======================

def load_degeneracy_annotations(annotation_files):
    """
    Load position-level degeneracy annotations.
    
    Input format (from describe_gene_positions_by_degeneracy.py):
        Chr  Gene  Position  Base  Codon_Position  Degeneracy  Ref_Codon  Amino_Acid  [Strand]
    
    Returns:
        {chrom: {pos: {'gene': str, 'degeneracy': str, 'base': str}}}
    """
    degeneracy = defaultdict(dict)
    
    for annot_file in annotation_files:
        print(f"Loading degeneracy: {annot_file}")
        with open(annot_file, 'r') as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                
                chrom = parts[0]
                gene_id = parts[1]
                pos = int(parts[2])
                base = parts[3]
                deg_class = parts[5]  # "0-fold", "2-fold", "3-fold", "4-fold"
                
                degeneracy[chrom][pos] = {
                    'gene': gene_id,
                    'degeneracy': deg_class,
                    'base': base
                }
    
    total = sum(len(v) for v in degeneracy.values())
    print(f"  Loaded {total:,} annotated positions")
    
    return degeneracy


# ====================== POSITION LOOKUP ======================

def binary_search_region(chrom, pos, regions):
    """
    Check if position falls within any region.
    regions: list of (start, end, ...) tuples, sorted by start
    
    Returns: matching region tuple or None
    """
    if chrom not in regions:
        return None
    
    region_list = regions[chrom]
    low, high = 0, len(region_list) - 1
    
    while low <= high:
        mid = (low + high) // 2
        start = region_list[mid][0]
        end = region_list[mid][1]
        
        if start <= pos < end:
            return region_list[mid]
        elif pos < start:
            high = mid - 1
        else:
            low = mid + 1
    
    return None


def classify_position(chrom, pos, exons, introns, intergenic, degeneracy):
    """
    Classify a position into one of the compartments.
    
    Returns:
        (compartment, details)
        compartment: 'intergenic' | 'intron' | 'first_exon_4fold' | 'nonfirst_exon_4fold' | None
        details: dict with additional info
    """
    # Check exons first (need 4-fold)
    exon_match = binary_search_region(chrom, pos, exons)
    if exon_match:
        ex_start, ex_end, gene_id, exon_num, strand = exon_match
        
        # Check if 4-fold degenerate
        if chrom in degeneracy and pos in degeneracy[chrom]:
            deg_info = degeneracy[chrom][pos]
            if deg_info['degeneracy'] == '4-fold':
                if exon_num == 1:
                    return 'first_exon_4fold', {'gene': gene_id, 'exon': exon_num}
                else:
                    return 'nonfirst_exon_4fold', {'gene': gene_id, 'exon': exon_num}
        return None, {}  # Exon but not 4-fold
    
    # Check introns
    intron_match = binary_search_region(chrom, pos, introns)
    if intron_match:
        in_start, in_end, gene_id, strand = intron_match
        return 'intron', {'gene': gene_id}
    
    # Check intergenic
    ig_match = binary_search_region(chrom, pos, intergenic)
    if ig_match:
        ig_start, ig_end = ig_match[:2]
        return 'intergenic', {'window': f"{chrom}:{ig_start}-{ig_end}"}
    
    return None, {}


# ====================== VCF PROCESSING ======================

def parse_vcf_line(line):
    """
    Parse a VCF line and extract genotype information.
    Returns: (chrom, pos, ref, alt, genotypes) or None
    """
    parts = line.strip().split('\t')
    if len(parts) < 10:
        return None
    
    chrom = parts[0]
    pos = int(parts[1]) - 1  # Convert to 0-based
    ref = parts[3]
    alt = parts[4]
    
    # Parse FORMAT field
    fmt = parts[8].split(':')
    try:
        gt_idx = fmt.index('GT')
        ad_idx = fmt.index('AD')
    except ValueError:
        return None
    
    genotypes = []
    for sample_field in parts[9:]:
        sample_parts = sample_field.split(':')
        if len(sample_parts) <= max(gt_idx, ad_idx):
            genotypes.append(('./.', 0, 0))
            continue
        
        gt = sample_parts[gt_idx]
        ad_str = sample_parts[ad_idx]
        
        try:
            if ',' in ad_str:
                ad_parts = ad_str.split(',')
                ref_count = int(ad_parts[0])
                alt_count = int(ad_parts[1]) if len(ad_parts) > 1 else 0
            else:
                ref_count = int(ad_str) if ad_str and ad_str != '.' else 0
                alt_count = 0
        except ValueError:
            ref_count, alt_count = 0, 0
        
        genotypes.append((gt, ref_count, alt_count))
    
    return chrom, pos, ref, alt, genotypes


def calculate_pi_site(genotypes):
    """
    Calculate π for a single site using homozygous calls only.
    Matches the logic in calculate_pi.py
    """
    ref_hom = 0
    alt_hom = 0
    
    for gt, ref_count, alt_count in genotypes:
        if gt == '0/0' and ref_count > MIN_DEPTH_RATIO * alt_count:
            ref_hom += 1
        elif gt == '1/1' and alt_count > MIN_DEPTH_RATIO * ref_count:
            alt_hom += 1
    
    if min(ref_hom, alt_hom) > 0:
        nx = float(ref_hom + alt_hom)
        px = float(ref_hom) / nx
        pi = 2.0 * nx * px * (1.0 - px) / (nx - 1.0)
        return True, pi, int(nx)
    else:
        return False, 0.0, 0


def init_compartment_stats():
    """
    Initialize stats dictionary with all compartments and nucleotide categories.
    Structure: {compartment: {nuc_category: {'sites': int, 'poly': int, 'pi_sum': float}}}
    """
    compartments = ['intergenic', 'intron', 'first_exon_4fold', 'nonfirst_exon_4fold']
    stats = {}
    for comp in compartments:
        stats[comp] = {}
        for nuc in NUC_CATEGORIES:
            stats[comp][nuc] = {'sites': 0, 'poly': 0, 'pi_sum': 0.0}
    return stats


def process_vcf(vcf_path, exons, introns, intergenic, degeneracy, stream=False, target_chrom=None):
    """
    Process VCF and calculate π for each compartment, with C/G polarization.
    
    Args:
        vcf_path: Path to VCF file
        exons, introns, intergenic, degeneracy: Annotation data structures
        stream: If True, read from stdin
        target_chrom: If set, only process this chromosome (for parallel processing)
    
    Returns:
        compartment_stats: {compartment: {nuc_category: {'sites': int, 'poly': int, 'pi_sum': float}}}
        
    Where nuc_category is 'all', 'C', 'G', or 'AT'
    """
    stats = init_compartment_stats()
    
    # Open input
    if stream:
        input_handle = sys.stdin
        chrom_filter = f" (filtering for {target_chrom})" if target_chrom else ""
        print(f"Reading VCF from stdin...{chrom_filter}")
    elif vcf_path.endswith('.gz'):
        input_handle = gzip.open(vcf_path, 'rt')
        print(f"Reading VCF: {vcf_path}")
    else:
        input_handle = open(vcf_path, 'r')
        print(f"Reading VCF: {vcf_path}")
    
    line_count = 0
    classified_count = 0
    
    for line in input_handle:
        if line.startswith('#'):
            continue
        
        line_count += 1
        if line_count % 1000000 == 0:
            print(f"  Processed {line_count:,} sites, classified {classified_count:,}...", 
                  file=sys.stderr)
        
        parsed = parse_vcf_line(line)
        if not parsed:
            continue
        
        chrom, pos, ref, alt, genotypes = parsed
        
        # Filter by target chromosome if specified
        if target_chrom and chrom != target_chrom:
            continue
        
        # Skip multi-allelic or indels
        is_invariant = (alt == '.' or alt == '<NON_REF>' or alt == '*')
        if not is_invariant and (len(ref) > 1 or len(alt) > 1 or ',' in alt):
            continue
        
        # Classify position
        compartment, details = classify_position(
            chrom, pos, exons, introns, intergenic, degeneracy
        )
        
        if compartment is None:
            continue
        
        classified_count += 1
        
        # Determine strand for polarization
        # For exons/introns, get strand from details or region match
        strand = '+'
        if compartment in ['first_exon_4fold', 'nonfirst_exon_4fold']:
            # Re-check exon to get strand
            exon_match = binary_search_region(chrom, pos, exons)
            if exon_match:
                strand = exon_match[4]  # (start, end, gene_id, exon_num, strand)
        elif compartment == 'intron':
            intron_match = binary_search_region(chrom, pos, introns)
            if intron_match:
                strand = intron_match[3]  # (start, end, gene_id, strand)
        # Intergenic: use '+' (no strand correction needed)
        
        # Get nucleotide categories for this site
        nuc_categories = get_site_nucleotide_category(ref, alt, strand, is_invariant)
        
        # Update stats
        if is_invariant:
            # Update 'all' category
            stats[compartment]['all']['sites'] += 1
            # Update specific nucleotide categories
            for nuc_cat in nuc_categories:
                stats[compartment][nuc_cat]['sites'] += 1
        else:
            is_poly, pi_val, n_samples = calculate_pi_site(genotypes)
            if is_poly:
                # Update 'all' category
                stats[compartment]['all']['sites'] += 1
                stats[compartment]['all']['poly'] += 1
                stats[compartment]['all']['pi_sum'] += pi_val
                # Update specific nucleotide categories
                for nuc_cat in nuc_categories:
                    stats[compartment][nuc_cat]['sites'] += 1
                    stats[compartment][nuc_cat]['poly'] += 1
                    stats[compartment][nuc_cat]['pi_sum'] += pi_val
            else:
                # Not polymorphic by our criteria
                stats[compartment]['all']['sites'] += 1
                for nuc_cat in nuc_categories:
                    stats[compartment][nuc_cat]['sites'] += 1
    
    if not stream and input_handle != sys.stdin:
        input_handle.close()
    
    print(f"  Total sites processed: {line_count:,}")
    print(f"  Sites classified: {classified_count:,}")
    
    return stats


# ====================== OUTPUT ======================

def write_summary(stats, output_file, chromosome=None):
    """Write summary statistics per compartment with C/G polarization."""
    
    chrom_str = chromosome if chromosome else "all"
    print(f"\nWriting summary to: {output_file}")
    
    compartments = ['intergenic', 'intron', 'first_exon_4fold', 'nonfirst_exon_4fold']
    
    with open(output_file, 'w') as out:
        # Header with chromosome column for merging across parallel jobs
        out.write("Chromosome\tCompartment\tNuc_Category\tSites\tPolymorphic\tPi_sum\tPi_mean\tPoly_fraction\n")
        
        for compartment in compartments:
            for nuc_cat in NUC_CATEGORIES:
                s = stats[compartment][nuc_cat]
                n_sites = s['sites']
                n_poly = s['poly']
                pi_sum = s['pi_sum']
                pi_mean = pi_sum / n_sites if n_sites > 0 else 0.0
                poly_frac = n_poly / n_sites if n_sites > 0 else 0.0
                
                out.write(f"{chrom_str}\t{compartment}\t{nuc_cat}\t{n_sites}\t{n_poly}\t{pi_sum:.6f}\t{pi_mean:.8f}\t{poly_frac:.6f}\n")
    
    # Also print to console - formatted nicely
    print("\n" + "=" * 90)
    print(f"SUMMARY: Nucleotide Diversity (π) by Genomic Compartment - {chrom_str}")
    print("With C/G/AT Polarization (strand-corrected)")
    print("=" * 90)
    
    # Print overall stats first
    print("\n--- OVERALL (all nucleotides) ---")
    print(f"{'Compartment':<25} {'Sites':>12} {'Poly':>10} {'π mean':>12}")
    print("-" * 60)
    for compartment in compartments:
        s = stats[compartment]['all']
        n_sites = s['sites']
        n_poly = s['poly']
        pi_mean = s['pi_sum'] / n_sites if n_sites > 0 else 0.0
        print(f"{compartment:<25} {n_sites:>12,} {n_poly:>10,} {pi_mean:>12.8f}")
    
    # Print C vs G vs AT comparison
    print("\n--- BY NUCLEOTIDE CATEGORY ---")
    print(f"{'Compartment':<25} {'π_C':>12} {'π_G':>12} {'π_AT':>12} {'C/AT ratio':>12}")
    print("-" * 75)
    for compartment in compartments:
        pi_C = stats[compartment]['C']['pi_sum'] / stats[compartment]['C']['sites'] if stats[compartment]['C']['sites'] > 0 else 0.0
        pi_G = stats[compartment]['G']['pi_sum'] / stats[compartment]['G']['sites'] if stats[compartment]['G']['sites'] > 0 else 0.0
        pi_AT = stats[compartment]['AT']['pi_sum'] / stats[compartment]['AT']['sites'] if stats[compartment]['AT']['sites'] > 0 else 0.0
        ratio_C_AT = pi_C / pi_AT if pi_AT > 0 else float('nan')
        print(f"{compartment:<25} {pi_C:>12.8f} {pi_G:>12.8f} {pi_AT:>12.8f} {ratio_C_AT:>12.4f}")
    
    print("=" * 90)
    print("\nNote: C and G categories are strand-corrected for coding regions.")
    print("      AT = sites with A or T alleles (weak nucleotides).")
    print("      Ratio > 1 indicates elevated π at C sites relative to AT.")


# ====================== MAIN ======================

def main():
    parser = argparse.ArgumentParser(
        description="Calculate π across genomic compartments (intergenic, intron, first/non-first exon 4-fold sites)"
    )
    parser.add_argument('--vcf', type=str, help='Path to VCF file (gzip supported)')
    parser.add_argument('--stream', action='store_true', help='Read VCF from stdin')
    parser.add_argument('--gff', type=str, required=True, help='GFF3 annotation file')
    parser.add_argument('--degeneracy', type=str, nargs='+', required=True,
                        help='Degeneracy annotation file(s) from describe_gene_positions_by_degeneracy.py')
    parser.add_argument('--chromosome', type=str, default=None,
                        help='Process only this chromosome (for parallel processing)')
    parser.add_argument('--output', type=str, default='pi_by_compartment.txt',
                        help='Output file (default: pi_by_compartment.txt)')
    
    args = parser.parse_args()
    
    if not args.vcf and not args.stream:
        print("Error: Must provide --vcf or --stream")
        sys.exit(1)
    
    target_chrom = args.chromosome
    if target_chrom:
        print(f"Filtering for chromosome: {target_chrom}")
    
    # 1. Parse GFF3 (filter by chromosome if specified)
    genes, exons, introns, gene_intervals = parse_gff3(args.gff, target_chrom=target_chrom)
    
    # 2. Calculate intergenic windows
    intergenic = calculate_intergenic_windows(gene_intervals)
    
    # 3. Load degeneracy annotations
    degeneracy = load_degeneracy_annotations(args.degeneracy)
    
    # 4. Process VCF
    stats = process_vcf(
        args.vcf if args.vcf else None,
        exons, introns, intergenic, degeneracy,
        stream=args.stream,
        target_chrom=target_chrom
    )
    
    # 5. Write output
    write_summary(stats, args.output, chromosome=target_chrom)
    
    print("\nDone!")


if __name__ == "__main__":
    main()
