#!/usr/bin/env python3
"""
Calculate nucleotide diversity (π) across different genomic compartments.

Compartments:
    1. intergenic - windows (50kb) far from genes - neutral baseline
    2. intergenic_upstream_2kb - 2kb upstream of genes (promoter region)
    3. intergenic_upstream_10kb - 10kb upstream of genes
    4. intron - introns (50bp trimmed from splice sites) - neutral within genes
    5. exon_all - ALL exonic sites (any degeneracy) - amino acid selection visible
    6. first_exon_4fold - first exons, 4-fold degenerate sites only
    7. nonfirst_exon_4fold - non-first exons, 4-fold degenerate sites only

Nucleotide categories (C, G, AT) are based on the REFERENCE allele only.
This ensures consistent site counting where each site belongs to exactly one
nucleotide category, and monomorphic sites are properly included.

Input:
    - VCF file with variant AND invariant sites
    - GFF3 annotation file
    - Degeneracy annotation file (from describe_gene_positions_by_degeneracy.py)

Output:
    - pi_by_compartment.txt: Summary statistics per compartment (Updated with Additive Components)
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

# Coordinate system notes:
# - VCF: 1-based positions
# - GFF3: 1-based, fully-closed intervals [start, end]
# - Internal storage: 0-based, half-open [start, end) for efficient binary search
# - Degeneracy annotations: 1-based (from describe_gene_positions_by_degeneracy.py)
#
# The binary_search_region function expects 0-based half-open intervals,
# while degeneracy lookups use 1-based positions directly.

INTRON_TRIM_BP = 50       # bp to trim from intron boundaries (splice site removal)
INTERGENIC_WINDOW = 50000 # 50kb windows for intergenic regions
UPSTREAM_2KB = 2000       # 2kb upstream of gene transcription start
UPSTREAM_10KB = 10000     # 10kb upstream of gene transcription start
MIN_SAMPLES = 10          # Minimum sample size for π calculation
MIN_DEPTH_RATIO = 5       # Depth ratio filter for homozygous calls (matches proc2.py)

# Nucleotide categories for polarized analysis
# Based on mutation bias literature (Monroe et al. 2022 Nature; Bird 1980)
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


def get_site_nucleotide_category(ref, strand):
    """
    Determine which nucleotide category a site belongs to based on REFERENCE allele only.
    
    Using only the reference allele ensures:
    1. Each site belongs to exactly ONE nucleotide category
    2. Monomorphic sites are correctly counted in their category
    3. sum(C + G + AT sites) = all sites
    4. π calculations are consistent across categories
    
    Biological rationale for C/G/AT categorization:
    - C sites: Subject to CpG methylation → deamination (C→T), leading to elevated
      mutation rates at methylated cytosines (Bird 1980; Monroe et al. 2022 Nature)
    - G sites: Complementary to C; also shows elevated diversity due to C→T on
      opposite strand appearing as G→A
    - AT sites: Lower mutation rate baseline; used as neutral reference for
      comparing C/G diversity (Boman et al. 2021 GBE)
    
    For coding regions on minus strand, we apply strand correction so that
    categories reflect the sense strand nucleotide.
    
    Args:
        ref: Reference allele from VCF
        strand: '+' or '-' for strand correction (use '+' for intergenic)
    
    Returns: Single category string: 'C', 'G', or 'AT'
    """
    # Apply strand correction for coding regions
    if strand == '-':
        ref_corrected = get_complement(ref)
    else:
        ref_corrected = ref.upper()
    
    return classify_nucleotide(ref_corrected)


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


def calculate_upstream_regions(genes, gene_intervals):
    """
    Calculate upstream regions for each gene.
    
    Upstream is defined as 5' of the transcription start site (TSS):
    - For + strand genes: upstream is BEFORE (lower coordinates) the gene start
    - For - strand genes: upstream is AFTER (higher coordinates) the gene end
    
    We avoid overlapping with other genes by clipping upstream regions.
    
    Returns:
        upstream_2kb: {chrom: [(start, end, gene_id), ...]} for 2kb upstream
        upstream_10kb: {chrom: [(start, end, gene_id), ...]} for 10kb upstream
    """
    upstream_2kb = defaultdict(list)
    upstream_10kb = defaultdict(list)
    
    # Build merged gene intervals to check for overlaps
    merged_genes = {}
    for chrom, intervals in gene_intervals.items():
        # Merge overlapping gene intervals
        sorted_intervals = sorted(intervals)
        merged = []
        for start, end in sorted_intervals:
            if merged and start <= merged[-1][1]:
                # Overlapping - extend
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append([start, end])
        merged_genes[chrom] = [(s, e) for s, e in merged]
    
    def find_upstream_boundary(chrom, pos, direction):
        """Find nearest gene boundary in the given direction."""
        if chrom not in merged_genes:
            return None
        for start, end in merged_genes[chrom]:
            if direction == 'left':  # Looking for boundary to the left of pos
                if end <= pos:
                    continue  # Keep looking
                if start >= pos:
                    return None  # No boundary found before hitting another gene
            elif direction == 'right':  # Looking for boundary to the right of pos
                if start >= pos:
                    return start  # Found next gene start
        return None
    
    for gene_id, info in genes.items():
        chrom = info['chrom']
        strand = info['strand']
        gene_start = info['start']
        gene_end = info['end']
        
        if strand == '+':
            # Upstream is before gene start (lower coordinates)
            up_end = gene_start
            up_start_2kb = max(0, gene_start - UPSTREAM_2KB)
            up_start_10kb = max(0, gene_start - UPSTREAM_10KB)
            
            # Check for overlapping genes - find the nearest gene end before this gene
            for g_start, g_end in merged_genes.get(chrom, []):
                if g_end <= gene_start and g_end > up_start_10kb:
                    # There's a gene ending in our upstream region
                    up_start_10kb = max(up_start_10kb, g_end)
                if g_end <= gene_start and g_end > up_start_2kb:
                    up_start_2kb = max(up_start_2kb, g_end)
            
            if up_start_2kb < up_end:
                upstream_2kb[chrom].append((up_start_2kb, up_end, gene_id))
            if up_start_10kb < up_end:
                upstream_10kb[chrom].append((up_start_10kb, up_end, gene_id))
                
        else:  # strand == '-'
            # Upstream is after gene end (higher coordinates)
            up_start = gene_end
            up_end_2kb = gene_end + UPSTREAM_2KB
            up_end_10kb = gene_end + UPSTREAM_10KB
            
            # Check for overlapping genes - find the nearest gene start after this gene
            for g_start, g_end in merged_genes.get(chrom, []):
                if g_start >= gene_end and g_start < up_end_10kb:
                    # There's a gene starting in our upstream region
                    up_end_10kb = min(up_end_10kb, g_start)
                if g_start >= gene_end and g_start < up_end_2kb:
                    up_end_2kb = min(up_end_2kb, g_start)
            
            if up_start < up_end_2kb:
                upstream_2kb[chrom].append((up_start, up_end_2kb, gene_id))
            if up_start < up_end_10kb:
                upstream_10kb[chrom].append((up_start, up_end_10kb, gene_id))
    
    # Sort for binary search
    for chrom in upstream_2kb:
        upstream_2kb[chrom].sort()
    for chrom in upstream_10kb:
        upstream_10kb[chrom].sort()
    
    total_2kb = sum(len(v) for v in upstream_2kb.values())
    total_10kb = sum(len(v) for v in upstream_10kb.values())
    print(f"  Upstream regions: {total_2kb} (2kb), {total_10kb} (10kb)")
    
    return upstream_2kb, upstream_10kb


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


def classify_position(chrom, pos, exons, introns, intergenic, degeneracy, 
                       upstream_2kb=None, upstream_10kb=None):
    """
    Classify a position into one of the compartments.
    
    Priority order (higher priority compartments checked first):
    1. Exons (returns exon_all, plus first_exon_4fold or nonfirst_exon_4fold if 4-fold)
    2. Introns
    3. Upstream (2kb, then 10kb - these overlap, so a site can be in both)
    4. Intergenic (far from genes)
    
    Returns:
        (compartments, details)
        compartments: list of matching compartments (a site can match multiple)
        details: dict with additional info including strand
    """
    compartments = []
    details = {'strand': '+'}  # Default strand for intergenic
    
    # Check exons first
    exon_match = binary_search_region(chrom, pos, exons)
    if exon_match:
        ex_start, ex_end, gene_id, exon_num, strand = exon_match
        details = {'gene': gene_id, 'exon': exon_num, 'strand': strand}
        
        # Always add to exon_all
        compartments.append('exon_all')
        
        # Check if 4-fold degenerate
        if chrom in degeneracy and pos in degeneracy[chrom]:
            deg_info = degeneracy[chrom][pos]
            if deg_info['degeneracy'] == '4-fold':
                if exon_num == 1:
                    compartments.append('first_exon_4fold')
                else:
                    compartments.append('nonfirst_exon_4fold')
        
        return compartments, details
    
    # Check introns
    intron_match = binary_search_region(chrom, pos, introns)
    if intron_match:
        in_start, in_end, gene_id, strand = intron_match
        details = {'gene': gene_id, 'strand': strand}
        return ['intron'], details
    
    # Check upstream regions (can be in both 2kb and 10kb)
    # 2kb is a subset of 10kb, but we track them separately
    if upstream_2kb:
        up2_match = binary_search_region(chrom, pos, upstream_2kb)
        if up2_match:
            up_start, up_end, gene_id = up2_match
            details = {'gene': gene_id, 'strand': '+'}  # Upstream uses + for consistent C/G
            compartments.append('intergenic_upstream_2kb')
    
    if upstream_10kb:
        up10_match = binary_search_region(chrom, pos, upstream_10kb)
        if up10_match:
            up_start, up_end, gene_id = up10_match
            details = {'gene': gene_id, 'strand': '+'}
            if 'intergenic_upstream_10kb' not in compartments:
                compartments.append('intergenic_upstream_10kb')
    
    if compartments:
        return compartments, details
    
    # Check intergenic (far from genes)
    ig_match = binary_search_region(chrom, pos, intergenic)
    if ig_match:
        ig_start, ig_end = ig_match[:2]
        details = {'window': f"{chrom}:{ig_start}-{ig_end}", 'strand': '+'}
        return ['intergenic'], details
    
    return [], {}


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
    pos = int(parts[1])  # Keep 1-based to match degeneracy annotations
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
    
    This implements the standard Nei (1987) estimator for nucleotide diversity:
        π = 2 * n * p * (1-p) / (n-1)
    
    where:
        n = sample size (number of clear homozygous calls)
        p = frequency of reference allele
    
    This is mathematically equivalent to pixy's approach (Korunes & Samuk 2021):
        π = count_diffs / count_comps
    
    For inbred lines, we only use homozygous calls with clear depth support
    (depth ratio > 5:1) to avoid heterozygous contamination artifacts.
    
    References:
        - Nei, M. (1987). Molecular Evolutionary Genetics. Columbia University Press.
        - Korunes & Samuk (2021). Mol Ecol Res. https://doi.org/10.1111/1755-0998.13326
    
    Returns:
        (is_polymorphic, pi_value, n_samples)
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
    compartments = [
        'intergenic',           # Far from genes (50kb windows)
        'intergenic_upstream_2kb',   # 2kb upstream of TSS
        'intergenic_upstream_10kb',  # 10kb upstream of TSS
        'intron',               # Introns (trimmed)
        'exon_all',             # All exonic sites (any degeneracy)
        'first_exon_4fold',     # First exon, 4-fold sites
        'nonfirst_exon_4fold'   # Non-first exons, 4-fold sites
    ]
    stats = {}
    for comp in compartments:
        stats[comp] = {}
        for nuc in NUC_CATEGORIES:
            stats[comp][nuc] = {'sites': 0, 'poly': 0, 'pi_sum': 0.0}
    return stats


def process_vcf(vcf_path, exons, introns, intergenic, degeneracy, 
                 upstream_2kb=None, upstream_10kb=None,
                 stream=False, target_chrom=None):
    """
    Process VCF and calculate π for each compartment, with C/G polarization.
    
    Args:
        vcf_path: Path to VCF file
        exons, introns, intergenic, degeneracy: Annotation data structures
        upstream_2kb, upstream_10kb: Upstream region data structures
        stream: If True, read from stdin
        target_chrom: If set, only process this chromosome (for parallel processing)
    
    Returns:
        compartment_stats: {compartment: {nuc_category: {'sites': int, 'poly': int, 'pi_sum': float}}}
        
    Where nuc_category is 'all', 'C', 'G', or 'AT'
    
    Note: Nucleotide categories are based on REFERENCE allele only.
    This ensures each site belongs to exactly one nucleotide category,
    and sum(C + G + AT sites) = all sites.
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
        
        # Classify position - now returns list of compartments
        compartments, details = classify_position(
            chrom, pos, exons, introns, intergenic, degeneracy,
            upstream_2kb, upstream_10kb
        )
        
        if not compartments:
            continue
        
        classified_count += 1
        
        # Get strand from details (for exon strand correction)
        strand = details.get('strand', '+')
        
        # Get nucleotide category based on REFERENCE allele only
        # This ensures each site belongs to exactly ONE category
        nuc_category = get_site_nucleotide_category(ref, strand)
        
        # Calculate pi for polymorphic sites
        is_poly = False
        pi_val = 0.0
        if not is_invariant:
            is_poly, pi_val, _ = calculate_pi_site(genotypes)
        
        # Update stats for ALL compartments this site belongs to
        for compartment in compartments:
            # Update 'all' category
            stats[compartment]['all']['sites'] += 1
            if is_poly:
                stats[compartment]['all']['poly'] += 1
                stats[compartment]['all']['pi_sum'] += pi_val
            
            # Update specific nucleotide category
            if nuc_category:
                stats[compartment][nuc_category]['sites'] += 1
                if is_poly:
                    stats[compartment][nuc_category]['poly'] += 1
                    stats[compartment][nuc_category]['pi_sum'] += pi_val
    
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
    
    # All compartments in logical order
    compartments = [
        'intergenic',
        'intergenic_upstream_10kb',
        'intergenic_upstream_2kb',
        'intron',
        'exon_all',
        'first_exon_4fold',
        'nonfirst_exon_4fold'
    ]
    
    with open(output_file, 'w') as out:
        # Header with chromosome column for merging across parallel jobs
        # [MODIFIED]: Added 'Pi_component' to the header
        out.write("Chromosome\tCompartment\tNuc_Category\tSites\tPolymorphic\tPi_sum\tPi_mean\tPi_component\tPoly_fraction\n")
        
        for compartment in compartments:
            
            if compartment not in stats:
                continue

            # [ADDED]: Get total sites for the compartment to use as denominator for additive components
            # This follows the decomposition framework where sum(Pi_component_X) = Pi_total
            total_compartment_sites = stats[compartment]['all']['sites']

            for nuc_cat in NUC_CATEGORIES:
                s = stats[compartment][nuc_cat]
                n_sites = s['sites']
                n_poly = s['poly']
                pi_sum = s['pi_sum']
                
                # 1. Conditional Mean: Diversity *within* this category (Sum_Pi / N_category)
                pi_mean = pi_sum / n_sites if n_sites > 0 else 0.0
                
                # [ADDED] 2. Additive Component: Contribution to total diversity (Sum_Pi / N_total_compartment)
                # This ensures that Pi_component(C) + Pi_component(G) + Pi_component(AT) = Pi_total
                pi_component = pi_sum / total_compartment_sites if total_compartment_sites > 0 else 0.0
                
                poly_frac = n_poly / n_sites if n_sites > 0 else 0.0
                
                # [MODIFIED]: Added pi_component to output
                out.write(f"{chrom_str}\t{compartment}\t{nuc_cat}\t{n_sites}\t{n_poly}\t{pi_sum:.6f}\t{pi_mean:.8f}\t{pi_component:.8f}\t{poly_frac:.6f}\n")
    
    # Also print to console - formatted nicely
    print("\n" + "=" * 95)
    print(f"SUMMARY: Nucleotide Diversity (π) by Genomic Compartment - {chrom_str}")
    print("Nucleotide categories based on REFERENCE allele (strand-corrected for coding regions)")
    print("=" * 95)
    
    # Print overall stats first
    print("\n--- OVERALL (all nucleotides) ---")
    print(f"{'Compartment':<30} {'Sites':>12} {'Poly':>10} {'π mean':>12}")
    print("-" * 65)
    for compartment in compartments:
        if compartment not in stats:
            continue
        s = stats[compartment]['all']
        n_sites = s['sites']
        n_poly = s['poly']
        if n_sites == 0:
            continue
        pi_mean = s['pi_sum'] / n_sites if n_sites > 0 else 0.0
        print(f"{compartment:<30} {n_sites:>12,} {n_poly:>10,} {pi_mean:>12.8f}")
    
    # Print C vs G vs AT comparison
    print("\n--- BY NUCLEOTIDE CATEGORY (Sites / Poly / π_mean) ---")
    print(f"{'Compartment':<22} {'Sites_C':>10} {'Poly_C':>8} {'π_C':>10} {'Sites_G':>10} {'Poly_G':>8} {'π_G':>10} {'Sites_AT':>10} {'Poly_AT':>8} {'π_AT':>10} {'C/AT':>6}")
    print("-" * 130)
    for compartment in compartments:
        if compartment not in stats:
            continue
        if stats[compartment]['all']['sites'] == 0:
            continue
        n_C = stats[compartment]['C']['sites']
        n_G = stats[compartment]['G']['sites']
        n_AT = stats[compartment]['AT']['sites']
        p_C = stats[compartment]['C']['poly']
        p_G = stats[compartment]['G']['poly']
        p_AT = stats[compartment]['AT']['poly']
        pi_C = stats[compartment]['C']['pi_sum'] / n_C if n_C > 0 else 0.0
        pi_G = stats[compartment]['G']['pi_sum'] / n_G if n_G > 0 else 0.0
        pi_AT = stats[compartment]['AT']['pi_sum'] / n_AT if n_AT > 0 else 0.0
        ratio_C_AT = pi_C / pi_AT if pi_AT > 0 else float('nan')
        print(f"{compartment:<22} {n_C:>10,} {p_C:>8,} {pi_C:>10.6f} {n_G:>10,} {p_G:>8,} {pi_G:>10.6f} {n_AT:>10,} {p_AT:>8,} {pi_AT:>10.6f} {ratio_C_AT:>6.2f}")
    
    # [ADDED]: Console verification of the additive components
    print("\n--- ADDITIVE COMPONENTS VERIFICATION ---")
    print("Checking decomposition: π_total ≈ π_comp(C) + π_comp(G) + π_comp(AT)")
    print(f"{'Compartment':<22} {'π_total':>10} {'Sum_Comps':>10} {'π_comp_C':>10} {'π_comp_G':>10} {'π_comp_AT':>10} {'Match':>6}")
    print("-" * 85)
    for compartment in compartments:
        if compartment not in stats or stats[compartment]['all']['sites'] == 0:
            continue
        
        n_all = stats[compartment]['all']['sites']
        pi_all = stats[compartment]['all']['pi_sum'] / n_all if n_all > 0 else 0.0
        
        # Calculate components using total site count as denominator
        comp_C = stats[compartment]['C']['pi_sum'] / n_all
        comp_G = stats[compartment]['G']['pi_sum'] / n_all
        comp_AT = stats[compartment]['AT']['pi_sum'] / n_all
        
        sum_comps = comp_C + comp_G + comp_AT
        
        match = "✓" if abs(pi_all - sum_comps) < 1e-10 else "✗"
        print(f"{compartment:<22} {pi_all:>10.6f} {sum_comps:>10.6f} {comp_C:>10.6f} {comp_G:>10.6f} {comp_AT:>10.6f} {match:>6}")
    
    print("=" * 95)
    print("\nNotes:")
    print("  - Nucleotide categories use REFERENCE allele only (each site in exactly one category)")
    print("  - C and G are strand-corrected for coding regions")
    print("  - AT = sites with A or T reference (weak nucleotides)")
    print("  - π_mean: Average diversity within the category (can be > π_total for C sites)")
    print("  - π_comp: Additive contribution to total diversity (π_comp values sum to π_total)")
    print("  - exon_all includes ALL exonic sites (selection on amino acids visible)")
    print("  - first/nonfirst_exon_4fold are subsets of exon_all (synonymous sites only)")


# ====================== MAIN ======================

def main():
    parser = argparse.ArgumentParser(
        description="Calculate π across genomic compartments with C/G/AT polarization"
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
    
    print("\n=== Step 1: Parse GFF3 annotations ===")
    genes, exons, introns, gene_intervals = parse_gff3(args.gff, target_chrom=target_chrom)
    
    print("\n=== Step 2: Calculate intergenic windows ===")
    intergenic = calculate_intergenic_windows(gene_intervals)
    
    print("\n=== Step 3: Calculate upstream regions ===")
    upstream_2kb, upstream_10kb = calculate_upstream_regions(genes, gene_intervals)
    
    print("\n=== Step 4: Load degeneracy annotations ===")
    degeneracy = load_degeneracy_annotations(args.degeneracy)
    
    print("\n=== Step 5: Process VCF ===")
    stats = process_vcf(
        args.vcf if args.vcf else None,
        exons, introns, intergenic, degeneracy,
        upstream_2kb=upstream_2kb,
        upstream_10kb=upstream_10kb,
        stream=args.stream,
        target_chrom=target_chrom
    )
    
    print("\n=== Step 6: Write output ===")
    write_summary(stats, args.output, chromosome=target_chrom)
    
    print("\nDone!")

if __name__ == "__main__":
    main()