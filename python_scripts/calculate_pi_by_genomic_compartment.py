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

Nucleotide categories (C, G, AT, CG) are based on ALLELE PRESENCE.
Sites are classified based on the set of alleles observed (Ref and Alt):
    - CG: Both C and G alleles present (e.g., C/G SNP)
    - C:  C allele present, no G (e.g., C/T or C/A SNP, or monomorphic C)
    - G:  G allele present, no C (e.g., G/A or G/T SNP, or monomorphic G)
    - AT: Only A and/or T alleles present (e.g., A/T SNP, or monomorphic A or T)

Input:
    - VCF file with variant AND invariant sites
    - GFF3 annotation file
    - Degeneracy annotation file (from describe_gene_positions_by_degeneracy.py)

Output:
    - pi_by_compartment.txt: Summary statistics per compartment
    - pi_by_window.txt: Per-window/region statistics

Author: Luis Javier Madrigal-Roca & GitHub Copilot
Date: 2026-01-22
"""

import sys
import os
import re
import argparse
import gzip
from collections import defaultdict
import math

# ====================== CONSTANTS ======================

INTRON_TRIM_BP = 50       # bp to trim from intron boundaries
INTERGENIC_WINDOW = 50000 # 50kb windows for intergenic regions
UPSTREAM_2KB = 2000       # 2kb upstream of gene transcription start
UPSTREAM_10KB = 10000     # 10kb upstream of gene transcription start
MIN_SAMPLES = 10          # Minimum sample size for π calculation
MIN_DEPTH_RATIO = 5       # Depth ratio filter for homozygous calls

# Nucleotide categories for polarized analysis
# Updated to include CG for sites where both C and G alleles segregate
NUC_CATEGORIES = ['all', 'C', 'G', 'AT', 'CG']


# ====================== NUCLEOTIDE POLARIZATION ======================

def get_complement(nuc):
    """Return complement of a nucleotide."""
    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', '.': '.'}
    return mapping.get(nuc.upper(), 'N')


def get_site_nucleotide_category(ref, alt, strand):
    """
    Determine which nucleotide category a site belongs to based on PRESENCE of alleles.
    
    Logic defined by advisor:
    1. Identify all alleles at site (Ref + Alt)
    2. Correct for strand if necessary (coding regions)
    3. Categorize:
       - If C and G present -> CG
       - Else if C present  -> C
       - Else if G present  -> G
       - Else (only A/T)    -> AT
       
    Args:
        ref: Reference allele from VCF
        alt: Alternate allele from VCF ('.' or '<NON_REF>' if invariant)
        strand: '+' or '-' for strand correction
    
    Returns: Single category string: 'C', 'G', 'AT', or 'CG'
    """
    # 1. Collect alleles
    alleles = set()
    alleles.add(ref.upper())
    
    # Check if variant
    is_variant = not (alt == '.' or alt == '<NON_REF>' or alt == '*')
    if is_variant:
        # Note: Multi-allelics (comma separated) are filtered out in process_vcf
        # so we assume alt is a single base here
        alleles.add(alt.upper())
    
    # 2. Strand correction
    final_alleles = set()
    if strand == '-':
        for base in alleles:
            final_alleles.add(get_complement(base))
    else:
        final_alleles = alleles
        
    # 3. Classification
    has_c = 'C' in final_alleles
    has_g = 'G' in final_alleles
    
    if has_c and has_g:
        return 'CG'
    elif has_c:
        return 'C'
    elif has_g:
        return 'G'
    else:
        # Default to AT if neither C nor G are present (implies only A, T, or N)
        return 'AT'


# ====================== GFF PARSING ======================

def parse_gff3(gff_file, target_chrom=None):
    """
    Parse GFF3 to extract gene boundaries, exons, and calculate introns.
    Returns genes, exons, introns, gene_intervals.
    """
    chrom_filter = f" (filtering for {target_chrom})" if target_chrom else ""
    print(f"Parsing GFF3: {gff_file}{chrom_filter}")
    
    genes = {}
    mrna_exons = defaultdict(list)
    mrna_to_gene = {}
    gene_info = {}
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            
            chrom, source, feature, start, end, score, strand, phase, attributes = parts
            if "Chr" not in chrom: continue
            if target_chrom and chrom != target_chrom: continue
            
            start, end = int(start) - 1, int(end)
            
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
                if mrna_id and parent: mrna_to_gene[mrna_id] = parent
            
            elif feature in ['CDS', 'exon']:
                parent = attr_dict.get('Parent', '')
                if parent: mrna_exons[parent].append((start, end))
    
    print(f"  Found {len(genes)} genes, {len(mrna_exons)} mRNAs with CDS/exon features")
    
    exons = defaultdict(list)
    introns = defaultdict(list)
    gene_intervals = defaultdict(list)
    
    for gene_id, info in genes.items():
        gene_intervals[info['chrom']].append((info['start'], info['end']))
    
    for chrom in gene_intervals:
        gene_intervals[chrom].sort()
    
    for mrna_id, exon_list in mrna_exons.items():
        gene_id = mrna_to_gene.get(mrna_id)
        if not gene_id or gene_id not in gene_info: continue
        
        chrom = gene_info[gene_id]['chrom']
        strand = gene_info[gene_id]['strand']
        sorted_exons = sorted(exon_list, key=lambda x: x[0])
        
        if strand == '-':
            for i, (ex_start, ex_end) in enumerate(reversed(sorted_exons)):
                exons[chrom].append((ex_start, ex_end, gene_id, i + 1, strand))
        else:
            for i, (ex_start, ex_end) in enumerate(sorted_exons):
                exons[chrom].append((ex_start, ex_end, gene_id, i + 1, strand))
        
        if len(sorted_exons) >= 2:
            for i in range(len(sorted_exons) - 1):
                intron_start = sorted_exons[i][1] + INTRON_TRIM_BP
                intron_end = sorted_exons[i + 1][0] - INTRON_TRIM_BP
                if intron_end > intron_start:
                    introns[chrom].append((intron_start, intron_end, gene_id, strand))
    
    for chrom in exons: exons[chrom].sort()
    for chrom in introns: introns[chrom].sort()
    
    return genes, exons, introns, gene_intervals


def calculate_intergenic_windows(gene_intervals):
    intergenic_windows = defaultdict(list)
    for chrom, intervals in gene_intervals.items():
        if not intervals: continue
        sorted_intervals = sorted(intervals)
        prev_end = 0
        
        for gene_start, gene_end in sorted_intervals:
            if gene_start > prev_end:
                ig_len = gene_start - prev_end
                if ig_len < INTERGENIC_WINDOW:
                    intergenic_windows[chrom].append((prev_end, gene_start))
                else:
                    for ws in range(prev_end, gene_start, INTERGENIC_WINDOW):
                        intergenic_windows[chrom].append((ws, min(ws + INTERGENIC_WINDOW, gene_start)))
            prev_end = max(prev_end, gene_end)
            
    print(f"  Intergenic: {sum(len(v) for v in intergenic_windows.values())} windows")
    return intergenic_windows


def calculate_upstream_regions(genes, gene_intervals):
    upstream_2kb = defaultdict(list)
    upstream_10kb = defaultdict(list)
    merged_genes = {}
    
    for chrom, intervals in gene_intervals.items():
        sorted_intervals = sorted(intervals)
        merged = []
        for s, e in sorted_intervals:
            if merged and s <= merged[-1][1]: merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else: merged.append([s, e])
        merged_genes[chrom] = merged
    
    for gene_id, info in genes.items():
        chrom, strand, g_start, g_end = info['chrom'], info['strand'], info['start'], info['end']
        
        if strand == '+':
            up_end = g_start
            up_start_2kb = max(0, g_start - UPSTREAM_2KB)
            up_start_10kb = max(0, g_start - UPSTREAM_10KB)
            for m_s, m_e in merged_genes.get(chrom, []):
                if m_e <= g_start:
                    if m_e > up_start_10kb: up_start_10kb = m_e
                    if m_e > up_start_2kb: up_start_2kb = m_e
            if up_start_2kb < up_end: upstream_2kb[chrom].append((up_start_2kb, up_end, gene_id))
            if up_start_10kb < up_end: upstream_10kb[chrom].append((up_start_10kb, up_end, gene_id))
        else:
            up_start = g_end
            up_end_2kb = g_end + UPSTREAM_2KB
            up_end_10kb = g_end + UPSTREAM_10KB
            for m_s, m_e in merged_genes.get(chrom, []):
                if m_s >= g_end:
                    if m_s < up_end_10kb: up_end_10kb = m_s
                    if m_s < up_end_2kb: up_end_2kb = m_s
            if up_start < up_end_2kb: upstream_2kb[chrom].append((up_start, up_end_2kb, gene_id))
            if up_start < up_end_10kb: upstream_10kb[chrom].append((up_start, up_end_10kb, gene_id))

    for c in upstream_2kb: upstream_2kb[c].sort()
    for c in upstream_10kb: upstream_10kb[c].sort()
    return upstream_2kb, upstream_10kb


# ====================== DEGENERACY & LOOKUP ======================

def load_degeneracy_annotations(annotation_files):
    degeneracy = defaultdict(dict)
    for annot_file in annotation_files:
        print(f"Loading degeneracy: {annot_file}")
        with open(annot_file, 'r') as f:
            f.readline()
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 8: continue
                degeneracy[parts[0]][int(parts[2])] = {'gene': parts[1], 'degeneracy': parts[5]}
    print(f"  Loaded annotated positions")
    return degeneracy

def binary_search_region(chrom, pos, regions):
    if chrom not in regions: return None
    lst = regions[chrom]
    l, h = 0, len(lst) - 1
    while l <= h:
        mid = (l + h) // 2
        s, e = lst[mid][:2]
        if s <= pos < e: return lst[mid]
        elif pos < s: h = mid - 1
        else: l = mid + 1
    return None

def classify_position(chrom, pos, exons, introns, intergenic, degeneracy, u2kb, u10kb):
    compartments = []
    details = {'strand': '+'}
    
    ex_match = binary_search_region(chrom, pos, exons)
    if ex_match:
        ex_s, ex_e, gid, ex_num, strand = ex_match
        details = {'gene': gid, 'exon': ex_num, 'strand': strand}
        compartments.append('exon_all')
        if chrom in degeneracy and pos in degeneracy[chrom]:
            if degeneracy[chrom][pos]['degeneracy'] == '4-fold':
                compartments.append('first_exon_4fold' if ex_num == 1 else 'nonfirst_exon_4fold')
        return compartments, details

    in_match = binary_search_region(chrom, pos, introns)
    if in_match:
        return ['intron'], {'gene': in_match[2], 'strand': in_match[3]}
    
    if u2kb:
        u2 = binary_search_region(chrom, pos, u2kb)
        if u2:
            compartments.append('intergenic_upstream_2kb')
            details = {'gene': u2[2], 'strand': '+'} # Upstream assumes sense orientation of upstream gene
    if u10kb:
        u10 = binary_search_region(chrom, pos, u10kb)
        if u10:
            if 'intergenic_upstream_10kb' not in compartments: compartments.append('intergenic_upstream_10kb')
            details = {'gene': u10[2], 'strand': '+'}

    if compartments: return compartments, details

    if binary_search_region(chrom, pos, intergenic):
        return ['intergenic'], {'strand': '+'}
    
    return [], {}


# ====================== VCF PROCESSING ======================

def parse_vcf_line(line):
    parts = line.strip().split('\t')
    if len(parts) < 10: return None
    chrom, pos, ref, alt = parts[0], int(parts[1]), parts[3], parts[4]
    
    fmt = parts[8].split(':')
    try:
        gt_idx, ad_idx = fmt.index('GT'), fmt.index('AD')
    except ValueError: return None
    
    genotypes = []
    for s in parts[9:]:
        sp = s.split(':')
        if len(sp) <= max(gt_idx, ad_idx):
            genotypes.append(('./.', 0, 0))
            continue
        ad_str = sp[ad_idx]
        try:
            if ',' in ad_str:
                adp = ad_str.split(',')
                rc, ac = int(adp[0]), (int(adp[1]) if len(adp) > 1 else 0)
            else: rc, ac = (int(ad_str) if ad_str != '.' else 0), 0
        except ValueError: rc, ac = 0, 0
        genotypes.append((sp[gt_idx], rc, ac))
    
    return chrom, pos, ref, alt, genotypes

def calculate_pi_site(genotypes):
    ref_hom, alt_hom = 0, 0
    for gt, rc, ac in genotypes:
        if gt == '0/0' and rc > MIN_DEPTH_RATIO * ac: ref_hom += 1
        elif gt == '1/1' and ac > MIN_DEPTH_RATIO * rc: alt_hom += 1
    
    if min(ref_hom, alt_hom) > 0:
        n = float(ref_hom + alt_hom)
        p = float(ref_hom) / n
        return True, 2.0 * n * p * (1.0 - p) / (n - 1.0), int(n)
    return False, 0.0, 0

def init_stats():
    comps = ['intergenic', 'intergenic_upstream_2kb', 'intergenic_upstream_10kb', 
             'intron', 'exon_all', 'first_exon_4fold', 'nonfirst_exon_4fold']
    stats = {}
    for c in comps:
        stats[c] = {}
        for n in NUC_CATEGORIES: stats[c][n] = {'sites': 0, 'poly': 0, 'pi_sum': 0.0}
    return stats

def process_vcf(vcf_path, exons, introns, intergenic, degeneracy, u2kb, u10kb, stream, target_chrom):
    stats = init_stats()
    
    if stream: f = sys.stdin
    elif vcf_path.endswith('.gz'): f = gzip.open(vcf_path, 'rt')
    else: f = open(vcf_path, 'r')
    
    lc, cc = 0, 0
    print("Processing VCF...")
    
    for line in f:
        if line.startswith('#'): continue
        lc += 1
        if lc % 1000000 == 0: sys.stderr.write(f"\rProcessed {lc:,} sites...")
        
        parsed = parse_vcf_line(line)
        if not parsed: continue
        chrom, pos, ref, alt, gts = parsed
        
        if target_chrom and chrom != target_chrom: continue
        
        # Skip multi-allelic or indels for now to keep logic simple
        is_inv = (alt == '.' or alt == '<NON_REF>' or alt == '*')
        if not is_inv and (len(ref) > 1 or len(alt) > 1 or ',' in alt): continue
        
        comps, det = classify_position(chrom, pos, exons, introns, intergenic, degeneracy, u2kb, u10kb)
        if not comps: continue
        cc += 1
        
        strand = det.get('strand', '+')
        # [MODIFIED]: Now passing 'alt' to the category function
        nuc_cat = get_site_nucleotide_category(ref, alt, strand)
        
        is_poly, pi_val = False, 0.0
        if not is_inv:
            is_poly, pi_val, _ = calculate_pi_site(gts)
            
        for c in comps:
            # Update 'all'
            stats[c]['all']['sites'] += 1
            if is_poly:
                stats[c]['all']['poly'] += 1
                stats[c]['all']['pi_sum'] += pi_val
            
            # Update specific category (C, G, AT, or CG)
            if nuc_cat:
                stats[c][nuc_cat]['sites'] += 1
                if is_poly:
                    stats[c][nuc_cat]['poly'] += 1
                    stats[c][nuc_cat]['pi_sum'] += pi_val
                    
    if not stream and f != sys.stdin: f.close()
    print(f"\nTotal sites: {lc:,}, Classified: {cc:,}")
    return stats


# ====================== OUTPUT ======================

def write_summary(stats, output_file, chromosome=None):
    chrom_str = chromosome if chromosome else "all"
    print(f"\nWriting summary to: {output_file}")
    
    compartments = ['intergenic', 'intergenic_upstream_10kb', 'intergenic_upstream_2kb',
                    'intron', 'exon_all', 'first_exon_4fold', 'nonfirst_exon_4fold']
    
    with open(output_file, 'w') as out:
        out.write("Chromosome\tCompartment\tNuc_Category\tSites\tPolymorphic\tPi_sum\tPi_mean\tPi_component\tPoly_fraction\n")
        
        for comp in compartments:
            if comp not in stats: continue
            
            # Denominator for additive components
            total_sites = stats[comp]['all']['sites']
            
            for nuc in NUC_CATEGORIES:
                s = stats[comp][nuc]
                ns, np, pis = s['sites'], s['poly'], s['pi_sum']
                
                pi_mean = pis / ns if ns > 0 else 0.0
                pi_comp = pis / total_sites if total_sites > 0 else 0.0
                pf = np / ns if ns > 0 else 0.0
                
                out.write(f"{chrom_str}\t{comp}\t{nuc}\t{ns}\t{np}\t{pis:.6f}\t{pi_mean:.8f}\t{pi_comp:.8f}\t{pf:.6f}\n")

    # Console Verification
    print("\n" + "="*95)
    print(f"SUMMARY: Nucleotide Diversity (π) by Genomic Compartment - {chrom_str}")
    print("Categories based on ALLELE PRESENCE (Strand Corrected)")
    print("="*95)
    
    print("\n--- ADDITIVE COMPONENTS VERIFICATION ---")
    print("Checking decomposition: π_total ≈ π_comp(C) + π_comp(G) + π_comp(AT) + π_comp(CG)")
    print(f"{'Compartment':<25} {'π_total':>10} {'Sum_Comps':>10} {'π_C':>8} {'π_G':>8} {'π_AT':>8} {'π_CG':>8} {'Match':>6}")
    print("-" * 100)
    
    for comp in compartments:
        if comp not in stats or stats[comp]['all']['sites'] == 0: continue
        
        n_all = stats[comp]['all']['sites']
        pi_all = stats[comp]['all']['pi_sum'] / n_all
        
        # Calculate components
        comp_C = stats[comp]['C']['pi_sum'] / n_all
        comp_G = stats[comp]['G']['pi_sum'] / n_all
        comp_AT = stats[comp]['AT']['pi_sum'] / n_all
        comp_CG = stats[comp]['CG']['pi_sum'] / n_all
        
        sum_comps = comp_C + comp_G + comp_AT + comp_CG
        match = "✓" if abs(pi_all - sum_comps) < 1e-9 else "✗"
        
        print(f"{comp:<25} {pi_all:>10.6f} {sum_comps:>10.6f} {comp_C:>8.6f} {comp_G:>8.6f} {comp_AT:>8.6f} {comp_CG:>8.6f} {match:>6}")


# ====================== MAIN ======================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf', type=str)
    parser.add_argument('--stream', action='store_true')
    parser.add_argument('--gff', type=str, required=True)
    parser.add_argument('--degeneracy', type=str, nargs='+', required=True)
    parser.add_argument('--chromosome', type=str)
    parser.add_argument('--output', type=str, default='pi_by_compartment.txt')
    args = parser.parse_args()
    
    if not args.vcf and not args.stream: sys.exit("Error: --vcf or --stream required")
    
    genes, exons, introns, intervals = parse_gff3(args.gff, args.chromosome)
    intergenic = calculate_intergenic_windows(intervals)
    u2kb, u10kb = calculate_upstream_regions(genes, intervals)
    degeneracy = load_degeneracy_annotations(args.degeneracy)
    
    stats = process_vcf(args.vcf, exons, introns, intergenic, degeneracy, u2kb, u10kb, args.stream, args.chromosome)
    write_summary(stats, args.output, args.chromosome)

if __name__ == "__main__":
    main()