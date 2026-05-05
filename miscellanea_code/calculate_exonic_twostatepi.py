#!/usr/bin/env python3
"""
Calculate two-allele π (nucleotide diversity) statistics at exonic (CDS) sites.
Mirrors calculate_intronic_twostatepi.py but for exonic regions.

Usage:
    python3 calculate_exonic_twostatepi.py <vcf> <gff3> <out.csv> [options]

Example:
    zcat data/variants.vcf.gz | python3 calculate_exonic_twostatepi.py /dev/stdin data/annotation.gff3 results/exonic_pi.csv
"""

import sys
import argparse
import csv
from collections import defaultdict

def parse_gff3_exons(gff3_file, buffer_bytes=16 * 1024 * 1024):
    """
    Parse GFF3 file and extract CDS (exonic) regions.
    Returns dict: exon_id -> {'chrom': ..., 'start': ..., 'end': ..., 'strand': ..., 'gene_id': ...}
    """
    exons = {}
    gene_to_strand = {}
    chrom_set = set()
    
    try:
        if gff3_file == '/dev/stdin' or gff3_file == '-':
            f = sys.stdin
        else:
            f = open(gff3_file, 'r', buffering=buffer_bytes)
    except Exception as e:
        sys.stderr.write(f"Error opening GFF3 file: {e}\n")
        sys.exit(1)
    
    exon_count = 0
    try:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            
            chrom, source, feature, start, end, score, strand, phase, attributes = fields[:9]
            
            # Only consider CDS features as exons
            if feature != 'CDS':
                continue
            
            start, end = int(start), int(end)
            chrom_set.add(chrom)
            
            # Parse attributes to extract gene_id
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, val = attr.split('=', 1)
                    attr_dict[key.strip()] = val.strip()
            
            parent = attr_dict.get('Parent', '')
            if not parent:
                continue
            
            # Normalize gene ID: strip 'MgIM767.' prefix and '.v2.1' suffix if present
            gene_id = parent
            if gene_id.startswith('MgIM767.'):
                gene_id = gene_id[8:]
            if gene_id.endswith('.v2.1'):
                gene_id = gene_id[:-5]
            
            # Record strand for normalization
            gene_to_strand[gene_id] = strand
            
            exon_id = f"{chrom}:{start}-{end}:{gene_id}"
            exons[exon_id] = {
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'gene_id': gene_id
            }
            exon_count += 1
    finally:
        if gff3_file != '/dev/stdin' and gff3_file != '-':
            f.close()
    
    sys.stderr.write(f"Parsed {len(chrom_set)} chromosomes with {exon_count} CDS features\n")
    return exons, gene_to_strand


def complement_base(base):
    """Return complement of a DNA base."""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return comp.get(base, 'N')


def reverse_complement(seq):
    """Reverse complement a DNA sequence."""
    return ''.join(complement_base(b) for b in reversed(seq))


def parse_vcf_line(line):
    """Parse a single VCF line. Returns (chrom, pos, ref, alt, genotypes, ads) or None."""
    fields = line.strip().split('\t')
    if len(fields) < 10:
        return None
    
    chrom, pos, vid, ref, alt, qual, filt, info, fmt = fields[:9]
    samples = fields[9:]
    
    try:
        pos = int(pos)
    except:
        return None
    
    # Skip multiallelic or missing ALT
    if ',' in alt or alt == '.':
        return (chrom, pos, ref, alt, [], [])
    
    # Parse FORMAT and samples
    fmt_keys = fmt.split(':')
    gt_idx = fmt_keys.index('GT') if 'GT' in fmt_keys else -1
    ad_idx = fmt_keys.index('AD') if 'AD' in fmt_keys else -1
    
    genotypes = []
    ads = []
    for sample in samples:
        vals = sample.split(':')
        gt = vals[gt_idx] if gt_idx >= 0 and gt_idx < len(vals) else './.'
        ad = vals[ad_idx] if ad_idx >= 0 and ad_idx < len(vals) else '.'
        genotypes.append(gt)
        ads.append(ad)
    
    return (chrom, pos, ref, alt, genotypes, ads)


def calc_hom_counts(genotypes, min_ratio=5):
    """
    Count reference and alternate homozygotes from genotypes.
    Returns (ref_hom_count, alt_hom_count).
    """
    ref_hom = 0
    alt_hom = 0
    
    for gt in genotypes:
        if gt == '0/0' or gt == '0|0':
            ref_hom += 1
        elif gt == '1/1' or gt == '1|1':
            alt_hom += 1
    
    return ref_hom, alt_hom


def process_vcf(vcf_file, exons, gene_to_strand, buffer_bytes=16 * 1024 * 1024):
    """
    Process VCF file and compute π for each exon.
    Returns dict: gene_id -> {'exon_count': ..., 'total_sites': ..., 'q_C': ..., 'pi_C': ..., 'q_G': ..., 'pi_G': ..., 'exon_lengths': [...]}
    """
    # Build interval tree for fast exon lookup
    exon_intervals = defaultdict(list)
    for exon_id, exon_info in exons.items():
        chrom = exon_info['chrom']
        start = exon_info['start']
        end = exon_info['end']
        exon_intervals[chrom].append((start, end, exon_id))
    
    # Sort intervals for faster lookup
    for chrom in exon_intervals:
        exon_intervals[chrom].sort()
    
    # Accumulate per-gene stats
    gene_stats = defaultdict(lambda: {
        'sites': 0,
        'q_C_sum': 0.0,
        'pi_C_sum': 0.0,
        'q_G_sum': 0.0,
        'pi_G_sum': 0.0,
        'exon_lengths': []
    })
    
    try:
        if vcf_file == '/dev/stdin' or vcf_file == '-':
            f = sys.stdin
        else:
            f = open(vcf_file, 'r', buffering=buffer_bytes)
    except Exception as e:
        sys.stderr.write(f"Error opening VCF file: {e}\n")
        sys.exit(1)
    
    sites_processed = 0
    sites_in_exons = 0
    
    try:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parsed = parse_vcf_line(line)
            if not parsed:
                continue
            
            chrom, pos, ref, alt, genotypes, ads = parsed
            sites_processed += 1
            
            # Find overlapping exons
            if chrom not in exon_intervals:
                continue
            
            overlapping = []
            for start, end, exon_id in exon_intervals[chrom]:
                if start <= pos <= end:
                    overlapping.append(exon_id)
            
            if not overlapping:
                continue
            
            sites_in_exons += 1
            
            # Process for each overlapping exon
            for exon_id in overlapping:
                exon_info = exons[exon_id]
                gene_id = exon_info['gene_id']
                strand = exon_info['strand']
                
                # Normalize alleles if on minus strand
                if strand == '-':
                    ref = reverse_complement(ref)
                    if alt != '.':
                        alt = reverse_complement(alt)
                
                # Count hom genotypes
                ref_hom, alt_hom = calc_hom_counts(genotypes)
                nx = ref_hom + alt_hom
                
                if nx == 0:
                    continue
                
                # Compute π for each nucleotide system
                # C system: ALT is C
                if alt == 'C':
                    if nx == ref_hom:
                        q_C = 0.0
                    elif nx == alt_hom:
                        q_C = 1.0
                    else:
                        q_C = alt_hom / nx
                    pi_C = 2.0 * ref_hom * alt_hom / (nx * (nx - 1)) if nx > 1 else 0.0
                else:
                    q_C = 0.0
                    pi_C = 0.0
                
                # G system: ALT is G
                if alt == 'G':
                    if nx == ref_hom:
                        q_G = 0.0
                    elif nx == alt_hom:
                        q_G = 1.0
                    else:
                        q_G = alt_hom / nx
                    pi_G = 2.0 * ref_hom * alt_hom / (nx * (nx - 1)) if nx > 1 else 0.0
                else:
                    q_G = 0.0
                    pi_G = 0.0
                
                # Accumulate
                gene_stats[gene_id]['sites'] += 1
                gene_stats[gene_id]['q_C_sum'] += q_C
                gene_stats[gene_id]['pi_C_sum'] += pi_C
                gene_stats[gene_id]['q_G_sum'] += q_G
                gene_stats[gene_id]['pi_G_sum'] += pi_G
                
                # Track exon length (only once per exon per gene)
                if len(gene_stats[gene_id]['exon_lengths']) == 0:
                    exon_len = exon_info['end'] - exon_info['start'] + 1
                    gene_stats[gene_id]['exon_lengths'].append(exon_len)
    
    finally:
        if vcf_file != '/dev/stdin' and vcf_file != '-':
            f.close()
    
    sys.stderr.write(f"Processed {sites_processed} sites; {sites_in_exons} fell in exons\n")
    
    return gene_stats


def write_output(gene_stats, output_path, buffer_bytes=16 * 1024 * 1024):
    """Write per-gene exonic π statistics to CSV."""
    try:
        f = open(output_path, 'w', buffering=buffer_bytes, newline='')
    except Exception as e:
        sys.stderr.write(f"Error opening output file: {e}\n")
        sys.exit(1)
    
    writer = csv.writer(f)
    writer.writerow(['Gene', 'n_sites', 'q_pref_C', 'pi_2allele_C', 'q_pref_G', 'pi_2allele_G', 'exon_length_bp'])
    
    for gene_id in sorted(gene_stats.keys()):
        stats = gene_stats[gene_id]
        n_sites = stats['sites']
        
        if n_sites == 0:
            continue
        
        q_C = stats['q_C_sum'] / n_sites
        pi_C = stats['pi_C_sum'] / n_sites
        q_G = stats['q_G_sum'] / n_sites
        pi_G = stats['pi_G_sum'] / n_sites
        
        # For exonic regions, report total exonic length
        total_exon_len = sum(stats['exon_lengths']) if stats['exon_lengths'] else 0
        
        writer.writerow([
            gene_id,
            n_sites,
            f"{q_C:.8f}",
            f"{pi_C:.8f}",
            f"{q_G:.8f}",
            f"{pi_G:.8f}",
            total_exon_len
        ])
    
    f.close()
    sys.stderr.write(f"Wrote {len(gene_stats)} genes to {output_path}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Calculate two-allele π at exonic (CDS) sites.'
    )
    parser.add_argument('vcf', help='VCF file (or /dev/stdin)')
    parser.add_argument('gff3', help='GFF3 file with CDS features')
    parser.add_argument('output', help='Output CSV file')
    parser.add_argument('--buffer-mb', type=int, default=16, help='Buffer size in MB (default: 16)')
    parser.add_argument('--workers', type=int, default=1, help='Number of workers (placeholder for multiprocessing)')
    parser.add_argument('--batch-size', type=int, default=20000, help='Batch size for processing')
    
    args = parser.parse_args()
    
    buffer_bytes = args.buffer_mb * 1024 * 1024
    
    sys.stderr.write(f"Parsing GFF3: {args.gff3}\n")
    exons, gene_to_strand = parse_gff3_exons(args.gff3, buffer_bytes=buffer_bytes)
    
    sys.stderr.write(f"Processing VCF: {args.vcf}\n")
    gene_stats = process_vcf(args.vcf, exons, gene_to_strand, buffer_bytes=buffer_bytes)
    
    sys.stderr.write(f"Writing output: {args.output}\n")
    write_output(gene_stats, args.output, buffer_bytes=buffer_bytes)
    
    sys.stderr.write("Done.\n")


if __name__ == '__main__':
    main()
