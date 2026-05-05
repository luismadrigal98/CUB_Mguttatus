#!/usr/bin/env python3
"""
Calculate nucleotide diversity (π) for exonic (CDS) regions across all sites.
Regular π (not two-state C/G systems).

Usage:
    python3 calculate_exonic_pi.py <vcf> <gff3> <out.csv> [options]

Example:
    zcat data/variants.vcf.gz | python3 calculate_exonic_pi.py /dev/stdin data/annotation.gff3 results/exonic_pi.csv
"""

import sys
import argparse
import csv
from collections import defaultdict

def parse_gff3_exons(gff3_file, buffer_bytes=16 * 1024 * 1024):
    """
    Parse GFF3 file and extract CDS (exonic) regions.
    Returns dict: exon_id -> {'chrom': ..., 'start': ..., 'end': ..., 'gene_id': ...}
    """
    exons = {}
    
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
            
            exon_id = f"{chrom}:{start}-{end}:{gene_id}"
            exons[exon_id] = {
                'chrom': chrom,
                'start': start,
                'end': end,
                'gene_id': gene_id
            }
            exon_count += 1
    finally:
        if gff3_file != '/dev/stdin' and gff3_file != '-':
            f.close()
    
    sys.stderr.write(f"Parsed {exon_count} CDS features\n")
    return exons


def parse_vcf_line(line):
    """Parse a single VCF line. Returns (chrom, pos, ref, alt, genotypes) or None."""
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
        return (chrom, pos, ref, alt, [])
    
    # Parse FORMAT and samples
    fmt_keys = fmt.split(':')
    gt_idx = fmt_keys.index('GT') if 'GT' in fmt_keys else -1
    
    genotypes = []
    for sample in samples:
        vals = sample.split(':')
        gt = vals[gt_idx] if gt_idx >= 0 and gt_idx < len(vals) else './.'
        genotypes.append(gt)
    
    return (chrom, pos, ref, alt, genotypes)


def calc_hom_counts(genotypes):
    """
    Count reference and alternate alleles from genotypes.
    Returns (n_ref, n_alt, n_total).
    """
    n_ref = 0
    n_alt = 0
    
    for gt in genotypes:
        if gt == '0/0' or gt == '0|0':
            n_ref += 2
        elif gt == '1/1' or gt == '1|1':
            n_alt += 2
        elif gt == '0/1' or gt == '1/0' or gt == '0|1' or gt == '1|0':
            n_ref += 1
            n_alt += 1
    
    return n_ref, n_alt, n_ref + n_alt


def calc_pi(n_ref, n_alt, n_total):
    """
    Calculate π (nucleotide diversity) for a site.
    π = 2 * n_ref * n_alt / (n_total * (n_total - 1))
    """
    if n_total < 2:
        return 0.0
    return 2.0 * n_ref * n_alt / (n_total * (n_total - 1))


def process_vcf(vcf_file, exons, buffer_bytes=16 * 1024 * 1024):
    """
    Process VCF file and compute π for each exon.
    Returns dict: exon_id -> {'sites': n_sites, 'pi_sum': sum_of_pi, 'pi_values': [list of per-site π]}
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
    
    # Accumulate per-exon stats
    exon_stats = defaultdict(lambda: {
        'sites': 0,
        'pi_sum': 0.0,
        'pi_values': []
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
            
            chrom, pos, ref, alt, genotypes = parsed
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
            
            # Compute π for this site
            n_ref, n_alt, n_total = calc_hom_counts(genotypes)
            pi_val = calc_pi(n_ref, n_alt, n_total)
            
            # Accumulate for each overlapping exon
            for exon_id in overlapping:
                exon_stats[exon_id]['sites'] += 1
                exon_stats[exon_id]['pi_sum'] += pi_val
                exon_stats[exon_id]['pi_values'].append(pi_val)
    
    finally:
        if vcf_file != '/dev/stdin' and vcf_file != '-':
            f.close()
    
    sys.stderr.write(f"Processed {sites_processed} sites; {sites_in_exons} fell in exons\n")
    
    return exon_stats


def write_output(exons, exon_stats, output_path, buffer_bytes=16 * 1024 * 1024):
    """Write per-exon π statistics to CSV."""
    try:
        f = open(output_path, 'w', buffering=buffer_bytes, newline='')
    except Exception as e:
        sys.stderr.write(f"Error opening output file: {e}\n")
        sys.exit(1)
    
    writer = csv.writer(f)
    writer.writerow(['exon_id', 'chrom', 'start', 'end', 'length_bp', 'gene_id', 'n_sites', 'pi'])
    
    for exon_id in sorted(exon_stats.keys()):
        stats = exon_stats[exon_id]
        n_sites = stats['sites']
        
        if n_sites == 0:
            continue
        
        exon_info = exons[exon_id]
        mean_pi = stats['pi_sum'] / n_sites
        
        writer.writerow([
            exon_id,
            exon_info['chrom'],
            exon_info['start'],
            exon_info['end'],
            exon_info['end'] - exon_info['start'] + 1,
            exon_info['gene_id'],
            n_sites,
            f"{mean_pi:.8f}"
        ])
    
    f.close()
    sys.stderr.write(f"Wrote {len(exon_stats)} exons to {output_path}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Calculate nucleotide diversity (π) for exonic regions.'
    )
    parser.add_argument('vcf', help='VCF file (or /dev/stdin)')
    parser.add_argument('gff3', help='GFF3 file with CDS features')
    parser.add_argument('output', help='Output CSV file')
    parser.add_argument('--buffer-mb', type=int, default=16, help='Buffer size in MB (default: 16)')
    
    args = parser.parse_args()
    
    buffer_bytes = args.buffer_mb * 1024 * 1024
    
    sys.stderr.write(f"Parsing GFF3: {args.gff3}\n")
    exons = parse_gff3_exons(args.gff3, buffer_bytes=buffer_bytes)
    
    sys.stderr.write(f"Processing VCF: {args.vcf}\n")
    exon_stats = process_vcf(args.vcf, exons, buffer_bytes=buffer_bytes)
    
    sys.stderr.write(f"Writing output: {args.output}\n")
    write_output(exons, exon_stats, args.output, buffer_bytes=buffer_bytes)
    
    sys.stderr.write("Done.\n")


if __name__ == '__main__':
    main()
