#!/usr/bin/env python3
"""
Calculate π at the individual feature level (per-intron, per-exon) with feature sizes.
Outputs: feature_id, feature_type, feature_size, q_C, pi_C, q_G, pi_G

Usage:
    python3 calculate_pi_by_feature.py <vcf> <gff3> <out.csv> [options]
"""

import sys
import argparse
import csv
from collections import defaultdict

def parse_gff3_features(gff3_file, buffer_bytes=16 * 1024 * 1024):
    """
    Parse GFF3 file and extract both introns (from CDS) and exons (CDS).
    Returns:
        features: dict feature_id -> {'type': 'exon'|'intron', 'chrom': ..., 'start': ..., 'end': ..., 'strand': ...}
        gene_to_strand: dict gene_id -> strand
    """
    features = {}
    gene_to_strand = {}
    cds_by_gene = defaultdict(list)  # To compute introns
    
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
            
            if feature != 'CDS':
                continue
            
            start, end = int(start), int(end)
            
            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, val = attr.split('=', 1)
                    attr_dict[key.strip()] = val.strip()
            
            parent = attr_dict.get('Parent', '')
            if not parent:
                continue
            
            # Normalize gene ID
            gene_id = parent
            if gene_id.startswith('MgIM767.'):
                gene_id = gene_id[8:]
            if gene_id.endswith('.v2.1'):
                gene_id = gene_id[:-5]
            
            gene_to_strand[gene_id] = strand
            
            # Add exon feature
            exon_id = f"{chrom}:{start}-{end}:{gene_id}:exon"
            features[exon_id] = {
                'type': 'exon',
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'gene_id': gene_id
            }
            exon_count += 1
            
            # Collect CDS for intron computation
            cds_by_gene[(chrom, gene_id)].append((start, end))
    
    finally:
        if gff3_file != '/dev/stdin' and gff3_file != '-':
            f.close()
    
    # Compute introns from CDS boundaries
    intron_count = 0
    for (chrom, gene_id), cds_list in cds_by_gene.items():
        if len(cds_list) < 2:
            continue
        
        cds_list.sort()
        strand = gene_to_strand.get((chrom, gene_id), '+')
        
        for i in range(len(cds_list) - 1):
            intron_start = cds_list[i][1] + 1
            intron_end = cds_list[i + 1][0] - 1
            
            if intron_start >= intron_end:
                continue
            
            intron_id = f"{chrom}:{intron_start}-{intron_end}:{gene_id}:intron"
            features[intron_id] = {
                'type': 'intron',
                'chrom': chrom,
                'start': intron_start,
                'end': intron_end,
                'strand': strand,
                'gene_id': gene_id
            }
            intron_count += 1
    
    sys.stderr.write(f"Extracted {exon_count} exons and {intron_count} introns\n")
    return features, gene_to_strand


def complement_base(base):
    """Return complement of a DNA base."""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return comp.get(base, 'N')


def reverse_complement(seq):
    """Reverse complement a DNA sequence."""
    return ''.join(complement_base(b) for b in reversed(seq))


def parse_vcf_line(line):
    """Parse a single VCF line."""
    fields = line.strip().split('\t')
    if len(fields) < 10:
        return None
    
    chrom, pos, vid, ref, alt, qual, filt, info, fmt = fields[:9]
    samples = fields[9:]
    
    try:
        pos = int(pos)
    except:
        return None
    
    if ',' in alt or alt == '.':
        return (chrom, pos, ref, alt, [], [])
    
    fmt_keys = fmt.split(':')
    gt_idx = fmt_keys.index('GT') if 'GT' in fmt_keys else -1
    ad_idx = fmt_keys.index('AD') if 'AD' in fmt_keys else -1
    
    genotypes = []
    for sample in samples:
        vals = sample.split(':')
        gt = vals[gt_idx] if gt_idx >= 0 and gt_idx < len(vals) else './.'
        genotypes.append(gt)
    
    return (chrom, pos, ref, alt, genotypes, [])


def calc_hom_counts(genotypes):
    """Count reference and alternate homozygotes."""
    ref_hom = 0
    alt_hom = 0
    
    for gt in genotypes:
        if gt == '0/0' or gt == '0|0':
            ref_hom += 1
        elif gt == '1/1' or gt == '1|1':
            alt_hom += 1
    
    return ref_hom, alt_hom


def process_vcf(vcf_file, features, buffer_bytes=16 * 1024 * 1024):
    """Process VCF and compute π for each feature."""
    # Build interval tree
    feature_intervals = defaultdict(list)
    for feature_id, feature_info in features.items():
        chrom = feature_info['chrom']
        start = feature_info['start']
        end = feature_info['end']
        feature_intervals[chrom].append((start, end, feature_id))
    
    for chrom in feature_intervals:
        feature_intervals[chrom].sort()
    
    # Accumulate per-feature stats
    feature_stats = defaultdict(lambda: {
        'sites': 0,
        'q_C_sum': 0.0,
        'pi_C_sum': 0.0,
        'q_G_sum': 0.0,
        'pi_G_sum': 0.0
    })
    
    try:
        if vcf_file == '/dev/stdin' or vcf_file == '-':
            f = sys.stdin
        else:
            f = open(vcf_file, 'r', buffering=buffer_bytes)
    except Exception as e:
        sys.stderr.write(f"Error opening VCF: {e}\n")
        sys.exit(1)
    
    sites_processed = 0
    sites_in_features = 0
    
    try:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parsed = parse_vcf_line(line)
            if not parsed:
                continue
            
            chrom, pos, ref, alt, genotypes, _ = parsed
            sites_processed += 1
            
            if chrom not in feature_intervals:
                continue
            
            # Find overlapping features
            overlapping = []
            for start, end, feature_id in feature_intervals[chrom]:
                if start <= pos <= end:
                    overlapping.append(feature_id)
            
            if not overlapping:
                continue
            
            sites_in_features += 1
            
            for feature_id in overlapping:
                feature_info = features[feature_id]
                strand = feature_info['strand']
                
                # Normalize alleles
                if strand == '-':
                    ref = reverse_complement(ref)
                    if alt != '.':
                        alt = reverse_complement(alt)
                
                ref_hom, alt_hom = calc_hom_counts(genotypes)
                nx = ref_hom + alt_hom
                
                if nx == 0:
                    continue
                
                # Compute π
                if alt == 'C':
                    q_C = 1.0 if nx == alt_hom else (0.0 if nx == ref_hom else alt_hom / nx)
                    pi_C = 2.0 * ref_hom * alt_hom / (nx * (nx - 1)) if nx > 1 else 0.0
                else:
                    q_C, pi_C = 0.0, 0.0
                
                if alt == 'G':
                    q_G = 1.0 if nx == alt_hom else (0.0 if nx == ref_hom else alt_hom / nx)
                    pi_G = 2.0 * ref_hom * alt_hom / (nx * (nx - 1)) if nx > 1 else 0.0
                else:
                    q_G, pi_G = 0.0, 0.0
                
                feature_stats[feature_id]['sites'] += 1
                feature_stats[feature_id]['q_C_sum'] += q_C
                feature_stats[feature_id]['pi_C_sum'] += pi_C
                feature_stats[feature_id]['q_G_sum'] += q_G
                feature_stats[feature_id]['pi_G_sum'] += pi_G
    
    finally:
        if vcf_file != '/dev/stdin' and vcf_file != '-':
            f.close()
    
    sys.stderr.write(f"Processed {sites_processed} sites; {sites_in_features} fell in features\n")
    return feature_stats


def write_output(features, feature_stats, output_path, buffer_bytes=16 * 1024 * 1024):
    """Write per-feature π to CSV."""
    try:
        f = open(output_path, 'w', buffering=buffer_bytes, newline='')
    except Exception as e:
        sys.stderr.write(f"Error opening output: {e}\n")
        sys.exit(1)
    
    writer = csv.writer(f)
    writer.writerow(['feature_id', 'feature_type', 'chrom', 'start', 'end', 'length_bp', 'gene_id', 'n_sites', 'q_C', 'pi_C', 'q_G', 'pi_G'])
    
    for feature_id in sorted(feature_stats.keys()):
        stats = feature_stats[feature_id]
        n_sites = stats['sites']
        
        if n_sites == 0:
            continue
        
        feature_info = features[feature_id]
        q_C = stats['q_C_sum'] / n_sites
        pi_C = stats['pi_C_sum'] / n_sites
        q_G = stats['q_G_sum'] / n_sites
        pi_G = stats['pi_G_sum'] / n_sites
        
        feature_len = feature_info['end'] - feature_info['start'] + 1
        
        writer.writerow([
            feature_id,
            feature_info['type'],
            feature_info['chrom'],
            feature_info['start'],
            feature_info['end'],
            feature_len,
            feature_info['gene_id'],
            n_sites,
            f"{q_C:.8f}",
            f"{pi_C:.8f}",
            f"{q_G:.8f}",
            f"{pi_G:.8f}"
        ])
    
    f.close()
    sys.stderr.write(f"Wrote output to {output_path}\n")


def main():
    parser = argparse.ArgumentParser(description='Calculate per-feature π for exons and introns.')
    parser.add_argument('vcf', help='VCF file (or /dev/stdin)')
    parser.add_argument('gff3', help='GFF3 file with CDS features')
    parser.add_argument('output', help='Output CSV')
    parser.add_argument('--buffer-mb', type=int, default=16)
    
    args = parser.parse_args()
    buffer_bytes = args.buffer_mb * 1024 * 1024
    
    sys.stderr.write(f"Parsing GFF3: {args.gff3}\n")
    features, gene_to_strand = parse_gff3_features(args.gff3, buffer_bytes=buffer_bytes)
    
    sys.stderr.write(f"Processing VCF: {args.vcf}\n")
    feature_stats = process_vcf(args.vcf, features, buffer_bytes=buffer_bytes)
    
    sys.stderr.write(f"Writing output: {args.output}\n")
    write_output(features, feature_stats, args.output, buffer_bytes=buffer_bytes)


if __name__ == '__main__':
    main()
