#!/usr/bin/env python3
"""
Calculate two-allele π statistics at intronic sites for C and G nucleotide systems.

For each gene, computes the preferred-allele frequency (q_pref) and per-site
heterozygosity (pi_2allele) treating C vs non-C (C system) and G vs non-G
(G system) as separate two-allele models.  Used to estimate intronic (U, V)
parameters for Wright's mutation-selection-drift model.

C system is relevant for 7/8 four-fold amino-acid families (Ala, Gly, Pro,
Thr, Leu_4, Ser_4, Arg_4 all prefer C-ending codons).  G system applies to
Val (GTG preferred).  Genome-wide intronic composition therefore calibrates
neutral V/U for each system independently.

Input:
    <vcf_file>   VCF with variant AND invariant sites (per chromosome)
    <gff3_file>  Genome annotation (GFF3 format)
    <output>     Output CSV path

Options:
    --trim-bp    Bases trimmed at each splice junction (default: 30)
    --min-width  Minimum trimmed intron width in bp (default: 86)
    --chrom      Restrict processing to this chromosome name

Output CSV columns:
    Gene         gene ID stripped of 'MgIM767.' prefix / '.v2.1' suffix
    n_sites      total intronic sites used
    q_pref_C     mean C-allele frequency across all sites (C system)
    pi_2allele_C per-site C-system heterozygosity (pi_C_sum / n_sites)
    q_pref_G     mean G-allele frequency across all sites (G system)
    pi_2allele_G per-site G-system heterozygosity (pi_G_sum / n_sites)

Author: Luis Javier Madrigal-Roca & John K. Kelly
"""

import sys
import argparse
import bisect
import re
import multiprocessing
from collections import defaultdict


# ---------------------------------------------------------------------------
# GFF3 parsing
# ---------------------------------------------------------------------------

def parse_gff3_introns(gff3_file, trim_bp=30, min_width=86, buffer_bytes=16 * 1024 * 1024):
    """
    Parse GFF3 and derive trimmed intronic intervals per chromosome.

    Introns are the gaps between consecutive exon-like intervals of the same mRNA
    (supports both 'exon' and 'CDS' features).
    Each intron is trimmed by trim_bp on both ends to exclude splice
    signals.  Introns whose original length < 2*trim_bp + min_width are
    discarded (matches get_intron_sequences() defaults in local_M_estimation.R).

    Returns:
        dict of chrom -> (starts, ends, genes) where starts/ends are
        sorted integer lists and genes is the parallel gene-ID list.
        Gene IDs are stripped of 'MgIM767.' prefix and '.v2.1' suffix.
    """
    mrna_to_gene   = {}
    mrna_to_exons  = defaultdict(list)
    gene_to_chrom  = {}
    gene_to_strand = {}

    with open(gff3_file, buffering=buffer_bytes) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            chrom, _, feat, start_s, end_s, _, _, _, attrs_s = cols
            start = int(start_s)
            end   = int(end_s)

            attr = {}
            for item in attrs_s.split(';'):
                item = item.strip()
                if '=' in item:
                    k, v = item.split('=', 1)
                    attr[k.strip()] = v.strip()

            if feat in ('mRNA', 'transcript'):
                mrna_id = attr.get('ID', '')
                parent  = attr.get('Parent', '')
                if mrna_id:
                    mrna_to_gene[mrna_id] = parent
                    gene_to_chrom[parent] = chrom
                    short = re.sub(r'^MgIM767\.', '', parent)
                    short = re.sub(r'\.v\d+\.\d+$', '', short)
                    gene_to_strand[short] = cols[6]

            elif feat in ('exon', 'CDS'):
                parent_raw = attr.get('Parent', '')
                for mrna_id in parent_raw.split(','):
                    mrna_id = mrna_id.strip()
                    if mrna_id:
                        mrna_to_exons[mrna_id].append((start, end))

    # Build intron intervals
    raw = defaultdict(list)   # chrom -> list of (t_start, t_end, gene_id)
    required_orig = 2 * trim_bp + min_width  # minimum original intron length

    for mrna_id, exons in mrna_to_exons.items():
        gene_id = mrna_to_gene.get(mrna_id, '')
        if not gene_id:
            continue
        chrom = gene_to_chrom.get(gene_id, '')
        if not chrom:
            continue

        short = re.sub(r'^MgIM767\.', '', gene_id)
        short = re.sub(r'\.v\d+\.\d+$', '', short)

        exons_sorted = sorted(exons)
        for i in range(len(exons_sorted) - 1):
            e1_end   = exons_sorted[i][1]
            e2_start = exons_sorted[i + 1][0]
            # Intron occupies [e1_end+1, e2_start-1] (GFF3 1-based inclusive)
            orig_start = e1_end + 1
            orig_end   = e2_start - 1
            orig_width = orig_end - orig_start + 1
            if orig_width < required_orig:
                continue
            t_start = orig_start + trim_bp
            t_end   = orig_end   - trim_bp
            raw[chrom].append((t_start, t_end, short))

    # Sort and split into parallel arrays for fast bisect lookup
    result = {}
    for chrom, ivs in raw.items():
        ivs.sort()
        starts = [iv[0] for iv in ivs]
        ends   = [iv[1] for iv in ivs]
        genes  = [iv[2] for iv in ivs]
        result[chrom] = (starts, ends, genes)

    return result, gene_to_strand


def complement_base(base):
    return {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}.get(base, base)


def find_gene(chrom_data, pos):
    """
    Return the gene_id of the first interval containing pos, or None.

    Performs O(log n) bisect followed by a short backward scan to handle
    intervals that end at or after pos.  Multiple overlapping intervals
    (different transcripts of the same gene, or overlapping loci) resolve
    to the first match found.
    """
    if chrom_data is None:
        return None
    starts, ends, genes = chrom_data

    # Find rightmost interval whose start <= pos
    idx = bisect.bisect_right(starts, pos) - 1

    # Scan backward; typical introns are short so this loop is very tight
    while idx >= 0:
        if ends[idx] >= pos:
            return genes[idx]
        # Once the start is more than 1 Mb behind pos, no intron can reach pos
        if pos - starts[idx] > 1_000_000:
            break
        idx -= 1

    return None


# ---------------------------------------------------------------------------
# VCF parsing (mirrors calculate_pi.py exactly)
# ---------------------------------------------------------------------------

def parse_vcf_line(line):
    cols = line.rstrip('\n').split('\t')
    if len(cols) < 10:
        return None
    chrom = cols[0]
    pos   = int(cols[1])
    ref   = cols[3].upper()
    alt   = cols[4].upper()
    genotypes = []
    for j in range(9, len(cols)):
        parts = cols[j].split(':')
        if len(parts) < 3:
            genotypes.append(('./.', 0, 0))
            continue
        gt = parts[0]
        try:
            ad = parts[2]
            if ',' in ad:
                rc, ac = int(ad.split(',')[0]), int(ad.split(',')[1])
            else:
                rc, ac = int(ad), 0
        except (ValueError, IndexError):
            rc = ac = 0
        if rc == 0 and ac == 0:
            genotypes.append(('./.', 0, 0))
            continue
        genotypes.append((gt, rc, ac))
    return chrom, pos, ref, alt, genotypes


def calc_hom_counts(genotypes, min_ratio=5):
    """
    Count strict homozygotes (mirrors calculate_pi.py criterion).

    Returns (ref_hom, alt_hom).
    """
    ref_hom = alt_hom = 0
    for gt, rc, ac in genotypes:
        if gt == '0/0' and rc > min_ratio * ac:
            ref_hom += 1
        elif gt == '1/1' and ac > min_ratio * rc:
            alt_hom += 1
    return ref_hom, alt_hom


# ---------------------------------------------------------------------------
# Per-gene accumulator
# ---------------------------------------------------------------------------

class GeneStats:
    __slots__ = ('n_sites', 'C_freq_sum', 'C_pi_sum', 'G_freq_sum', 'G_pi_sum')

    def __init__(self):
        self.n_sites    = 0
        self.C_freq_sum = 0.0
        self.C_pi_sum   = 0.0
        self.G_freq_sum = 0.0
        self.G_pi_sum   = 0.0


# Worker globals for multiprocessing mode
WORKER_INTERVALS = None
WORKER_GENE_TO_STRAND = None
WORKER_TARGET_CHROM = None


def worker_init(intervals, gene_to_strand, target_chrom):
    global WORKER_INTERVALS, WORKER_GENE_TO_STRAND, WORKER_TARGET_CHROM
    WORKER_INTERVALS = intervals
    WORKER_GENE_TO_STRAND = gene_to_strand
    WORKER_TARGET_CHROM = target_chrom


def add_site_to_stats(stats, gene_id, ref, alt, genotypes):
    if gene_id not in stats:
        stats[gene_id] = [0, 0.0, 0.0, 0.0, 0.0]
    rec = stats[gene_id]
    rec[0] += 1

    # ---- invariant site ----------------------------------------
    if alt == '.':
        if ref == 'C':
            rec[1] += 1.0
        elif ref == 'G':
            rec[3] += 1.0
        # A/T sites contribute 0.0 to both sums (implicit)
        return

    # ---- variant site ------------------------------------------
    ref_hom, alt_hom = calc_hom_counts(genotypes)
    nx = ref_hom + alt_hom
    is_poly = nx >= 2 and min(ref_hom, alt_hom) > 0

    if is_poly:
        # shared pi formula: 2*nx*p*(1-p)/(nx-1)
        p_ref = ref_hom / nx
        pi_val = 2.0 * nx * p_ref * (1.0 - p_ref) / (nx - 1.0)
    else:
        pi_val = 0.0

    # C system: C is preferred (vs non-C)
    if ref == 'C':
        if is_poly:
            q_C = p_ref            # C allele frequency = ref_hom/nx
            rec[2] += pi_val
        else:
            q_C = 1.0 if ref_hom > 0 else 0.0
        rec[1] += q_C
    elif alt == 'C':
        if is_poly:
            q_C = 1.0 - p_ref      # C allele frequency = alt_hom/nx
            rec[2] += pi_val
        else:
            q_C = 1.0              # fixed for alt C allele
        rec[1] += q_C

    # G system: G is preferred (vs non-G)
    if ref == 'G':
        if is_poly:
            q_G = p_ref
            rec[4] += pi_val
        else:
            q_G = 1.0 if ref_hom > 0 else 0.0
        rec[3] += q_G
    elif alt == 'G':
        if is_poly:
            q_G = 1.0 - p_ref
            rec[4] += pi_val
        else:
            q_G = 1.0              # fixed for alt G allele
        rec[3] += q_G


def merge_partial_stats(stats, partial_stats):
    for gene_id, rec in partial_stats.items():
        if gene_id not in stats:
            gs = GeneStats()
            stats[gene_id] = gs
        else:
            gs = stats[gene_id]
        gs.n_sites    += rec[0]
        gs.C_freq_sum += rec[1]
        gs.C_pi_sum   += rec[2]
        gs.G_freq_sum += rec[3]
        gs.G_pi_sum   += rec[4]


def process_batch(lines):
    local_stats = {}
    for line in lines:
        if line.startswith('#'):
            continue

        parsed = parse_vcf_line(line)
        if parsed is None:
            continue
        chrom, pos, ref, alt, genotypes = parsed

        if WORKER_TARGET_CHROM and chrom != WORKER_TARGET_CHROM:
            continue

        chrom_data = WORKER_INTERVALS.get(chrom)
        if chrom_data is None:
            continue

        gene_id = find_gene(chrom_data, pos)
        if gene_id is None:
            continue

        if WORKER_GENE_TO_STRAND.get(gene_id, '+') == '-':
            ref = complement_base(ref)
            if alt != '.':
                alt = ','.join(complement_base(a) for a in alt.split(','))

        add_site_to_stats(local_stats, gene_id, ref, alt, genotypes)

    return local_stats


# ---------------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------------

def process_vcf(vcf_file, intervals, gene_to_strand, target_chrom=None,
                buffer_bytes=16 * 1024 * 1024, workers=1, batch_size=20000):
    """
    Stream VCF, accumulate two-allele stats at intronic sites.

    For each site that falls in a trimmed intron:
      - Invariant (alt == '.'): C/G frequency = 1.0 if REF == 'C'/'G', else 0.0
      - Variant: compute ref_hom / alt_hom; if both > 0 it is polymorphic.
        The C allele is REF when REF=='C', else ALT when ALT=='C', else absent.
        Likewise for G.  A C/G polymorphic site contributes to BOTH systems.

    Returns dict of gene_id -> GeneStats.
    """
    stats = {}
    line_count = 0

    if workers <= 1:
        with open(vcf_file, buffering=buffer_bytes) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                line_count += 1
                if line_count % 1_000_000 == 0:
                    print(f"  {line_count:,} sites processed...", file=sys.stderr)

                parsed = parse_vcf_line(line)
                if parsed is None:
                    continue
                chrom, pos, ref, alt, genotypes = parsed

                if target_chrom and chrom != target_chrom:
                    continue

                chrom_data = intervals.get(chrom)
                if chrom_data is None:
                    continue

                gene_id = find_gene(chrom_data, pos)
                if gene_id is None:
                    continue

                if gene_to_strand.get(gene_id, '+') == '-':
                    ref = complement_base(ref)
                    if alt != '.':
                        alt = ','.join(complement_base(a) for a in alt.split(','))

                partial = {}
                add_site_to_stats(partial, gene_id, ref, alt, genotypes)
                merge_partial_stats(stats, partial)

        print(f"  Total VCF sites read: {line_count:,}", file=sys.stderr)
        return stats

    print(f"  Parallel mode enabled: workers={workers}, batch_size={batch_size}", file=sys.stderr)

    with multiprocessing.Pool(processes=workers,
                              initializer=worker_init,
                              initargs=(intervals, gene_to_strand, target_chrom)) as pool:
        pending = []
        batch = []
        with open(vcf_file, buffering=buffer_bytes) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                line_count += 1
                if line_count % 1_000_000 == 0:
                    print(f"  {line_count:,} sites processed...", file=sys.stderr)

                batch.append(line)
                if len(batch) >= batch_size:
                    pending.append(pool.apply_async(process_batch, (batch,)))
                    batch = []

                    # Apply light backpressure to bound memory during long runs
                    if len(pending) >= workers * 8:
                        drain_n = workers * 4
                        for _ in range(drain_n):
                            merge_partial_stats(stats, pending.pop(0).get())

            if batch:
                pending.append(pool.apply_async(process_batch, (batch,)))

        for job in pending:
            merge_partial_stats(stats, job.get())

    print(f"  Total VCF sites read: {line_count:,}", file=sys.stderr)
    return stats


def write_output(stats, output_path, buffer_bytes=16 * 1024 * 1024):
    with open(output_path, 'w', buffering=buffer_bytes) as out:
        out.write("Gene,n_sites,q_pref_C,pi_2allele_C,q_pref_G,pi_2allele_G\n")
        for gene_id in sorted(stats):
            gs = stats[gene_id]
            n = gs.n_sites
            if n == 0:
                continue
            q_C  = gs.C_freq_sum / n
            pi_C = gs.C_pi_sum   / n
            q_G  = gs.G_freq_sum / n
            pi_G = gs.G_pi_sum   / n
            out.write(f"{gene_id},{n},{q_C:.8f},{pi_C:.8f},{q_G:.8f},{pi_G:.8f}\n")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Two-allele intronic π for C and G systems (Wright UV calibration)"
    )
    parser.add_argument('vcf_file',    help="VCF with invariant+variant sites")
    parser.add_argument('gff3_file',   help="Genome annotation (GFF3)")
    parser.add_argument('output',      help="Output CSV path")
    parser.add_argument('--trim-bp',   type=int, default=30,
                        help="Bases trimmed at each splice junction (default: 30)")
    parser.add_argument('--min-width', type=int, default=86,
                        help="Minimum trimmed intron width in bp (default: 86)")
    parser.add_argument('--chrom',     default=None,
                        help="Restrict to this chromosome (default: all)")
    parser.add_argument('--buffer-mb', type=int, default=16,
                        help="Buffered I/O size in MB for GFF/VCF/CSV (default: 16)")
    parser.add_argument('--workers',   type=int, default=1,
                        help="Worker processes for parallel VCF batches (default: 1)")
    parser.add_argument('--batch-size', type=int, default=20000,
                        help="VCF lines per worker batch in parallel mode (default: 20000)")
    args = parser.parse_args()

    buffer_bytes = max(1, args.buffer_mb) * 1024 * 1024

    print(f"Parsing GFF3: {args.gff3_file}", file=sys.stderr)
    print(f"  trim_bp={args.trim_bp}, min_width={args.min_width}", file=sys.stderr)
    print(f"  io_buffer={args.buffer_mb} MB", file=sys.stderr)
    intervals, gene_to_strand = parse_gff3_introns(args.gff3_file,
                                                   trim_bp=args.trim_bp,
                                                   min_width=args.min_width,
                                                   buffer_bytes=buffer_bytes)
    n_ivs = sum(len(v[0]) for v in intervals.values())
    print(f"  {n_ivs:,} trimmed intron intervals across {len(intervals)} chromosomes",
          file=sys.stderr)

    print(f"Processing VCF: {args.vcf_file}", file=sys.stderr)
    if args.chrom:
        print(f"  Restricting to chromosome: {args.chrom}", file=sys.stderr)
    print(f"  workers={args.workers}, batch_size={args.batch_size}", file=sys.stderr)
    stats = process_vcf(args.vcf_file, intervals, gene_to_strand,
                        target_chrom=args.chrom,
                        buffer_bytes=buffer_bytes,
                        workers=max(1, args.workers),
                        batch_size=max(1000, args.batch_size))
    print(f"  {len(stats):,} genes with intronic data", file=sys.stderr)

    print(f"Writing output: {args.output}", file=sys.stderr)
    write_output(stats, args.output, buffer_bytes=buffer_bytes)
    print("Done.", file=sys.stderr)


if __name__ == '__main__':
    main()
