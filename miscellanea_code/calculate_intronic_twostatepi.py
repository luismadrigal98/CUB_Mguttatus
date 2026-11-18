#!/usr/bin/env python3
"""
Two-allele statistics and a per-gene site frequency spectrum at intronic sites,
for the C system (C vs non-C) and the G system (G vs non-G).

The introns of each gene are the no-selection reference for its 4-fold sites:
same chromatin, same transcription, same linked selection. Three analyses in
main_CUB.R read this script's output:
  * the exon/intron diversity ratio against Wright's prediction (Referee 1);
  * the preference-polarised intronic SFS by expression region, the null for
    the 4-fold asymmetry test (Referee 2);
  * the Branch A pi validation.

INPUT MUST BE THE ALL-SITES VCF
-------------------------------
Every intronic site has to be present, monomorphic ones included, because the
per-site denominators and the composition both come from them. The
`*_ready*.vcf` files have had their invariant sites removed by vcf_cleaner.py
(277 M of 315 M lines). Run on one of those, this script used to report a C
frequency of ~0.32 (the composition of variable sites, not of introns) over
~19% of intronic positions. It now refuses: if fewer than half the intronic
sites it reads are invariant it exits with status 2 and writes nothing.
Use Lines187.to.IM767v2.1.vcf.gz.

Intron definition
-----------------
Primary transcript only (`longest=1`), with every exonic feature it has (CDS,
UTRs, exon) merged before taking the gaps, so UTR introns count and a gap is
never part of another isoform's exon. Each intron is trimmed by --trim-bp at
both ends to drop splice signals and kept only if the trimmed width is at
least --min-width. The same definition is implemented in R by
intronic_reference_composition() (src/intronic_controls.R), so the population
data and the reference-sequence composition describe the same bases.

Site filter
-----------
Biallelic SNPs and invariant sites only (indels and multi-allelic records are
skipped). A line is counted when it is a strict homozygote with read support
(0/0 with ref depth > 5x alt depth, or 1/1 the other way round); the site is
kept if at least --min-called lines qualify. Heterozygous calls are residual
heterozygosity in inbred lines and are ignored.

Outputs
-------
<output> (CSV), one row per gene:
    Gene          ID without the 'MgIM767.' prefix and '.v2.1' suffix
    n_sites       intronic sites passing the filter
    q_pref_C      mean C frequency over all n_sites (composition)
    pi_2allele_C  C-vs-non-C heterozygosity per site, 2k(n-k)/(n(n-1)) summed / n_sites
    q_pref_G, pi_2allele_G   the same for G
    n_seg_C, n_fixed_C       segregating sites, and sites fixed for C, in the C system
    n_seg_G, n_fixed_G       the same for G
    pi_total      total heterozygosity per site, 2 r a / (n (n-1)) summed / n_sites
--mask NAME=BED (repeatable): the same table again as <output stem>_<NAME>.csv,
    restricted to the BED's positions. Used with the intronic sites aligned to an
    outgroup (miscellanea_code/divergence_by_gene.py), so that intronic pi and
    divergence are measured on one set of sites. All masks share the one VCF pass.
--outgroup NAME=FASTA (repeatable): like --mask, with the mask taken from the
    outgroup base track that divergence_by_gene.py writes
    (intron_outgroup_<NAME>.fa.gz), and five more columns on <output stem>_<NAME>.csv:
    n_out        sites with an outgroup base (= n_sites)
    out_C        sites where the outgroup carries C (transcript orientation)
    sub_ref      reference base != outgroup base
    fixed_sub    outgroup base absent from every called line (a fixed difference)
    fixed_subC   C-vs-other fixed difference: outgroup C and no line C, or
                 outgroup non-C and every line C
    Counting against the reference picks up polymorphism (at 4-fold sites 37% of
    JUNG1 and 72% of KONZA reference differences fall on segregating sites), so
    the R side uses the fixed counts.
--sfs-out (CSV), segregating sites only:
    Gene, System (C or G), n, k (count of the focal base), count
    Monomorphic sites are the per-gene columns above: fixed for the focal base
    = n_fixed_X, absent = n_sites - n_seg_X - n_fixed_X.

Author: Luis Javier Madrigal-Roca & John K. Kelly
"""

import sys
import argparse
import bisect
import re
import multiprocessing
from collections import defaultdict

INVARIANT_ALTS = {'.', '<NON_REF>', '<*>', '*'}
BASES = {'A', 'C', 'G', 'T'}
COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
EXONIC_FEATURES = {'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR'}
MIN_HOM_RATIO = 5

# Per-gene accumulator layout (a plain list keeps pickling cheap)
N_SITES, C_FREQ, C_PI, G_FREQ, G_PI, C_SEG, C_FIX, G_SEG, G_FIX, TOT_PI = range(10)
N_FIELDS = 10
# Outgroup tables carry five more fields
D_N, D_OUTC, D_SUBREF, D_FIX, D_FIXC = range(10, 15)
N_FIELDS_OUT = 15

TALLY_KEYS = ('intronic_lines', 'used', 'used_invariant', 'skip_nonsnp',
              'skip_format', 'skip_lowcov')


def short_gene_id(gene_id):
    s = re.sub(r'^MgIM767\.', '', gene_id)
    return re.sub(r'\.v\d+\.\d+$', '', s)


# ---------------------------------------------------------------------------
# GFF3 -> trimmed intron intervals
# ---------------------------------------------------------------------------

def parse_gff3_introns(gff3_file, trim_bp=30, min_width=86, primary_only=True,
                       buffer_bytes=16 * 1024 * 1024):
    """
    Returns (intervals, gene_strand, summary).

    intervals: chrom -> (starts, ends, max_end, genes), sorted by start, where
    max_end[i] = max(ends[:i+1]) makes the containment lookup exact and fast
    even when introns of different genes overlap.
    """
    tx_gene, tx_chrom, tx_strand, tx_primary = {}, {}, {}, {}
    tx_parts = defaultdict(list)

    with open(gff3_file, buffering=buffer_bytes) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            chrom, feat, strand = cols[0], cols[2], cols[6]
            if feat not in EXONIC_FEATURES and feat not in ('mRNA', 'transcript'):
                continue
            attr = {}
            for item in cols[8].split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    attr[k.strip()] = v.strip()
            start, end = int(cols[3]), int(cols[4])
            if feat in ('mRNA', 'transcript'):
                tid = attr.get('ID', '')
                if tid:
                    tx_gene[tid] = attr.get('Parent', '')
                    tx_chrom[tid] = chrom
                    tx_strand[tid] = strand
                    tx_primary[tid] = attr.get('longest', '') == '1'
            else:
                for tid in attr.get('Parent', '').split(','):
                    tid = tid.strip()
                    if tid:
                        tx_parts[tid].append((start, end))

    use_primary = primary_only and any(tx_primary.values())
    required = 2 * trim_bp + min_width
    raw = defaultdict(list)
    gene_strand = {}
    n_tx = 0

    for tid, parts in tx_parts.items():
        if tid not in tx_gene or (use_primary and not tx_primary.get(tid, False)):
            continue
        n_tx += 1
        gid = short_gene_id(tx_gene[tid])
        gene_strand[gid] = tx_strand[tid]
        # CDS and UTR pieces of one exon abut; merge them into exons first.
        merged = []
        for s, e in sorted(parts):
            if merged and s <= merged[-1][1] + 1:
                if e > merged[-1][1]:
                    merged[-1][1] = e
            else:
                merged.append([s, e])
        for (_, e1), (s2, _) in zip(merged, merged[1:]):
            o_start, o_end = e1 + 1, s2 - 1          # GFF3 1-based inclusive
            if o_end - o_start + 1 < required:
                continue
            raw[tx_chrom[tid]].append((o_start + trim_bp, o_end - trim_bp, gid))

    intervals, n_introns, union_bp = {}, 0, 0
    for chrom, ivs in raw.items():
        ivs.sort()
        starts = [iv[0] for iv in ivs]
        ends = [iv[1] for iv in ivs]
        genes = [iv[2] for iv in ivs]
        max_end, m = [], -1
        for e in ends:
            m = e if e > m else m
            max_end.append(m)
        intervals[chrom] = (starts, ends, max_end, genes)
        n_introns += len(ivs)
        cur_s, cur_e = None, None                    # union length, for coverage
        for s, e, _ in ivs:
            if cur_e is None or s > cur_e:
                if cur_e is not None:
                    union_bp += cur_e - cur_s + 1
                cur_s, cur_e = s, e
            elif e > cur_e:
                cur_e = e
        if cur_e is not None:
            union_bp += cur_e - cur_s + 1

    summary = {'transcripts': n_tx, 'primary_only': use_primary,
               'introns': n_introns, 'genes': len(set(g for v in intervals.values() for g in v[3])),
               'union_bp': union_bp}
    return intervals, gene_strand, summary


def find_gene(chrom_data, pos):
    """Gene of the first trimmed intron containing pos, or None."""
    starts, ends, max_end, genes = chrom_data
    idx = bisect.bisect_right(starts, pos) - 1
    while idx >= 0 and max_end[idx] >= pos:
        if ends[idx] >= pos:
            return genes[idx]
        idx -= 1
    return None


# ---------------------------------------------------------------------------
# Per-line processing
# ---------------------------------------------------------------------------

def strict_homozygote_counts(fmt, samples):
    """(ref_hom, alt_hom), or None when GT or AD is absent from FORMAT."""
    fields = fmt.split(':')
    try:
        gi, ai = fields.index('GT'), fields.index('AD')
    except ValueError:
        return None
    need = max(gi, ai)
    ref_hom = alt_hom = 0
    for s in samples:
        p = s.split(':')
        if len(p) <= need:
            continue
        ad = p[ai]
        if ad in ('.', ''):
            continue
        try:
            if ',' in ad:
                d = ad.split(',')
                rc = int(d[0]) if d[0] != '.' else 0
                ac = int(d[1]) if d[1] != '.' else 0
            else:
                rc, ac = int(ad), 0
        except ValueError:
            continue
        if rc == 0 and ac == 0:
            continue
        gt = p[gi]
        if gt in ('0/0', '0|0') and rc > MIN_HOM_RATIO * ac:
            ref_hom += 1
        elif gt in ('1/1', '1|1') and ac > MIN_HOM_RATIO * rc:
            alt_hom += 1
    return ref_hom, alt_hom


def add_site(rec, k_c, k_g, n, het_total):
    """One site's contribution to a per-gene record (both two-allele systems
    plus total heterozygosity)."""
    rec[N_SITES] += 1
    for k, i_freq, i_pi, i_seg, i_fix in ((k_c, C_FREQ, C_PI, C_SEG, C_FIX),
                                         (k_g, G_FREQ, G_PI, G_SEG, G_FIX)):
        rec[i_freq] += k / n
        if 0 < k < n:
            rec[i_pi] += 2.0 * k * (n - k) / (n * (n - 1.0))
            rec[i_seg] += 1
        elif k == n:
            rec[i_fix] += 1
    rec[TOT_PI] += het_total


def load_bed_mask(path):
    """BED (0-based half-open) -> {chrom: (starts, ends)} of merged intervals."""
    by_chrom = defaultdict(list)
    with open(path) as fh:
        for ln in fh:
            f = ln.split()
            if len(f) >= 3 and not ln.startswith(('#', 'track', 'browser')):
                by_chrom[f[0]].append((int(f[1]), int(f[2])))
    out = {}
    for chrom, iv in by_chrom.items():
        iv.sort()
        starts, ends = [], []
        for s0, e0 in iv:
            if starts and s0 <= ends[-1]:
                ends[-1] = max(ends[-1], e0)
            else:
                starts.append(s0); ends.append(e0)
        out[chrom] = (starts, ends)
    return out


def load_outgroup_track(path):
    """divergence_by_gene.py intron_outgroup_<T>.fa.gz -> {chrom: bytes}, one byte
    per reference position (0-based index), the outgroup base or '-'."""
    import gzip
    opener = gzip.open if path.endswith('.gz') else open
    out, name, buf = {}, None, []
    with opener(path, 'rb') as fh:
        for ln in fh:
            if ln[:1] == b'>':
                if name is not None:
                    out[name] = b''.join(buf)
                name, buf = ln[1:].split()[0].decode(), []
            else:
                buf.append(ln.rstrip(b'\r\n'))
    if name is not None:
        out[name] = b''.join(buf)
    return out


def add_divergence(rec, out_base, ref, alt, invariant, ref_hom, alt_hom, k_c, n):
    """One site's outgroup comparison (transcript orientation throughout)."""
    rec[D_N] += 1
    if out_base == 'C':
        rec[D_OUTC] += 1
    if out_base != ref:
        rec[D_SUBREF] += 1
    if out_base == ref:
        out_count = ref_hom
    elif not invariant and out_base == alt:
        out_count = alt_hom
    else:
        out_count = 0
    if out_count == 0:
        rec[D_FIX] += 1
    if (out_base == 'C' and k_c == 0) or (out_base != 'C' and k_c == n):
        rec[D_FIXC] += 1


_TRACK_BASES = {ord(b): b for b in 'ACGT'}


def in_mask(mask, chrom, pos1):
    d = mask.get(chrom)
    if d is None:
        return False
    starts, ends = d
    i = bisect.bisect_right(starts, pos1 - 1) - 1
    return i >= 0 and pos1 - 1 < ends[i]


def process_line(line, intervals, gene_strand, min_called, stats, sfs, tally,
                 target_chrom=None, masks=None, mstats=None, outgroups=None):
    if line.startswith('#'):
        return
    f = line.split('\t', 5)
    if len(f) < 6:
        return
    chrom = f[0]
    if target_chrom and chrom != target_chrom:
        return
    cd = intervals.get(chrom)
    if cd is None:
        return
    try:
        pos = int(f[1])
    except ValueError:
        return
    gene_id = find_gene(cd, pos)
    if gene_id is None:
        return
    tally['intronic_lines'] += 1

    ref, alt = f[3].upper(), f[4].upper()
    invariant = alt in INVARIANT_ALTS
    if ref not in BASES or (not invariant and alt not in BASES):
        tally['skip_nonsnp'] += 1                    # indel, multi-allelic, N
        return

    rest = f[5].rstrip('\n').split('\t')             # QUAL FILTER INFO FORMAT samples...
    if len(rest) < 5:
        tally['skip_format'] += 1
        return
    counts = strict_homozygote_counts(rest[3], rest[4:])
    if counts is None:
        tally['skip_format'] += 1
        return
    ref_hom, alt_hom = counts
    if invariant:
        alt_hom = 0
    n = ref_hom + alt_hom
    if n < min_called:
        tally['skip_lowcov'] += 1
        return

    if gene_strand.get(gene_id, '+') == '-':         # transcript orientation
        ref = COMPLEMENT[ref]
        if not invariant:
            alt = COMPLEMENT[alt]

    tally['used'] += 1
    if invariant:
        tally['used_invariant'] += 1

    def focal_count(base):
        if ref == base:
            return ref_hom
        if not invariant and alt == base:
            return alt_hom
        return 0

    k_c, k_g = focal_count('C'), focal_count('G')
    het_total = 2.0 * ref_hom * alt_hom / (n * (n - 1.0))
    rec = stats.get(gene_id)
    if rec is None:
        rec = stats[gene_id] = [0] * N_FIELDS
    add_site(rec, k_c, k_g, n, het_total)
    if 0 < k_c < n:
        sfs[(gene_id, 'C', n, k_c)] += 1
    if 0 < k_g < n:
        sfs[(gene_id, 'G', n, k_g)] += 1

    # Same site, counted again for every mask that contains it (e.g. the
    # intronic positions aligned to an outgroup, so pi and d share their sites).
    if masks:
        for name, mask in masks.items():
            if in_mask(mask, chrom, pos):
                mrec = mstats[name].get(gene_id)
                if mrec is None:
                    mrec = mstats[name][gene_id] = [0] * N_FIELDS
                add_site(mrec, k_c, k_g, n, het_total)
    # Outgroup tracks: the site is in the table if the outgroup has a base here.
    if outgroups:
        minus = gene_strand.get(gene_id, '+') == '-'
        for name, track in outgroups.items():
            t = track.get(chrom)
            if t is None or pos > len(t):
                continue
            b = _TRACK_BASES.get(t[pos - 1])
            if b is None:
                continue
            if minus:
                b = COMPLEMENT[b]
            mrec = mstats[name].get(gene_id)
            if mrec is None:
                mrec = mstats[name][gene_id] = [0] * N_FIELDS_OUT
            add_site(mrec, k_c, k_g, n, het_total)
            add_divergence(mrec, b, ref, alt, invariant, ref_hom, alt_hom, k_c, n)


# ---------------------------------------------------------------------------
# Parallel driver
# ---------------------------------------------------------------------------

WORKER = {}


def worker_init(intervals, gene_strand, min_called, target_chrom, masks=None, outgroups=None):
    WORKER.update(intervals=intervals, gene_strand=gene_strand,
                  min_called=min_called, target_chrom=target_chrom, masks=masks or {},
                  outgroups=outgroups or {})


def process_batch(lines):
    stats, sfs = {}, defaultdict(int)
    tally = dict.fromkeys(TALLY_KEYS, 0)
    mstats = {name: {} for name in list(WORKER['masks']) + list(WORKER['outgroups'])}
    for line in lines:
        process_line(line, WORKER['intervals'], WORKER['gene_strand'],
                     WORKER['min_called'], stats, sfs, tally, WORKER['target_chrom'],
                     WORKER['masks'], mstats, WORKER['outgroups'])
    return stats, dict(sfs), tally, mstats


def _merge_stats(stats, p_stats):
    for g, rec in p_stats.items():
        a = stats.get(g)
        if a is None:
            stats[g] = rec
        else:
            for i in range(len(rec)):
                a[i] += rec[i]


def merge_into(acc, part):
    stats, sfs, tally, mstats = acc
    p_stats, p_sfs, p_tally, p_mstats = part
    _merge_stats(stats, p_stats)
    for key, c in p_sfs.items():
        sfs[key] += c
    for k, v in p_tally.items():
        tally[k] += v
    for name, ms in p_mstats.items():
        _merge_stats(mstats.setdefault(name, {}), ms)


def process_vcf(vcf_file, intervals, gene_strand, min_called, target_chrom=None,
                buffer_bytes=16 * 1024 * 1024, workers=1, batch_size=20000, masks=None,
                outgroups=None):
    masks = masks or {}
    outgroups = outgroups or {}
    acc = ({}, defaultdict(int), dict.fromkeys(TALLY_KEYS, 0),
           {name: {} for name in list(masks) + list(outgroups)})
    line_count = 0

    if workers <= 1:
        worker_init(intervals, gene_strand, min_called, target_chrom, masks, outgroups)
        batch = []
        with open(vcf_file, buffering=buffer_bytes) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                line_count += 1
                batch.append(line)
                if len(batch) >= batch_size:
                    merge_into(acc, process_batch(batch))
                    batch = []
                if line_count % 10_000_000 == 0:
                    print(f"  {line_count:,} lines read...", file=sys.stderr)
            if batch:
                merge_into(acc, process_batch(batch))
        return acc, line_count

    print(f"  Parallel mode: workers={workers}, batch_size={batch_size}", file=sys.stderr)
    with multiprocessing.Pool(processes=workers, initializer=worker_init,
                              initargs=(intervals, gene_strand, min_called,
                                        target_chrom, masks, outgroups)) as pool:
        pending, batch = [], []
        with open(vcf_file, buffering=buffer_bytes) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                line_count += 1
                batch.append(line)
                if len(batch) >= batch_size:
                    pending.append(pool.apply_async(process_batch, (batch,)))
                    batch = []
                    if len(pending) >= workers * 8:          # bound memory
                        for _ in range(workers * 4):
                            merge_into(acc, pending.pop(0).get())
                if line_count % 10_000_000 == 0:
                    print(f"  {line_count:,} lines read...", file=sys.stderr)
            if batch:
                pending.append(pool.apply_async(process_batch, (batch,)))
        for job in pending:
            merge_into(acc, job.get())
    return acc, line_count


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_gene_table(stats, path, buffer_bytes, outgroup=False):
    with open(path, 'w', buffering=buffer_bytes) as out:
        out.write("Gene,n_sites,q_pref_C,pi_2allele_C,q_pref_G,pi_2allele_G,"
                  "n_seg_C,n_fixed_C,n_seg_G,n_fixed_G,pi_total"
                  + (",n_out,out_C,sub_ref,fixed_sub,fixed_subC" if outgroup else "") + "\n")
        for g in sorted(stats):
            r = stats[g]
            n = r[N_SITES]
            if n == 0:
                continue
            out.write(f"{g},{n},{r[C_FREQ] / n:.8f},{r[C_PI] / n:.8f},"
                      f"{r[G_FREQ] / n:.8f},{r[G_PI] / n:.8f},"
                      f"{r[C_SEG]},{r[C_FIX]},{r[G_SEG]},{r[G_FIX]},{r[TOT_PI] / n:.8f}"
                      + (f",{r[D_N]},{r[D_OUTC]},{r[D_SUBREF]},{r[D_FIX]},{r[D_FIXC]}"
                         if outgroup else "") + "\n")


def write_sfs(sfs, path, buffer_bytes):
    with open(path, 'w', buffering=buffer_bytes) as out:
        out.write("Gene,System,n,k,count\n")
        for (g, system, n, k) in sorted(sfs):
            out.write(f"{g},{system},{n},{k},{sfs[(g, system, n, k)]}\n")


def main():
    ap = argparse.ArgumentParser(
        description="Intronic two-allele statistics and per-gene SFS (C and G systems)")
    ap.add_argument('vcf_file', help="ALL-SITES VCF (invariant + variant); /dev/stdin to stream")
    ap.add_argument('gff3_file', help="Genome annotation (GFF3)")
    ap.add_argument('output', help="Per-gene CSV")
    ap.add_argument('--sfs-out', default=None,
                    help="Per-gene segregating-site SFS CSV (default: <output stem>_sfs_by_gene.csv)")
    ap.add_argument('--trim-bp', type=int, default=30)
    ap.add_argument('--min-width', type=int, default=86,
                    help="Minimum TRIMMED intron width; original must be >= 2*trim + this")
    ap.add_argument('--min-called', type=int, default=10,
                    help="Minimum strict-homozygote lines for a site to count (default 10)")
    ap.add_argument('--all-isoforms', action='store_true',
                    help="Use every mRNA, not only the primary (longest=1) transcript")
    ap.add_argument('--chrom', default=None, help="Restrict to one chromosome")
    ap.add_argument('--buffer-mb', type=int, default=16)
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--batch-size', type=int, default=20000)
    ap.add_argument('--mask', action='append', default=[], metavar='NAME=BED',
                    help="Also write <output stem>_<NAME>.csv restricted to the BED's "
                         "positions (0-based half-open), e.g. the intronic sites aligned "
                         "to an outgroup. Repeatable; all masks share one VCF pass.")
    ap.add_argument('--outgroup', action='append', default=[], metavar='NAME=FASTA',
                    help="Outgroup base track from divergence_by_gene.py "
                         "(intron_outgroup_<NAME>.fa.gz): writes <output stem>_<NAME>.csv "
                         "over the sites the outgroup aligns, with fixed-difference "
                         "counts. Repeatable; shares the one VCF pass.")
    ap.add_argument('--allow-variants-only', action='store_true',
                    help="Skip the invariant-site check (results are then NOT per-site)")
    args = ap.parse_args()

    buffer_bytes = max(1, args.buffer_mb) * 1024 * 1024
    sfs_out = args.sfs_out or re.sub(r'\.csv$', '', args.output) + '_sfs_by_gene.csv'

    print(f"Parsing GFF3: {args.gff3_file}", file=sys.stderr)
    intervals, gene_strand, gsum = parse_gff3_introns(
        args.gff3_file, trim_bp=args.trim_bp, min_width=args.min_width,
        primary_only=not args.all_isoforms, buffer_bytes=buffer_bytes)
    print(f"  {gsum['introns']:,} trimmed introns in {gsum['genes']:,} genes from "
          f"{gsum['transcripts']:,} transcripts (primary only: {gsum['primary_only']}); "
          f"{gsum['union_bp']:,} bp", file=sys.stderr)

    masks = {}
    for spec in args.mask:
        name, path = spec.split('=', 1)
        masks[name] = load_bed_mask(path)
        n_bp = sum(e - s for st, en in masks[name].values() for s, e in zip(st, en))
        print(f"  mask {name}: {n_bp:,} bp from {path}", file=sys.stderr)
    outgroups = {}
    for spec in args.outgroup:
        name, path = spec.split('=', 1)
        if name in masks:
            sys.exit(f"ERROR: {name} given as both --mask and --outgroup")
        outgroups[name] = load_outgroup_track(path)
        n_bp = sum(len(t) - t.count(b'-') for t in outgroups[name].values())
        print(f"  outgroup {name}: {n_bp:,} aligned intronic bp from {path}", file=sys.stderr)

    print(f"Processing VCF: {args.vcf_file}  (min_called = {args.min_called})", file=sys.stderr)
    (stats, sfs, tally, mstats), n_lines = process_vcf(
        args.vcf_file, intervals, gene_strand, args.min_called,
        target_chrom=args.chrom, buffer_bytes=buffer_bytes,
        workers=max(1, args.workers), batch_size=max(1000, args.batch_size), masks=masks,
        outgroups=outgroups)

    used = tally['used']
    inv_share = tally['used_invariant'] / used if used else 0.0
    n_tot = sum(r[N_SITES] for r in stats.values()) or 1
    q_c = sum(r[C_FREQ] for r in stats.values()) / n_tot
    q_g = sum(r[G_FREQ] for r in stats.values()) / n_tot
    seg_c = sum(r[C_SEG] for r in stats.values())
    print(f"  VCF data lines read:        {n_lines:,}", file=sys.stderr)
    print(f"  intronic lines:             {tally['intronic_lines']:,}", file=sys.stderr)
    print(f"    skipped indel/multi/N:    {tally['skip_nonsnp']:,}", file=sys.stderr)
    print(f"    skipped FORMAT:           {tally['skip_format']:,}", file=sys.stderr)
    print(f"    skipped < min_called:     {tally['skip_lowcov']:,}", file=sys.stderr)
    print(f"  sites used:                 {used:,} "
          f"({100 * used / max(1, gsum['union_bp']):.1f}% of trimmed intronic bp)", file=sys.stderr)
    print(f"  invariant share of used:    {inv_share:.3f}", file=sys.stderr)
    print(f"  intronic C frequency:       {q_c:.4f}   G: {q_g:.4f}   "
          f"(expect ~0.15 and ~0.17)", file=sys.stderr)
    print(f"  C-system segregating sites: {seg_c:,}", file=sys.stderr)

    if inv_share < 0.5 and not args.allow_variants_only:
        print("ERROR: fewer than half of the intronic sites are invariant, so this is "
              "not an all-sites VCF (the *_ready*.vcf files have had invariant sites "
              "removed). Nothing written. Rerun on Lines187.to.IM767v2.1.vcf.gz.",
              file=sys.stderr)
        sys.exit(2)
    if not 0.10 <= q_c <= 0.22:
        print(f"WARNING: intronic C frequency {q_c:.3f} is outside 0.10-0.22; "
              "check the strand handling and the VCF before using this output.",
              file=sys.stderr)

    print(f"Writing {args.output} ({len(stats):,} genes) and {sfs_out}", file=sys.stderr)
    write_gene_table(stats, args.output, buffer_bytes)
    write_sfs(sfs, sfs_out, buffer_bytes)
    for name, ms in mstats.items():
        mpath = re.sub(r'\.csv$', '', args.output) + f'_{name}.csv'
        n_ms = sum(r[N_SITES] for r in ms.values())
        print(f"Writing {mpath} ({len(ms):,} genes, {n_ms:,} sites inside mask {name})",
              file=sys.stderr)
        write_gene_table(ms, mpath, buffer_bytes, outgroup=name in outgroups)
        if name in outgroups:
            nd = sum(r[D_N] for r in ms.values()) or 1
            print(f"  {name}: divergence {sum(r[D_SUBREF] for r in ms.values()) / nd:.4f} "
                  f"by reference, {sum(r[D_FIX] for r in ms.values()) / nd:.4f} fixed",
                  file=sys.stderr)
    print("Done.", file=sys.stderr)


if __name__ == '__main__':
    main()
