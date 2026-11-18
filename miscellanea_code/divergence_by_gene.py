#!/usr/bin/env python3
"""
Per-gene divergence from IM767 to outgroup assemblies, from AnchorWave alignments,
at 4-fold degenerate sites, at 0-fold sites and at trimmed intronic sites.

Divergence is the polymorphism-free counterpart of the CUB pipeline's diversity
(main_CUB.R, Section 12.1i): pi / d at the same sites is invariant to the
mutational input (rate or spectrum), rises with S under weak selection, and falls
under linked selection. Every site class here is defined exactly as in the
pipeline, so pi and d can be computed on one set of sites:

  4-fold sites   3rd positions of codons whose first two coding-strand bases are
                 GC, GG, CC, AC, GT, CT, TC or CG (the eight four-fold families),
                 primary transcript (longest=1), reference reading frame.
  intronic sites the trimmed introns of calculate_intronic_twostatepi.py, imported
                 from it: primary transcript, exonic features merged, --trim-bp at
                 each end, trimmed width >= --min-width.

Alignments: AnchorWave proAli MAFs, one per chromosome and taxon, query aligned to
the IM767 reference (`<aln-dir>/<chrom>/<taxon>/aln.maf`). Reference columns that
are gaps in the query, or N, count as unaligned.

Divergence is reference vs outgroup, in transcript orientation:
  sub    reference base != outgroup base (both ACGT)
  subC   (reference is C) != (outgroup is C): the C-vs-other two-allele system of
         Wright's model, used with the two-allele pi.

Outputs (in --out-dir):
  fourfold_sites.tsv.gz          Gene Chrom Pos Strand Ref <T>_base ...
                                 one row per 4-fold site; Pos 1-based (the 3rd
                                 entry of Genomic_Positions in the codon file);
                                 bases on the coding strand, '.' if unaligned.
  intron_divergence_by_gene.tsv  Gene n_sites <T>_n_aln <T>_sub <T>_subC ...
  zerofold_divergence_by_gene.tsv Gene n_sites <T>_n_aln <T>_sub ...
                                 0-fold sites: every other base changes the amino
                                 acid (stop counts as one), the definition of
                                 describe_gene_positions_by_degeneracy.py, so these
                                 match the pipeline's Sites_0fold / Pi_sum_0fold.
  intron_outgroup_<T>.fa.gz      the outgroup's base at every aligned intronic
                                 position, '-' elsewhere: one line per chromosome,
                                 IM767 coordinates, + strand. The --outgroup input of
                                 calculate_intronic_twostatepi.py, which counts FIXED
                                 intronic differences (the outgroup base absent from
                                 every line) in the same VCF pass as pi.
  intron_aligned_<T>.bed         merged aligned intronic positions per taxon,
                                 0-based half-open; the --mask for the intronic
                                 polymorphism script, so pi uses the same sites.
  summary_by_chrom.tsv           aligned fractions and raw divergence, per
                                 chromosome and taxon (the reference-bias check:
                                 KONZA should be less diverged than JUNG1).

Author: Luis Javier Madrigal-Roca & John K. Kelly
"""
import argparse
import gzip
import os
import re
import sys

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from calculate_intronic_twostatepi import parse_gff3_introns, short_gene_id  # noqa: E402

FOURFOLD_PREFIX = {b"GC", b"GG", b"CC", b"AC", b"GT", b"CT", b"TC", b"CG"}
COMP = bytes.maketrans(b"ACGTN", b"TGCAN")
ACGT = np.frombuffer(b"ACGT", dtype=np.uint8)
GAP, DOT = ord("-"), ord(".")


def load_genome(path, wanted):
    seqs, name, buf = {}, None, []
    with open(path, "rb") as fh:
        for ln in fh:
            if ln[:1] == b">":
                if name in wanted:
                    seqs[name] = b"".join(buf).upper()
                name, buf = ln[1:].split()[0].decode(), []
            elif name in wanted:
                buf.append(ln.strip())
    if name in wanted:
        seqs[name] = b"".join(buf).upper()
    missing = set(wanted) - set(seqs)
    if missing:
        raise SystemExit(f"chromosomes missing from genome: {sorted(missing)}")
    return seqs


def primary_cds(gff_path):
    """{chrom: [(gene_short, strand, [(start1, end1, phase), ...]), ...]}"""
    tx, cds = {}, {}
    with open(gff_path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] == "mRNA" and "longest=1" in f[8]:
                tid = re.search(r"ID=([^;]+)", f[8]).group(1)
                par = re.search(r"Parent=([^;]+)", f[8]).group(1)
                tx[tid] = (f[0], short_gene_id(par), f[6])
            elif f[2] == "CDS":
                for par in re.search(r"Parent=([^;]+)", f[8]).group(1).split(","):
                    cds.setdefault(par, []).append((int(f[3]), int(f[4]), f[7]))
    out = {}
    for tid, (chrom, gene, strand) in tx.items():
        if tid in cds:
            out.setdefault(chrom, []).append((gene, strand, cds[tid]))
    return out


_BASES = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
GENETIC_CODE = {(a + b + c).encode(): _AA[16 * i + 4 * j + k]
                for i, a in enumerate(_BASES) for j, b in enumerate(_BASES)
                for k, c in enumerate(_BASES)}


def _zerofold_offsets(codon):
    """Codon positions (0-2) at which every other base changes the amino acid
    (stop counts as an amino acid), as in describe_gene_positions_by_degeneracy.py."""
    aa = GENETIC_CODE[codon]
    out = []
    for k in range(3):
        if all(GENETIC_CODE[codon[:k] + b + codon[k + 1:]] != aa
               for b in (b"A", b"C", b"G", b"T") if b != codon[k:k + 1]):
            out.append(k)
    return tuple(out)


ZEROFOLD_OFFSETS = {c: _zerofold_offsets(c) for c in GENETIC_CODE}


def degenerate_positions(ref, strand, parts):
    """0-based genomic positions of 4-fold 3rd codon positions and of 0-fold
    positions, each in coding order."""
    parts = sorted(parts, key=lambda c: c[0], reverse=(strand == "-"))
    ph = parts[0][2]
    phase = int(ph) if ph in ("0", "1", "2") else 0
    coord = []
    for s1, e1, _ in parts:
        coord.extend(range(s1 - 1, e1) if strand == "+" else range(e1 - 1, s1 - 2, -1))
    coding = bytes(ref[p] for p in coord)
    if strand == "-":
        coding = coding.translate(COMP)
    four, zero = [], []
    for i in range(phase, len(coding) - 2, 3):
        if coding[i:i + 2] in FOURFOLD_PREFIX:
            four.append(coord[i + 2])
        for k in ZEROFOLD_OFFSETS.get(coding[i:i + 3], ()):
            zero.append(coord[i + k])
    return four, zero


def maf_query_bases(maf_path, ref_len):
    """uint8 array over the reference: aligned query base (ACGT) or 0."""
    q_at = np.zeros(ref_len, dtype=np.uint8)
    up = bytes.maketrans(b"acgtn", b"ACGTN")
    r_start = r_seq = q_seq = None

    def flush():
        r = np.frombuffer(r_seq.translate(up), dtype=np.uint8)
        q = np.frombuffer(q_seq.translate(up), dtype=np.uint8)
        is_ref = r != GAP
        pos = r_start + np.cumsum(is_ref) - 1
        p, qb = pos[is_ref], q[is_ref]
        keep = (p < ref_len) & np.isin(qb, ACGT)
        q_at[p[keep]] = qb[keep]

    with open(maf_path, "rb") as fh:
        for line in fh:
            c = line[:1]
            if c == b"a":
                if r_seq is not None and q_seq is not None:
                    flush()
                r_start = r_seq = q_seq = None
            elif c == b"s":
                f = line.split()
                if r_seq is None:
                    r_start, r_seq = int(f[2]), f[6]
                elif q_seq is None:
                    q_seq = f[6]
    if r_seq is not None and q_seq is not None:
        flush()
    return q_at


def merged_bed(chrom, positions):
    """Sorted 0-based positions -> list of (chrom, start, end) half-open runs."""
    if positions.size == 0:
        return []
    brk = np.where(np.diff(positions) != 1)[0]
    starts = np.r_[positions[0], positions[brk + 1]]
    ends = np.r_[positions[brk], positions[-1]] + 1
    return [(chrom, int(s), int(e)) for s, e in zip(starts, ends)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--genome", required=True, help="IM767 genome FASTA (the alignment reference)")
    ap.add_argument("--gff", required=True, help="IM767 v2.1 gene GFF3 (the pipeline's)")
    ap.add_argument("--aln-dir", required=True,
                    help="anchorwave_proali_out/chrom (holds <chrom>/<taxon>/aln.maf)")
    ap.add_argument("--taxa", default="JUNG1,KONZA")
    ap.add_argument("--chroms", default="all", help="comma list, or 'all' (Chr_01..Chr_14)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--trim-bp", type=int, default=30)
    ap.add_argument("--min-width", type=int, default=86)
    a = ap.parse_args()

    taxa = a.taxa.split(",")
    chroms = [f"Chr_{i:02d}" for i in range(1, 15)] if a.chroms == "all" else a.chroms.split(",")
    os.makedirs(a.out_dir, exist_ok=True)

    print("Parsing GFF ...", file=sys.stderr)
    cds_by_chrom = primary_cds(a.gff)
    introns, gene_strand, isum = parse_gff3_introns(a.gff, trim_bp=a.trim_bp, min_width=a.min_width)
    print(f"  {isum['introns']:,} trimmed introns, {isum['genes']:,} genes", file=sys.stderr)
    genome = load_genome(a.genome, set(chroms))

    ff = gzip.open(os.path.join(a.out_dir, "fourfold_sites.tsv.gz"), "wt")
    ff.write("Gene\tChrom\tPos\tStrand\tRef\t" + "\t".join(f"{t}_base" for t in taxa) + "\n")
    intr = {}                                       # gene -> [n_sites, (n_aln, sub, subC) x taxa]
    zf = {}                                         # gene -> [n_sites, (n_aln, sub) x taxa]
    beds = {t: [] for t in taxa}
    tracks = {t: gzip.open(os.path.join(a.out_dir, f"intron_outgroup_{t}.fa.gz"), "wb",
                           compresslevel=6) for t in taxa}
    summary = []

    for chrom in chroms:
        ref = genome[chrom]
        n = len(ref)
        ref_arr = np.frombuffer(ref, dtype=np.uint8)
        q = {}
        for t in taxa:
            path = os.path.join(a.aln_dir, chrom, t, "aln.maf")
            if not os.path.exists(path):
                raise SystemExit(f"missing alignment: {path}")
            q[t] = maf_query_bases(path, n)

        # ---- 4-fold and 0-fold sites -----------------------------------------
        n4 = n0 = 0
        al4 = {t: 0 for t in taxa}
        sub4 = {t: 0 for t in taxa}
        al0 = {t: 0 for t in taxa}
        sub0 = {t: 0 for t in taxa}
        for gene, strand, parts in cds_by_chrom.get(chrom, []):
            four, zero = degenerate_positions(ref, strand, parts)
            z = np.asarray(zero, dtype=np.int64)
            z = z[np.isin(ref_arr[z], ACGT)] if z.size else z
            rec0 = zf.setdefault(gene, [0] + [0] * (2 * len(taxa)))
            rec0[0] += int(z.size)
            n0 += int(z.size)
            for j, t in enumerate(taxa):
                qq = q[t][z]
                al = qq != 0
                s = int((al & (qq != ref_arr[z])).sum())
                rec0[1 + 2 * j] += int(al.sum())
                rec0[2 + 2 * j] += s
                al0[t] += int(al.sum())
                sub0[t] += s
            for p in four:
                rb = ref[p]
                if rb not in b"ACGT":
                    continue
                n4 += 1
                cols = []
                for t in taxa:
                    qb = int(q[t][p])
                    if qb:
                        al4[t] += 1
                        sub4[t] += qb != rb
                    cols.append(qb)
                if strand == "-":
                    rb = COMP[rb]
                    cols = [COMP[c] if c else 0 for c in cols]
                ff.write(f"{gene}\t{chrom}\t{p + 1}\t{strand}\t{chr(rb)}\t"
                         + "\t".join(chr(c) if c else "." for c in cols) + "\n")

        # ---- intronic sites ----------------------------------------------------
        nI = 0
        alI = {t: 0 for t in taxa}
        subI = {t: 0 for t in taxa}
        aligned_pos = {t: [] for t in taxa}
        if chrom in introns:
            starts, ends, _mx, genes = introns[chrom]
            for s1, e1, gene in zip(starts, ends, genes):
                lo, hi = s1 - 1, min(e1, n)          # 1-based inclusive -> 0-based half-open
                if hi <= lo:
                    continue
                r = ref_arr[lo:hi]
                ok = np.isin(r, ACGT)
                minus = gene_strand.get(gene, "+") == "-"
                rC = r == (ord("G") if minus else ord("C"))   # C on the transcript strand
                rec = intr.setdefault(gene, [0] + [0] * (3 * len(taxa)))
                rec[0] += int(ok.sum())
                nI += int(ok.sum())
                for j, t in enumerate(taxa):
                    qq = q[t][lo:hi]
                    al = ok & (qq != 0)
                    qC = qq == (ord("G") if minus else ord("C"))
                    rec[1 + 3 * j] += int(al.sum())
                    rec[2 + 3 * j] += int((al & (qq != r)).sum())
                    rec[3 + 3 * j] += int((al & (qC != rC)).sum())
                    alI[t] += int(al.sum())
                    subI[t] += int((al & (qq != r)).sum())
                    aligned_pos[t].append(np.nonzero(al)[0] + lo)
        for t in taxa:
            pos = np.unique(np.concatenate(aligned_pos[t])) if aligned_pos[t] else np.array([], dtype=np.int64)
            beds[t].extend(merged_bed(chrom, pos))
            trk = np.full(n, ord("-"), dtype=np.uint8)
            trk[pos] = q[t][pos]
            tracks[t].write(f">{chrom}\n".encode() + trk.tobytes() + b"\n")
            summary.append((chrom, t, n4, al4[t], sub4[t], nI, alI[t], subI[t], n0, al0[t], sub0[t]))
            print(f"  {chrom} {t:6s} 4-fold aligned {100 * al4[t] / max(n4, 1):5.1f}% "
                  f"p {100 * sub4[t] / max(al4[t], 1):5.2f}% | 0-fold aligned "
                  f"{100 * al0[t] / max(n0, 1):5.1f}% p {100 * sub0[t] / max(al0[t], 1):5.2f}% | intron aligned "
                  f"{100 * alI[t] / max(nI, 1):5.1f}% p {100 * subI[t] / max(alI[t], 1):5.2f}%",
                  file=sys.stderr)
        del q

    ff.close()
    for t in taxa:
        tracks[t].close()
    with open(os.path.join(a.out_dir, "intron_divergence_by_gene.tsv"), "w") as fh:
        fh.write("Gene\tn_sites\t" + "\t".join(f"{t}_n_aln\t{t}_sub\t{t}_subC" for t in taxa) + "\n")
        for g in sorted(intr):
            fh.write(g + "\t" + "\t".join(map(str, intr[g])) + "\n")
    with open(os.path.join(a.out_dir, "zerofold_divergence_by_gene.tsv"), "w") as fh:
        fh.write("Gene\tn_sites\t" + "\t".join(f"{t}_n_aln\t{t}_sub" for t in taxa) + "\n")
        for g in sorted(zf):
            fh.write(g + "\t" + "\t".join(map(str, zf[g])) + "\n")
    for t in taxa:
        with open(os.path.join(a.out_dir, f"intron_aligned_{t}.bed"), "w") as fh:
            for c, s, e in beds[t]:
                fh.write(f"{c}\t{s}\t{e}\n")
    with open(os.path.join(a.out_dir, "summary_by_chrom.tsv"), "w") as fh:
        fh.write("Chrom\tTaxon\tn_4fold\tn_4fold_aln\tsub_4fold\tn_intron\tn_intron_aln\tsub_intron"
                 "\tn_0fold\tn_0fold_aln\tsub_0fold\n")
        for row in summary:
            fh.write("\t".join(map(str, row)) + "\n")
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
