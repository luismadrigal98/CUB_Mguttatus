#!/usr/bin/env python3
"""
Calculate polymorphism statistics (π, Tajima's D) per gene from VCF data.

Input:
    - VCF file with variant AND invariant sites
    - <chrom>.genic_bases.annotated.txt (from describe_gene_positions_by_degeneracy.py)
    - Number of samples (optional, default: 187)

Output:
    - <chrom>.bygene.pi.txt: Per-gene diversity statistics
    
Format:
    Chr  Gene  
    Sites_0fold  Poly_0fold  Pi_sum_0fold  Pi_mean_0fold  ThetaW_0fold  TajimaD_0fold
    Sites_2fold  Poly_2fold  Pi_sum_2fold  Pi_mean_2fold  ThetaW_2fold  TajimaD_2fold
    Sites_3fold  Poly_3fold  Pi_sum_3fold  Pi_mean_3fold  ThetaW_3fold  TajimaD_3fold
    Sites_4fold  Poly_4fold  Pi_sum_4fold  Pi_mean_4fold  ThetaW_4fold  TajimaD_4fold
    Sites_all    Poly_all    Pi_sum_all    Pi_mean_all    ThetaW_all    TajimaD_all

Author: Luis Javier Madrigal-Roca & John K. Kelly
"""

import sys
import math

def load_annotated_sites(annotation_file, chrom):
    """
    Load position annotations from describe_gene_positions_by_degeneracy.py output.
    
    Returns:
        gsites: dict of position -> [gene_id, base, degeneracy_class, ref_codon, amino_acid]
    """
    gsites = {}
    
    with open(annotation_file, 'r') as f:
        header = f.readline()  # Skip header: Chr Gene Position Base Codon_Position Degeneracy Ref_Codon Amino_Acid [Strand]
        
        for line in f:
            cols = line.strip().split('\t')
            # Format: Chr_06	06G000100	29638	A	1	0-fold	ATG	M	[+/-]
            if len(cols) < 8:
                continue
            
            # Use indexing for compatibility with both old (8-col) and new (9-col) formats
            chr_name = cols[0]
            gene_id = cols[1]
            pos = cols[2]
            base = cols[3]
            # codon_pos = cols[4]  # Not used in this script
            degeneracy = cols[5]
            ref_codon = cols[6]
            amino_acid = cols[7]
            # strand = cols[8] if len(cols) > 8 else '+'  # Not used in this script
            
            if chr_name == chrom:
                gsites[int(pos)] = [gene_id, base, degeneracy, ref_codon, amino_acid]
    
    return gsites

def parse_vcf_line(line):
    """
    Parse a single VCF line EXACTLY as proc2.py does.
    
    Format: GT:PL:AD (e.g., "0/0:0:4" or "0/0:0,36,255:12,0")
    - vv[0] = GT (genotype: 0/0, 0/1, 1/1)
    - vv[2] = AD (allele depths: "ref,alt")
    
    Returns:
        (chrom, pos, ref, alt, genotypes)
        genotypes: list of (gt, ref_count, alt_count) for each sample
    """
    cols = line.strip().split('\t')
    
    if len(cols) < 10:
        return None
    
    chrom = cols[0]
    pos = int(cols[1])
    ref = cols[3]
    alt = cols[4]
    
    # Parse genotypes (starting at column 9)
    # Match proc2.py: for j in range(9,len(cols)):
    genotypes = []
    for j in range(9, len(cols)):
        gt_field = cols[j]
        # vv=cols[j].split(":")
        parts = gt_field.split(':')
        
        if len(parts) < 3:
            # Missing data - skip this sample
            genotypes.append(("./.", 0, 0))
            continue
        
        # vv[0] = GT
        gt = parts[0]
        
        # RA=vv[2].split(",")
        try:
            ad_field = parts[2]
            ad_parts = ad_field.split(',')
            ref_count = int(ad_parts[0])
            alt_count = int(ad_parts[1]) if len(ad_parts) > 1 else 0
        except (ValueError, IndexError):
            ref_count, alt_count = 0, 0
        
        genotypes.append((gt, ref_count, alt_count))
    
    return chrom, pos, ref, alt, genotypes

def calculate_pi_site(genotypes, min_depth_ratio=5):
    """
    Calculate nucleotide diversity (π) for a single polymorphic site.
    
    EXACT implementation from proc2.py:
        - Count lines with clear homozygous ref: GT="0/0" AND ref_depth > 5 * alt_depth
        - Count lines with clear homozygous alt: GT="1/1" AND alt_depth > 5 * ref_depth
        - Calculate π = 2 * nx * px * (1-px) / (nx-1)
          where nx = total homozygotes, px = ref frequency
    
    Returns:
        (is_polymorphic, pi_value, n_samples)
    """
    ref_hom = 0  # RefAlt[0] in proc2.py
    alt_hom = 0  # RefAlt[1] in proc2.py
    
    for gt, ref_count, alt_count in genotypes:
        # Match proc2.py logic EXACTLY:
        # if vv[0]=="0/0" and int(RA[0])>5*int(RA[1]):
        if gt == "0/0" and ref_count > min_depth_ratio * alt_count:
            ref_hom += 1
        # elif vv[0]=="1/1" and int(RA[1])>5*int(RA[0]):
        elif gt == "1/1" and alt_count > min_depth_ratio * ref_count:
            alt_hom += 1
    
    # if min(RefAlt)>0:
    if min(ref_hom, alt_hom) > 0:
        # nx = float(sum(RefAlt))
        nx = float(ref_hom + alt_hom)
        # px = float(RefAlt[0])/nx
        px = float(ref_hom) / nx
        # pix = 2*nx*(px*(1-px))/(nx-1.0)
        pi = 2.0 * nx * px * (1.0 - px) / (nx - 1.0)
        return True, pi, int(nx)
    else:
        return False, 0.0, 0

def process_vcf(vcf_file, gsites, chrom):
    """
    Process VCF file and calculate diversity per gene per degeneracy class.
    
    Returns:
        counts: dict of gene -> {degeneracy_class: [invariant_sites, polymorphic_sites, sum_pi]}
    """
    counts = {}
    
    print(f"Processing VCF for {chrom}...")
    
    line_count = 0
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            line_count += 1
            if line_count % 1000000 == 0:
                print(f"  Processed {line_count:,} sites...", file=sys.stderr)
            
            parsed = parse_vcf_line(line)
            if not parsed:
                continue
            
            vcf_chrom, pos, ref, alt, genotypes = parsed
            
            if vcf_chrom != chrom:
                continue
            
            # Check if this position is in a gene
            if pos not in gsites:
                continue
            
            gene_id, expected_base, degeneracy, ref_codon, amino_acid = gsites[pos]
            
            # Initialize gene if not seen
            if gene_id not in counts:
                counts[gene_id] = {
                    "0-fold": [0, 0, 0.0],      # [invariant, polymorphic, sum_pi]
                    "2-fold": [0, 0, 0.0],
                    "3-fold": [0, 0, 0.0],
                    "4-fold": [0, 0, 0.0]
                }
            
            # Check if site is polymorphic
            if alt == '.':
                # Invariant site
                counts[gene_id][degeneracy][0] += 1
            else:
                # Potentially polymorphic
                is_poly, pi_value, n_samples = calculate_pi_site(genotypes)
                
                if is_poly:
                    counts[gene_id][degeneracy][1] += 1
                    counts[gene_id][degeneracy][2] += pi_value
                else:
                    # Not really polymorphic by our criteria
                    counts[gene_id][degeneracy][0] += 1
    
    print(f"  Total sites processed: {line_count:,}")
    return counts

def calculate_theta_w(n_poly, n_total, n_samples):
    """
    Calculate Watterson's theta (θ_W).
    θ_W = S / a_n
    where S = number of segregating sites, a_n = sum(1/i for i in 1..n-1)
    """
    if n_total == 0 or n_samples < 2:
        return 0.0
    
    # Calculate harmonic number a_n
    a_n = sum(1.0 / i for i in range(1, n_samples))
    
    theta_w_per_locus = float(n_poly) / a_n if a_n > 0 else 0.0
    
    theta_w_per_site = theta_w_per_locus / n_total  # Scale to total sites

    return theta_w_per_site

def calculate_tajimas_d(pi, theta_w, n_poly, n_samples):
    """
    Calculate Tajima's D statistic.
    D = (π - θ_W) / sqrt(Var(π - θ_W))
    """
    if n_samples < 4 or theta_w == 0:
        return float('nan')
    
    n = float(n_samples)
    S = float(n_poly)
    
    # Calculate variance terms (Tajima 1989)
    a1 = sum(1.0 / i for i in range(1, int(n)))
    a2 = sum(1.0 / (i * i) for i in range(1, int(n)))
    
    b1 = (n + 1.0) / (3.0 * (n - 1.0))
    b2 = 2.0 * (n * n + n + 3.0) / (9.0 * n * (n - 1.0))
    
    c1 = b1 - 1.0 / a1
    c2 = b2 - (n + 2.0) / (a1 * n) + a2 / (a1 * a1)
    
    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)
    
    var_d = e1 * S + e2 * S * (S - 1.0)
    
    if var_d <= 0:
        return float('nan')
    
    D = (pi - theta_w) / math.sqrt(var_d)
    
    return D

def write_output(counts, chrom, output_file, n_samples):
    """Write per-gene diversity statistics to output file."""
    
    with open(output_file, 'w') as out:
        # Header - NEW DEGENERACY CATEGORIES + overall metrics + neutrality stats for all classes
        out.write("Chr\tGene\t")
        out.write("Sites_0fold\tPoly_0fold\tPi_sum_0fold\tPi_mean_0fold\tThetaW_0fold\tTajimaD_0fold\t")
        out.write("Sites_2fold\tPoly_2fold\tPi_sum_2fold\tPi_mean_2fold\tThetaW_2fold\tTajimaD_2fold\t")
        out.write("Sites_3fold\tPoly_3fold\tPi_sum_3fold\tPi_mean_3fold\tThetaW_3fold\tTajimaD_3fold\t")
        out.write("Sites_4fold\tPoly_4fold\tPi_sum_4fold\tPi_mean_4fold\tThetaW_4fold\tTajimaD_4fold\t")
        out.write("Sites_all\tPoly_all\tPi_sum_all\tPi_mean_all\tThetaW_all\tTajimaD_all\n")
        
        for gene_id in sorted(counts.keys()):
            gene_counts = counts[gene_id]
            
            # Extract counts for each degeneracy class (NEW CATEGORIES)
            pos_0fold = gene_counts["0-fold"]
            pos_2fold = gene_counts["2-fold"]
            pos_3fold = gene_counts["3-fold"]
            pos_4fold = gene_counts["4-fold"]
            
            # Calculate overall totals
            total_sites = pos_0fold[0] + pos_2fold[0] + pos_3fold[0] + pos_4fold[0]
            total_sites += pos_0fold[1] + pos_2fold[1] + pos_3fold[1] + pos_4fold[1]  # Include polymorphic sites
            total_poly = pos_0fold[1] + pos_2fold[1] + pos_3fold[1] + pos_4fold[1]
            total_pi_sum = pos_0fold[2] + pos_2fold[2] + pos_3fold[2] + pos_4fold[2]
            total_pi_mean = total_pi_sum / total_sites if total_sites > 0 else 0.0
            
            # Calculate overall neutrality stats
            theta_w_all = calculate_theta_w(total_poly, total_sites, n_samples)
            tajima_d_all = calculate_tajimas_d(total_pi_mean, theta_w_all, total_poly, n_samples)
            
            # Write output row
            out.write(f"{chrom}\t{gene_id}\t")
            
            # Per-degeneracy class: Sites, Poly, Pi_sum, Pi_mean, ThetaW, TajimaD
            for deg_class in [pos_0fold, pos_2fold, pos_3fold, pos_4fold]:
                n_sites = deg_class[0] + deg_class[1]
                n_poly = deg_class[1]
                pi_sum = deg_class[2]
                pi_mean = pi_sum / n_sites if n_sites > 0 else 0.0
                
                # Calculate neutrality statistics for this degeneracy class
                theta_w = calculate_theta_w(n_poly, n_sites, n_samples)
                tajima_d = calculate_tajimas_d(pi_mean, theta_w, n_poly, n_samples)
                
                out.write(f"{n_sites}\t{n_poly}\t{pi_sum:.6f}\t{pi_mean:.6f}\t{theta_w:.6f}\t{tajima_d:.4f}\t")
            
            # Overall metrics (like proc2.py total)
            out.write(f"{total_sites}\t{total_poly}\t{total_pi_sum:.6f}\t{total_pi_mean:.6f}\t{theta_w_all:.6f}\t{tajima_d_all:.4f}\n")

def main():
    if len(sys.argv) not in [4, 5]:
        print("Usage: python calculate_pi.py <chromosome> <vcf_file> <annotation_file> [n_samples]")
        print("Example: python calculate_pi.py Chr_01 Included_snp.sites.txt Chr_01.genic_bases.annotated.txt 187")
        print("  n_samples: Number of samples in VCF (default: 187)")
        sys.exit(1)
    
    chrom = sys.argv[1]
    vcf_file = sys.argv[2]
    annotation_file = sys.argv[3]
    n_samples = int(sys.argv[4]) if len(sys.argv) == 5 else 187
    
    print(f"Loading annotated sites for {chrom}...")
    gsites = load_annotated_sites(annotation_file, chrom)
    print(f"  Loaded {len(gsites):,} annotated positions")
    
    print(f"Processing VCF file...")
    print(f"  Number of samples: {n_samples}")
    counts = process_vcf(vcf_file, gsites, chrom)
    print(f"  Found data for {len(counts)} genes")
    
    output_file = f"{chrom}.bygene.pi.txt"
    print(f"Writing output to {output_file}...")
    write_output(counts, chrom, output_file, n_samples)
    
    print("Done!")

if __name__ == "__main__":
    main()
