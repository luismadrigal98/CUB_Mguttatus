#!/usr/bin/env python3
"""
Calculate site frequency spectrum for preferred vs non-preferred states
at each degeneracy class.

This script analyzes allele frequencies at nucleotide sites, categorized by:
1. Degeneracy class (0-fold, 2-fold, 3-fold, 4-fold)
2. Whether the site changes preferred→non-preferred or vice versa

Input:
    - VCF file with variant sites
    - <chrom>.genic_bases.annotated.txt (from describe_gene_positions_by_degeneracy.py)
    - preferred_codons.txt (list of optimal codons from CAI analysis)

Output:
    - <chrom>.site_freq_by_preference.txt: Site-level frequencies per gene
    
Format:
    Gene  Degeneracy  Sites_Preferred_Ref  Sites_NonPref_Ref  
    Poly_Pref_to_NonPref  Poly_NonPref_to_Pref  
    Pi_Pref_Ref  Pi_NonPref_Ref

Author: Luis Javier Madrigal-Roca & John K. Kelly
"""

import sys
from collections import defaultdict

# Genetic code
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def load_preferred_codons(preferred_file):
    """Load list of preferred codons from file."""
    preferred = set()
    
    with open(preferred_file, 'r') as f:
        for line in f:
            codon = line.strip().upper()
            if len(codon) == 3 and codon in GENETIC_CODE:
                preferred.add(codon)
    
    print(f"  Loaded {len(preferred)} preferred codons:")
    # Group by amino acid
    aa_to_pref = defaultdict(list)
    for codon in preferred:
        aa = GENETIC_CODE[codon]
        aa_to_pref[aa].append(codon)
    
    for aa, codons in sorted(aa_to_pref.items()):
        print(f"    {aa}: {', '.join(sorted(codons))}")
    
    return preferred

def complement(base):
    """Return complement of a DNA base."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return comp.get(base.upper(), 'N')

def load_annotated_sites(annotation_file, chrom):
    """
    Load position annotations from describe_gene_positions_by_degeneracy.py output.
    
    Returns:
        sites: dict of position -> {
            'gene': str,
            'base': str,
            'codon_pos': int,
            'degeneracy': str,
            'ref_codon': str,
            'amino_acid': str,
            'strand': str
        }
    """
    sites = {}
    
    with open(annotation_file, 'r') as f:
        header = f.readline().strip().split('\t')
        
        # Check if strand column exists (for backwards compatibility)
        has_strand = 'Strand' in header
        if not has_strand:
            print("  WARNING: Annotation file lacks Strand column. Minus strand genes may be miscategorized.")
            print("           Re-run describe_gene_positions_by_degeneracy.py to regenerate annotations.")
        
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 8:
                continue
            
            chr_name = cols[0]
            gene_id = cols[1]
            pos = cols[2]
            base = cols[3]
            codon_pos = cols[4]
            degeneracy = cols[5]
            ref_codon = cols[6]
            amino_acid = cols[7]
            strand = cols[8] if len(cols) > 8 else '+'  # Default to + if missing
            
            if chr_name == chrom:
                sites[int(pos)] = {
                    'gene': gene_id,
                    'base': base,
                    'codon_pos': int(codon_pos),
                    'degeneracy': degeneracy,
                    'ref_codon': ref_codon,
                    'amino_acid': amino_acid,
                    'strand': strand
                }
    
    return sites

def is_preferred_state(codon, preferred_codons):
    """Check if a codon is in preferred state."""
    return codon in preferred_codons

def change_base_in_codon(codon, position, new_base):
    """
    Change base at position (0, 1, or 2) in codon.
    
    Args:
        codon: str, 3-letter codon
        position: int, 0=first, 1=second, 2=third
        new_base: str, new base to substitute
    
    Returns:
        str: modified codon
    """
    codon_list = list(codon)
    codon_list[position] = new_base
    return ''.join(codon_list)

def categorize_variant(site_info, ref_base, alt_base, preferred_codons):
    """
    Categorize a variant as preferred→non-preferred or vice versa.
    
    IMPORTANT: VCF reports variants on reference strand, but ref_codon is in
    gene orientation (5'→3'). For minus strand genes, we must complement the
    alt_base before substituting into the codon.
    
    Returns:
        str: 'pref_to_nonpref', 'nonpref_to_pref', 'pref_to_pref', 'nonpref_to_nonpref', or 'non_synonymous'
    """
    ref_codon = site_info['ref_codon']
    codon_pos = site_info['codon_pos'] - 1  # Convert to 0-indexed
    amino_acid = site_info['amino_acid']
    strand = site_info.get('strand', '+')
    
    # CRITICAL: Convert alt_base from reference strand to gene orientation
    # VCF always reports on + strand; for minus strand genes, complement the base
    if strand == '-':
        alt_base_gene = complement(alt_base)
    else:
        alt_base_gene = alt_base
    
    # Create alternate codon (now both are in gene orientation)
    alt_codon = change_base_in_codon(ref_codon, codon_pos, alt_base_gene)
    
    # Check if synonymous
    if alt_codon not in GENETIC_CODE:
        return 'invalid'
    
    alt_aa = GENETIC_CODE[alt_codon]
    
    if alt_aa != amino_acid:
        return 'non_synonymous'
    
    # Check preferred status
    ref_is_pref = ref_codon in preferred_codons
    alt_is_pref = alt_codon in preferred_codons
    
    if ref_is_pref and not alt_is_pref:
        return 'pref_to_nonpref'
    elif not ref_is_pref and alt_is_pref:
        return 'nonpref_to_pref'
    elif ref_is_pref and alt_is_pref:
        return 'pref_to_pref'
    else:
        return 'nonpref_to_nonpref'

def calculate_pi_site(genotypes, min_depth_ratio=5):
    """
    Calculate nucleotide diversity (π) for a single polymorphic site.
    
    Same implementation as calculate_pi.py.
    
    Returns:
        (is_polymorphic, pi_value, n_samples)
    """
    ref_hom = 0
    alt_hom = 0
    
    for gt, ref_count, alt_count in genotypes:
        if gt == "0/0" and ref_count > min_depth_ratio * alt_count:
            ref_hom += 1
        elif gt == "1/1" and alt_count > min_depth_ratio * ref_count:
            alt_hom += 1
    
    if min(ref_hom, alt_hom) > 0:
        nx = float(ref_hom + alt_hom)
        px = float(ref_hom) / nx
        pi = 2.0 * nx * px * (1.0 - px) / (nx - 1.0)
        return True, pi, int(nx)
    else:
        return False, 0.0, 0

def parse_vcf_line(line):
    """
    Parse a single VCF line.
    
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
    
    genotypes = []
    for j in range(9, len(cols)):
        gt_field = cols[j]
        parts = gt_field.split(':')
        
        if len(parts) < 3:
            genotypes.append(("./.", 0, 0))
            continue
        
        gt = parts[0]
        
        try:
            ad_field = parts[2]
            ad_parts = ad_field.split(',')
            ref_count = int(ad_parts[0])
            alt_count = int(ad_parts[1]) if len(ad_parts) > 1 else 0
        except (ValueError, IndexError):
            ref_count, alt_count = 0, 0
        
        genotypes.append((gt, ref_count, alt_count))
    
    return chrom, pos, ref, alt, genotypes

def process_vcf(vcf_file, sites, preferred_codons, chrom):
    """
    Process VCF file and categorize sites by degeneracy and preference.
    
    Returns:
        counts: dict of gene -> {
            degeneracy_class: {
                'sites_pref_ref': int,
                'sites_nonpref_ref': int,
                'poly_pref_to_nonpref': int,
                'poly_nonpref_to_pref': int,
                'pi_pref_ref': float,
                'pi_nonpref_ref': float
            }
        }
    """
    counts = defaultdict(lambda: {
        '0-fold': {'sites_pref_ref': 0, 'sites_nonpref_ref': 0, 
                   'poly_pref_to_nonpref': 0, 'poly_nonpref_to_pref': 0,
                   'pi_pref_ref': 0.0, 'pi_nonpref_ref': 0.0},
        '2-fold': {'sites_pref_ref': 0, 'sites_nonpref_ref': 0,
                   'poly_pref_to_nonpref': 0, 'poly_nonpref_to_pref': 0,
                   'pi_pref_ref': 0.0, 'pi_nonpref_ref': 0.0},
        '3-fold': {'sites_pref_ref': 0, 'sites_nonpref_ref': 0,
                   'poly_pref_to_nonpref': 0, 'poly_nonpref_to_pref': 0,
                   'pi_pref_ref': 0.0, 'pi_nonpref_ref': 0.0},
        '4-fold': {'sites_pref_ref': 0, 'sites_nonpref_ref': 0,
                   'poly_pref_to_nonpref': 0, 'poly_nonpref_to_pref': 0,
                   'pi_pref_ref': 0.0, 'pi_nonpref_ref': 0.0}
    })
    
    print(f"Processing VCF for {chrom}...")
    
    # First pass: count invariant sites by reference state
    for pos, site_info in sites.items():
        gene = site_info['gene']
        degeneracy = site_info['degeneracy']
        ref_codon = site_info['ref_codon']
        
        # Check if reference codon is preferred
        ref_is_pref = ref_codon in preferred_codons
        
        if ref_is_pref:
            counts[gene][degeneracy]['sites_pref_ref'] += 1
        else:
            counts[gene][degeneracy]['sites_nonpref_ref'] += 1
    
    print(f"  Counted {len(sites):,} total sites by reference state")
    
    # Second pass: process polymorphic sites
    line_count = 0
    poly_count = 0
    
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            line_count += 1
            if line_count % 1000000 == 0:
                print(f"  Processed {line_count:,} VCF lines...", file=sys.stderr)
            
            parsed = parse_vcf_line(line)
            if not parsed:
                continue
            
            vcf_chrom, pos, ref, alt, genotypes = parsed
            
            if vcf_chrom != chrom:
                continue
            
            if alt == '.':
                continue  # Invariant
            
            if pos not in sites:
                continue  # Not in CDS
            
            site_info = sites[pos]
            gene = site_info['gene']
            degeneracy = site_info['degeneracy']
            
            # Calculate π for this site
            is_poly, pi_value, n_samples = calculate_pi_site(genotypes)
            
            if not is_poly:
                continue
            
            poly_count += 1
            
            # Categorize variant
            category = categorize_variant(site_info, ref, alt, preferred_codons)
            
            if category == 'pref_to_nonpref':
                counts[gene][degeneracy]['poly_pref_to_nonpref'] += 1
                counts[gene][degeneracy]['pi_pref_ref'] += pi_value
            elif category == 'nonpref_to_pref':
                counts[gene][degeneracy]['poly_nonpref_to_pref'] += 1
                counts[gene][degeneracy]['pi_nonpref_ref'] += pi_value
            elif category == 'pref_to_pref':
                # Stays preferred - count as preferred ref
                counts[gene][degeneracy]['pi_pref_ref'] += pi_value
            elif category == 'nonpref_to_nonpref':
                # Stays non-preferred - count as non-preferred ref
                counts[gene][degeneracy]['pi_nonpref_ref'] += pi_value
    
    print(f"  Total sites processed: {line_count:,}")
    print(f"  Polymorphic sites found: {poly_count:,}")
    
    return counts

def write_output(counts, chrom, output_file):
    """Write per-gene site frequency statistics to output file."""
    
    with open(output_file, 'w') as out:
        # Header
        out.write("Chr\tGene\tDegeneracy\t")
        out.write("Sites_Pref_Ref\tSites_NonPref_Ref\t")
        out.write("Poly_Pref_to_NonPref\tPoly_NonPref_to_Pref\t")
        out.write("Pi_Pref_Ref\tPi_NonPref_Ref\n")
        
        for gene_id in sorted(counts.keys()):
            gene_counts = counts[gene_id]
            
            for degeneracy in ['0-fold', '2-fold', '3-fold', '4-fold']:
                data = gene_counts[degeneracy]
                
                # Calculate average π by dividing by number of sites
                sites_pref = data['sites_pref_ref']
                sites_nonpref = data['sites_nonpref_ref']
                
                # Average π per site (divide sum by number of sites)
                avg_pi_pref = data['pi_pref_ref'] / sites_pref if sites_pref > 0 else 0.0
                avg_pi_nonpref = data['pi_nonpref_ref'] / sites_nonpref if sites_nonpref > 0 else 0.0
                
                out.write(f"{chrom}\t{gene_id}\t{degeneracy}\t")
                out.write(f"{sites_pref}\t{sites_nonpref}\t")
                out.write(f"{data['poly_pref_to_nonpref']}\t{data['poly_nonpref_to_pref']}\t")
                out.write(f"{avg_pi_pref:.6f}\t{avg_pi_nonpref:.6f}\n")

def main():
    if len(sys.argv) != 5:
        print("Usage: python calculate_site_freq_by_preference.py <chromosome> <vcf_file> <annotation_file> <preferred_codons>")
        print("Example: python calculate_site_freq_by_preference.py Chr_01 variants.vcf Chr_01.genic_bases.annotated.txt preferred_codons.txt")
        sys.exit(1)
    
    chrom = sys.argv[1]
    vcf_file = sys.argv[2]
    annotation_file = sys.argv[3]
    preferred_file = sys.argv[4]
    
    print(f"Loading preferred codons from {preferred_file}...")
    preferred_codons = load_preferred_codons(preferred_file)
    
    print(f"Loading annotated sites for {chrom}...")
    sites = load_annotated_sites(annotation_file, chrom)
    print(f"  Loaded {len(sites):,} annotated positions")
    
    print(f"Processing VCF file...")
    counts = process_vcf(vcf_file, sites, preferred_codons, chrom)
    print(f"  Found data for {len(counts)} genes")
    
    output_file = f"{chrom}.site_freq_by_preference.txt"
    print(f"Writing output to {output_file}...")
    write_output(counts, chrom, output_file)
    
    print("Done!")

if __name__ == "__main__":
    main()
