#!/usr/bin/env python3
"""
Calculate the frequency of preferred vs non-preferred codons across lines.
This script extracts codon genotypes per line to test if preferred codons
show different allele frequency spectra (evidence of selection).

Input:
    - VCF file with variant AND invariant sites
    - <chrom>.genic_bases.annotated.txt
    - preferred_codons.txt (list of optimal codons from CAI analysis)
    - GFF3 file
    - Genome FASTA
    - CDS FASTA

Output:
    - <chrom>.codon_frequencies.txt: Codon genotypes per gene per position
    - <chrom>.preferred_codon_stats.txt: Summary statistics

Format of codon_frequencies.txt:
    Gene  Codon_Pos  AA  Ref_Codon  Is_Preferred  Line1_Codon  Line2_Codon  ...  Freq_Preferred

Author: Luis Javier Madrigal-Roca & John K. Kelly
"""

import sys
from collections import defaultdict, Counter

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

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
            if len(codon) == 3:
                preferred.add(codon)
    
    return preferred

def complement(base):
    """Return complement of a DNA base."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return comp.get(base.upper(), 'N')

def load_gene_cds_info(gff3_file, cds_fasta, chrom):
    """
    Load gene CDS information.
    
    Returns:
        gene_info: dict of gene_id -> {
            'strand': +/-,
            'cds_regions': [(start, end), ...],
            'cds_seq': DNA sequence
        }
    """
    from collections import defaultdict
    
    # Parse GFF3
    genes = defaultdict(lambda: {
        'strand': '+',
        'cds_regions': [],
        'cds_seq': ''
    })
    
    with open(gff3_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            
            chr_name, _, feature, start, end, _, strand, _, attrs = cols
            
            if chr_name != chrom or feature != 'CDS':
                continue
            
            # Extract gene ID from Parent attribute
            attr_dict = {}
            for attr in attrs.split(';'):
                if '=' in attr:
                    key, val = attr.split('=', 1)
                    attr_dict[key] = val
            
            if 'Parent' not in attr_dict:
                continue
            
            parent = attr_dict['Parent']
            
            # Only process primary transcript (.1)
            # Example: MgIM767.01G000100.1.v2.1 or MgIM767.01G000100.2.v2.1
            parts = parent.split('.')
            if len(parts) >= 3:
                transcript_num = parts[2]
                # Skip non-primary transcripts
                if transcript_num != '1':
                    continue
                gene_id = parts[1]  # Extract "01G000100"
            else:
                gene_id = parent.split('.')[1] if '.' in parent else parent
            
            genes[gene_id]['cds_regions'].append((int(start), int(end)))
            genes[gene_id]['strand'] = strand
    
    # Sort CDS regions
    for gene_id in genes:
        genes[gene_id]['cds_regions'].sort()
    
    # Load CDS sequences
    with open(cds_fasta, 'r') as f:
        current_id = None
        current_seq = []
        
        for line in f:
            if line.startswith('>'):
                # Save previous
                if current_id:
                    genes[current_id]['cds_seq'] = ''.join(current_seq)
                
                # Parse new header
                header = line.strip()[1:].split()[0]
                gene_id = header.split('.')[1] if '.' in header else header
                current_id = gene_id
                current_seq = []
            else:
                current_seq.append(line.strip().upper())
        
        # Save last
        if current_id:
            genes[current_id]['cds_seq'] = ''.join(current_seq)
    
    return genes

def get_codon_from_positions(genomic_pos, gene_info):
    """
    Given a genomic position, determine which codon it belongs to.
    
    Returns:
        (codon_number, codon_position, codon_start_pos, codon_end_pos)
        or None if position not in CDS
    """
    cds_regions = gene_info['cds_regions']
    
    # Find cumulative CDS position
    cds_pos = 0
    for start, end in cds_regions:
        if start <= genomic_pos <= end:
            cds_pos += (genomic_pos - start)
            break
        elif genomic_pos < start:
            return None
        else:
            cds_pos += (end - start + 1)
    else:
        return None
    
    # Determine codon
    codon_number = cds_pos // 3
    codon_position = cds_pos % 3
    
    # Find genomic positions of all 3 bases in this codon
    codon_cds_start = codon_number * 3
    codon_positions_genomic = []
    
    cumulative = 0
    for start, end in cds_regions:
        region_length = end - start + 1
        
        for i in range(3):
            target_cds_pos = codon_cds_start + i
            
            if cumulative <= target_cds_pos < cumulative + region_length:
                genomic_pos_i = start + (target_cds_pos - cumulative)
                codon_positions_genomic.append(genomic_pos_i)
        
        cumulative += region_length
        
        if len(codon_positions_genomic) == 3:
            break
    
    if len(codon_positions_genomic) != 3:
        return None
    
    return codon_number, codon_position, codon_positions_genomic

def extract_codon_genotypes(vcf_file, gene_info, preferred_codons, chrom, n_samples):
    """
    Extract codon genotypes for each line across all codon positions.
    
    Args:
        vcf_file: Path to VCF file
        gene_info: Dict of gene information including strand
        preferred_codons: Set of preferred codons
        chrom: Chromosome name
        n_samples: Number of samples/lines in VCF (starts at column 10)
    
    Returns:
        codon_data: dict of (gene_id, codon_number) -> {
            'ref_codon': str,
            'aa': str,
            'is_preferred': bool,
            'line_codons': [codon1, codon2, ...],
            'genomic_positions': [pos1, pos2, pos3]
        }
    """
    codon_data = defaultdict(lambda: {
        'ref_codon': '',
        'aa': '',
        'is_preferred': False,
        'line_codons': [],
        'genomic_positions': []
    })
    
    # Build position -> gene mapping
    pos_to_gene = {}
    for gene_id, info in gene_info.items():
        for start, end in info['cds_regions']:
            for pos in range(start, end + 1):
                pos_to_gene[pos] = gene_id
    
    # Store variants per position
    variants = {}  # position -> {line_idx: alt_base}
    missing = defaultdict(set)  # position -> set of line_idx with missing data
    
    print(f"Scanning VCF for variants...")
    line_count = 0
    
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            line_count += 1
            if line_count % 1000000 == 0:
                print(f"  {line_count:,} sites...", file=sys.stderr)
            
            cols = line.strip().split('\t')
            if len(cols) < 10:
                continue
            
            vcf_chrom = cols[0]
            if vcf_chrom != chrom:
                continue
            
            pos = int(cols[1])
            ref = cols[3]
            alt = cols[4]
            
            if alt == '.':
                continue  # Invariant
            
            if pos not in pos_to_gene:
                continue  # Not in CDS
            
            # Parse genotypes
            variants[pos] = {}
            for i in range(9, len(cols)):
                line_idx = i - 9
                gt_field = cols[i]
                parts = gt_field.split(':')
                
                if len(parts) < 3:
                    continue
                
                gt = parts[0]
                ad = parts[2]
                
                try:
                    # Handle both formats:
                    # Variant sites: "12,5" (ref,alt)
                    # Invariant sites: "2" (ref only)
                    if ',' in ad:
                        ref_count, alt_count = map(int, ad.split(','))
                    else:
                        # Invariant site - only ref depth
                        ref_count = int(ad)
                        alt_count = 0
                except (ValueError, IndexError):
                    continue
                
                # Track if no coverage (missing data)
                if ref_count == 0 and alt_count == 0:
                    missing[pos].add(line_idx)
                    continue
                
                # Determine genotype (simple: use majority allele with depth threshold)
                if alt_count > 5 * ref_count and alt_count > 0:
                    variants[pos][line_idx] = alt
                elif ref_count > 5 * alt_count and ref_count > 0:
                    variants[pos][line_idx] = ref
    
    print(f"  Found {len(variants):,} polymorphic sites in CDS")
    
    # Build codon genotypes
    print("Building codon genotypes...")
    
    for gene_id, info in gene_info.items():
        if not info['cds_seq']:
            continue
        
        cds_seq = info['cds_seq']
        strand = info['strand']
        n_codons = len(cds_seq) // 3
        
        # Build ALL CDS positions first (same logic as describe_gene_positions_by_degeneracy.py)
        all_cds_positions = []
        for start, end in info['cds_regions']:
            for pos in range(start, end + 1):
                all_cds_positions.append(pos)
        
        # For minus strand, reverse positions so CDS[0] maps to highest genomic coordinate
        if strand == '-':
            all_cds_positions = all_cds_positions[::-1]
        
        # Now extract codon positions from the (possibly reversed) list
        for codon_idx in range(n_codons):
            # Get reference codon
            ref_codon = cds_seq[codon_idx * 3:(codon_idx + 1) * 3]
            
            if len(ref_codon) != 3:
                continue
            
            aa = GENETIC_CODE.get(ref_codon, 'X')
            
            # Get genomic positions for this codon from the pre-built list
            cds_pos_start = codon_idx * 3
            if cds_pos_start + 3 > len(all_cds_positions):
                continue
            
            genomic_positions = all_cds_positions[cds_pos_start:cds_pos_start + 3]
            
            if len(genomic_positions) != 3:
                continue
            
            # Build codon genotypes for each line
            line_codons = []
            
            for line_idx in range(n_samples):
                # Check if ANY position in this codon is missing data for this sample
                has_missing = any(gpos in missing and line_idx in missing[gpos] 
                                 for gpos in genomic_positions)
                
                if has_missing:
                    # Skip this sample entirely - don't create fake reference codon
                    line_codons.append('NNN')  # Mark as missing
                    continue
                
                codon_bases = []
                
                for pos_idx, gpos in enumerate(genomic_positions):
                    if gpos in variants and line_idx in variants[gpos]:
                        # Base from VCF (always reference strand)
                        base = variants[gpos][line_idx]
                        
                        # If reverse strand, complement to get gene orientation
                        if strand == '-':
                            base = complement(base)
                    else:
                        # Use reference codon base (already in gene orientation)
                        base = ref_codon[pos_idx]
                    
                    codon_bases.append(base)
                
                # Codon bases are now all in gene orientation (5'→3')
                # No need to reverse - they're already correct!
                line_codon = ''.join(codon_bases)
                line_codons.append(line_codon)
            
            # Store
            key = (gene_id, codon_idx)
            codon_data[key] = {
                'ref_codon': ref_codon,
                'aa': aa,
                'is_preferred': ref_codon in preferred_codons,
                'line_codons': line_codons,
                'genomic_positions': genomic_positions
            }
    
    return codon_data

def calculate_preferred_frequency_spectrum(codon_data, preferred_codons):
    """
    Calculate allele frequency spectrum for preferred vs non-preferred codons.
    
    Returns:
        stats: dict with summary statistics
    """
    preferred_freqs = []
    nonpreferred_freqs = []
    
    for (gene_id, codon_idx), data in codon_data.items():
        aa = data['aa']
        ref_codon = data['ref_codon']
        line_codons = data['line_codons']
        
        if aa == '*' or aa == 'X':
            continue
        
        # Get all synonymous codons for this amino acid
        synonymous = [c for c, a in GENETIC_CODE.items() if a == aa]
        
        if len(synonymous) <= 1:
            continue  # Not synonymous
        
        # Count preferred vs non-preferred
        codon_counts = Counter(line_codons)
        
        # Remove missing data marker
        if 'NNN' in codon_counts:
            del codon_counts['NNN']
        
        total_lines = sum(codon_counts.values())
        
        if total_lines == 0:
            continue  # All samples missing
        
        for codon, count in codon_counts.items():
            if codon not in GENETIC_CODE or GENETIC_CODE[codon] != aa:
                continue
            
            freq = float(count) / total_lines
            
            if codon in preferred_codons:
                preferred_freqs.append(freq)
            else:
                nonpreferred_freqs.append(freq)
    
    return {
        'preferred_freqs': preferred_freqs,
        'nonpreferred_freqs': nonpreferred_freqs
    }

def write_output(codon_data, output_file):
    """Write codon genotype data to file."""
    
    with open(output_file, 'w') as out:
        # Header
        out.write("Gene\tCodon_Pos\tAA\tRef_Codon\tIs_Preferred\t")
        out.write("Genomic_Positions\t")
        
        # Codon counts
        out.write("Codon_Variants\tFrequencies\n")
        
        for (gene_id, codon_idx), data in sorted(codon_data.items()):
            ref_codon = data['ref_codon']
            aa = data['aa']
            is_pref = data['is_preferred']
            positions = ','.join(map(str, data['genomic_positions']))
            line_codons = data['line_codons']
            
            # Count codon variants
            codon_counts = Counter(line_codons)
            total = len(line_codons)
            
            codon_str = ';'.join(f"{codon}:{count}" for codon, count in codon_counts.most_common())
            freq_str = ';'.join(f"{codon}:{count/total:.3f}" for codon, count in codon_counts.most_common())
            
            out.write(f"{gene_id}\t{codon_idx}\t{aa}\t{ref_codon}\t{is_pref}\t")
            out.write(f"{positions}\t")
            out.write(f"{codon_str}\t{freq_str}\n")

def main():
    if len(sys.argv) not in [7, 8]:
        print("Usage: python calculate_freq_preferred_codon.py <chrom> <vcf> <gff3> <cds_fasta> <genome_fasta> <preferred_codons> [n_samples]")
        print("  n_samples: Number of samples in VCF (default: 187)")
        sys.exit(1)
    
    chrom = sys.argv[1]
    vcf_file = sys.argv[2]
    gff3_file = sys.argv[3]
    cds_fasta = sys.argv[4]
    # genome_fasta = sys.argv[5]  # Not currently used
    preferred_file = sys.argv[6]
    n_samples = int(sys.argv[7]) if len(sys.argv) == 8 else 187
    
    print(f"Configuration:")
    print(f"  Chromosome: {chrom}")
    print(f"  Number of samples: {n_samples}")
    
    print("Loading preferred codons...")
    preferred_codons = load_preferred_codons(preferred_file)
    print(f"  {len(preferred_codons)} preferred codons")
    
    print("Loading gene CDS information...")
    gene_info = load_gene_cds_info(gff3_file, cds_fasta, chrom)
    print(f"  {len(gene_info)} genes")
    
    print("Extracting codon genotypes from VCF...")
    codon_data = extract_codon_genotypes(vcf_file, gene_info, preferred_codons, chrom, n_samples)
    print(f"  {len(codon_data)} codon positions")
    
    output_file = f"{chrom}.codon_frequencies.txt"
    print(f"Writing output to {output_file}...")
    write_output(codon_data, output_file)
    
    print("Calculating preferred codon frequency spectrum...")
    stats = calculate_preferred_frequency_spectrum(codon_data, preferred_codons)
    
    print("\n=== Summary Statistics ===")
    print(f"Preferred codon observations: {len(stats['preferred_freqs'])}")
    print(f"Non-preferred codon observations: {len(stats['nonpreferred_freqs'])}")
    
    if stats['preferred_freqs']:
        if HAS_NUMPY:
            print(f"Mean frequency (preferred): {np.mean(stats['preferred_freqs']):.4f}")
        else:
            mean_pref = sum(stats['preferred_freqs']) / len(stats['preferred_freqs'])
            print(f"Mean frequency (preferred): {mean_pref:.4f}")
    
    if stats['nonpreferred_freqs']:
        if HAS_NUMPY:
            print(f"Mean frequency (non-preferred): {np.mean(stats['nonpreferred_freqs']):.4f}")
        else:
            mean_nonpref = sum(stats['nonpreferred_freqs']) / len(stats['nonpreferred_freqs'])
            print(f"Mean frequency (non-preferred): {mean_nonpref:.4f}")
    
    print("\nDone!")

if __name__ == "__main__":
    main()
