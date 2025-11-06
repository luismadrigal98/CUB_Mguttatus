#!/usr/bin/env python3
"""
Annotate each genomic position with:
- Gene ID
- Codon position (1, 2, 3)
- Degeneracy class (0-fold, 2-fold, 3-fold, 4-fold) - TRUE DEGENERACY
- Reference codon
- Amino acid

This matches the R classification scheme with accurate degeneracy testing.

Usage:
    python describe_gene_positions_by_degeneracy.py <chromosome> <gff3_file> <genome_fasta> <cds_fasta>

Output:
    <chromosome>.genic_bases.annotated.txt
    Format: Chr\tGene\tPosition\tBase\tCodon_Position\tDegeneracy\tRef_Codon\tAmino_Acid

Author: Luis Javier Madrigal-Roca & John K. Kelly
"""

import sys
from collections import defaultdict

# Genetic code for translation
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

def complement(base):
    """Return complement of a DNA base."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return comp.get(base.upper(), 'N')

def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    return ''.join(complement(b) for b in reversed(seq))

def translate_codon(codon):
    """Translate a codon to amino acid."""
    if len(codon) != 3:
        return 'X'
    return GENETIC_CODE.get(codon.upper(), 'X')

def calculate_degeneracy_level(codon, position):
    """
    Calculate TRUE degeneracy level by testing all 4 nucleotides.
    This matches the R function calculate_degeneracy_level.
    
    Args:
        codon: 3-letter codon string
        position: position within codon (0, 1, or 2 for 1st, 2nd, 3rd)
    
    Returns:
        Number of synonymous nucleotides at this position (1, 2, 3, or 4)
    """
    if len(codon) != 3 or position < 0 or position > 2:
        return 0
    
    # Get original amino acid
    original_aa = translate_codon(codon)
    if original_aa == 'X':
        return 0
    
    # Test all 4 possible nucleotides at this position
    nucleotides = ['A', 'T', 'G', 'C']
    synonymous_count = 0
    
    for nuc in nucleotides:
        # Create test codon
        test_codon = list(codon)
        test_codon[position] = nuc
        test_codon_str = ''.join(test_codon)
        
        # Translate
        test_aa = translate_codon(test_codon_str)
        
        # Count if synonymous
        if test_aa == original_aa:
            synonymous_count += 1
    
    return synonymous_count

def classify_degeneracy(codon, codon_pos):
    """
    Classify degeneracy using TRUE degeneracy testing.
    Matches R classification scheme.
    
    Args:
        codon: 3-letter codon string
        codon_pos: position within codon (0, 1, or 2 for 1st, 2nd, 3rd)
    
    Returns:
        "0-fold", "2-fold", "3-fold", or "4-fold"
    """
    degeneracy_level = calculate_degeneracy_level(codon, codon_pos)
    
    if degeneracy_level == 1:
        return "0-fold"
    elif degeneracy_level == 2:
        return "2-fold"
    elif degeneracy_level == 3:
        return "3-fold"
    elif degeneracy_level == 4:
        return "4-fold"
    else:
        return "unknown"

def load_genome_sequence(fasta_file, chrom):
    """Load chromosome sequence from genome FASTA."""
    seq = []
    in_target = False
    found_chroms = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Check if this is our chromosome
                header = line.strip()[1:].split()[0]
                found_chroms.append(header)
                in_target = (header == chrom)
                if in_target:
                    print(f"  Found chromosome: {header}", file=sys.stderr)
            elif in_target:
                seq.append(line.strip().upper())
    
    result = ''.join(seq)
    
    # Error checking
    if len(result) == 0:
        print(f"\nERROR: Chromosome '{chrom}' not found in {fasta_file}!", file=sys.stderr)
        print(f"Available chromosomes (first 5):", file=sys.stderr)
        for i, c in enumerate(found_chroms[:5]):
            print(f"  {c}", file=sys.stderr)
        if len(found_chroms) > 5:
            print(f"  ... and {len(found_chroms) - 5} more", file=sys.stderr)
        print(f"\nTip: Check that chromosome name matches exactly (case-sensitive)", file=sys.stderr)
        sys.exit(1)
    
    return result

def load_cds_sequences(cds_fasta):
    """Load CDS sequences from FASTA file."""
    cds_seqs = {}
    current_id = None
    current_seq = []
    
    with open(cds_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Save previous sequence
                if current_id:
                    cds_seqs[current_id] = ''.join(current_seq)
                
                # Parse header: >MgIM767.01G000100.1 pacid=64873285 ...
                header = line.strip()[1:].split()[0]
                # Extract gene ID (remove transcript suffix .1, .2, etc.)
                gene_id = '.'.join(header.split('.')[:2])  # MgIM767.01G000100
                # Simplify to match your format: 01G000100
                gene_id = gene_id.split('.')[1]  # 01G000100
                current_id = gene_id
                current_seq = []
            else:
                current_seq.append(line.strip().upper())
        
        # Save last sequence
        if current_id:
            cds_seqs[current_id] = ''.join(current_seq)
    
    return cds_seqs

def parse_gff3(gff3_file, chrom):
    """
    Parse GFF3 to extract gene features.
    
    Returns:
        genes: dict of gene_id -> {strand, cds_regions, utr5, utr3, introns}
    """
    genes = defaultdict(lambda: {
        'strand': '+',
        'cds': [],
        'utr5': [],
        'utr3': [],
        'gene_range': None
    })
    
    found_chroms = set()
    
    with open(gff3_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            
            chr_name, _, feature, start, end, _, strand, _, attrs = cols
            
            # Track all chromosomes seen
            found_chroms.add(chr_name)
            
            if chr_name != chrom:
                continue
            
            start = int(start)
            end = int(end)
            
            # Parse attributes
            attr_dict = {}
            for attr in attrs.split(';'):
                if '=' in attr:
                    key, val = attr.split('=', 1)
                    attr_dict[key] = val
            
            # Extract gene ID and check if primary transcript
            if 'Parent' in attr_dict:
                parent = attr_dict['Parent']
                # Only process primary transcript (.1)
                # Example: MgIM767.01G000100.1.v2.1 or MgIM767.01G000100.2.v2.1
                parts = parent.split('.')
                if len(parts) >= 3:
                    transcript_num = parts[2]
                    # Skip non-primary transcripts
                    if transcript_num != '1':
                        continue
                    gene_id = parts[1]  # Extract 01G000100
                else:
                    gene_id = parent.split('.')[1] if '.' in parent else parent
            elif 'ID' in attr_dict:
                # For gene features, just extract gene ID
                parts = attr_dict['ID'].split('.')
                if len(parts) >= 2:
                    gene_id = parts[1]
                else:
                    gene_id = attr_dict['ID']
            else:
                continue
            
            # Store features
            if feature == 'gene':
                genes[gene_id]['gene_range'] = (start, end)
                genes[gene_id]['strand'] = strand
            elif feature == 'CDS':
                genes[gene_id]['cds'].append((start, end))
                genes[gene_id]['strand'] = strand
            elif feature == 'five_prime_UTR':
                genes[gene_id]['utr5'].append((start, end))
            elif feature == 'three_prime_UTR':
                genes[gene_id]['utr3'].append((start, end))
    
    # Sort CDS regions by position
    for gene_id in genes:
        genes[gene_id]['cds'].sort()
        genes[gene_id]['utr5'].sort()
        genes[gene_id]['utr3'].sort()
    
    # Error checking
    if len(genes) == 0:
        print(f"\nERROR: No genes found for chromosome '{chrom}' in GFF3!", file=sys.stderr)
        print(f"Chromosomes found in GFF3 (first 10):", file=sys.stderr)
        for i, c in enumerate(sorted(found_chroms)[:10]):
            print(f"  {c}", file=sys.stderr)
        if len(found_chroms) > 10:
            print(f"  ... and {len(found_chroms) - 10} more", file=sys.stderr)
        print(f"\nTip: Check that chromosome name matches exactly (case-sensitive)", file=sys.stderr)
        sys.exit(1)
    
    return genes

def validate_cds_sequences(genes, genome_seq, cds_seqs, max_genes_to_check=20):
    """
    Validate CDS sequences by reconstructing from genome and comparing to provided CDS FASTA.
    
    Args:
        genes: Gene annotations from GFF3
        genome_seq: Chromosome sequence
        cds_seqs: CDS sequences from FASTA
        max_genes_to_check: Maximum number of genes to validate (for speed)
    
    Returns:
        validation_results: dict with validation statistics
    """
    print("\n=== CDS Validation ===")
    print(f"Validating up to {max_genes_to_check} genes...")
    
    validation_results = {
        'total_checked': 0,
        'perfect_matches': 0,
        'length_mismatches': 0,
        'sequence_mismatches': 0,
        'missing_in_fasta': 0,
        'warnings': []
    }
    
    genes_checked = 0
    for gene_id, gene_info in genes.items():
        if genes_checked >= max_genes_to_check:
            break
        
        cds_regions = gene_info['cds']
        strand = gene_info['strand']
        
        if not cds_regions:
            continue
        
        # Check if gene exists in CDS FASTA
        if gene_id not in cds_seqs:
            validation_results['missing_in_fasta'] += 1
            continue
        
        genes_checked += 1
        validation_results['total_checked'] += 1
        
        # Reconstruct CDS from genome
        reconstructed_cds = []
        for start, end in cds_regions:
            # Extract region (convert to 0-indexed)
            region_seq = genome_seq[start-1:end]
            reconstructed_cds.append(region_seq)
        
        reconstructed_cds = ''.join(reconstructed_cds)
        
        # Reverse complement if on minus strand
        if strand == '-':
            reconstructed_cds = reverse_complement(reconstructed_cds)
        
        # Get provided CDS
        provided_cds = cds_seqs[gene_id]
        
        # Compare
        if len(reconstructed_cds) != len(provided_cds):
            validation_results['length_mismatches'] += 1
            validation_results['warnings'].append(
                f"Gene {gene_id}: Length mismatch (reconstructed={len(reconstructed_cds)}, provided={len(provided_cds)})"
            )
        elif reconstructed_cds != provided_cds:
            validation_results['sequence_mismatches'] += 1
            # Find first difference
            for i, (r, p) in enumerate(zip(reconstructed_cds, provided_cds)):
                if r != p:
                    validation_results['warnings'].append(
                        f"Gene {gene_id}: Sequence mismatch at position {i} (reconstructed={r}, provided={p})"
                    )
                    break
        else:
            validation_results['perfect_matches'] += 1
    
    # Print summary
    print(f"\nValidation Results:")
    print(f"  Genes checked: {validation_results['total_checked']}")
    print(f"  Perfect matches: {validation_results['perfect_matches']}")
    print(f"  Length mismatches: {validation_results['length_mismatches']}")
    print(f"  Sequence mismatches: {validation_results['sequence_mismatches']}")
    print(f"  Missing in CDS FASTA: {validation_results['missing_in_fasta']}")
    
    if validation_results['warnings']:
        print(f"\nWarnings (showing first 10):")
        for warning in validation_results['warnings'][:10]:
            print(f"  {warning}")
    
    # Assess overall quality
    if validation_results['total_checked'] > 0:
        match_rate = validation_results['perfect_matches'] / validation_results['total_checked']
        print(f"\nMatch rate: {match_rate:.1%}")
        
        if match_rate >= 0.95:
            print("✓ CDS sequences are highly consistent with genome annotation")
        elif match_rate >= 0.80:
            print("⚠ Some CDS mismatches detected - proceeding with caution")
        else:
            print("✗ WARNING: Many CDS mismatches detected!")
            print("  This may indicate:")
            print("  - Genome and CDS FASTA are from different annotation versions")
            print("  - GFF3 and CDS FASTA gene IDs don't match properly")
            print("  - Errors in the reconstruction logic")
            print("  Recommend checking your input files carefully.")
    
    print("=" * 50 + "\n")
    
    return validation_results

def annotate_positions(chrom, genes, genome_seq, cds_seqs):
    """
    Annotate each position in coding regions with degeneracy class.
    
    Returns:
        annotations: list of (chr, gene, pos, base, codon_position, degeneracy, ref_codon, amino_acid)
    """
    annotations = []
    
    for gene_id, gene_info in genes.items():
        strand = gene_info['strand']
        cds_regions = gene_info['cds']
        
        if not cds_regions:
            continue
        
        # Get CDS sequence
        if gene_id not in cds_seqs:
            print(f"Warning: Gene {gene_id} not found in CDS FASTA", file=sys.stderr)
            continue
        
        cds_seq = cds_seqs[gene_id]
        
        # Build position -> (codon, codon_position) mapping
        cds_positions = []
        for start, end in cds_regions:
            for pos in range(start, end + 1):
                cds_positions.append(pos)
        
        # For minus strand, reverse positions so CDS[0] maps to highest genomic coordinate
        # This is necessary because:
        # 1. CDS FASTA sequences are already in gene orientation (5'→3')
        # 2. Minus strand genes: START codon (CDS[0:3]) is at highest genomic coordinates
        # 3. After reversal: cds_positions[0] = highest coordinate = START codon position
        if strand == '-':
            cds_positions = cds_positions[::-1]
        
        # Check if CDS length matches
        if len(cds_positions) != len(cds_seq):
            # This is expected for genes with overlapping isoforms or phase issues
            # Skip these genes - they're not suitable for codon analysis
            if len(cds_positions) == len(cds_seq) * 2:
                # Likely double-counting due to isoforms in GFF3
                print(f"Info: Gene {gene_id} skipped - likely has multiple isoforms (genomic={len(cds_positions)}, fasta={len(cds_seq)})", file=sys.stderr)
            else:
                print(f"Warning: Gene {gene_id} CDS length mismatch: genomic={len(cds_positions)}, fasta={len(cds_seq)}", file=sys.stderr)
            continue
        
        # Annotate each position
        for i, genomic_pos in enumerate(cds_positions):
            codon_idx = i // 3
            codon_pos = i % 3
            
            # Get codon from CDS sequence
            codon_start = codon_idx * 3
            codon_end = codon_start + 3
            
            if codon_end > len(cds_seq):
                break
            
            codon = cds_seq[codon_start:codon_end]
            
            # Get genomic base from reference genome (always from + strand)
            # This matches what VCF reports
            genomic_base = genome_seq[genomic_pos - 1]  # Convert to 0-indexed
            
            # Note: genomic_base is kept as reference strand (for VCF matching)
            # The codon from CDS FASTA is already in gene orientation (5'→3')
            
            # Classify degeneracy using TRUE degeneracy testing
            degeneracy = classify_degeneracy(codon, codon_pos)
            
            # Get amino acid
            amino_acid = translate_codon(codon)
            
            # 1-indexed codon position for output
            codon_position_label = codon_pos + 1
            
            annotations.append((chrom, gene_id, genomic_pos, genomic_base, 
                              codon_position_label, degeneracy, codon, amino_acid))
    
    return annotations

def main():
    if len(sys.argv) != 5:
        print("ERROR: Incorrect number of arguments", file=sys.stderr)
        print("\nUsage: python describe_gene_positions_by_degeneracy.py <chromosome> <gff3> <genome_fasta> <cds_fasta>", file=sys.stderr)
        print("\nExample:", file=sys.stderr)
        print("  python describe_gene_positions_by_degeneracy.py Chr_01 genes.gff3 genome.fa cds.fa", file=sys.stderr)
        print("\nIMPORTANT:", file=sys.stderr)
        print("  - Process ONE chromosome at a time", file=sys.stderr)
        print("  - Use exact chromosome name from FASTA (case-sensitive)", file=sys.stderr)
        print("  - For multiple chromosomes, run this script multiple times or use a loop", file=sys.stderr)
        sys.exit(1)
    
    chrom = sys.argv[1]
    gff3_file = sys.argv[2]
    genome_fasta = sys.argv[3]
    cds_fasta = sys.argv[4]
    
    # Check for spaces in chromosome name (common error)
    if ' ' in chrom:
        print(f"ERROR: Chromosome name contains spaces: '{chrom}'", file=sys.stderr)
        print("This script processes ONE chromosome at a time.", file=sys.stderr)
        print("Did you mean to run this in a loop?", file=sys.stderr)
        print("\nExample loop:", file=sys.stderr)
        print('  for CHR in Chr_01 Chr_02 Chr_03; do', file=sys.stderr)
        print('    python describe_gene_positions_by_degeneracy.py $CHR genes.gff3 genome.fa cds.fa', file=sys.stderr)
        print('  done', file=sys.stderr)
        sys.exit(1)
    
    # Check if files exist
    import os
    for filepath, name in [(gff3_file, "GFF3"), (genome_fasta, "Genome FASTA"), (cds_fasta, "CDS FASTA")]:
        if not os.path.exists(filepath):
            print(f"ERROR: {name} file not found: {filepath}", file=sys.stderr)
            sys.exit(1)
    
    print(f"Processing chromosome {chrom}...")
    
    # Load data
    print("Loading genome sequence...")
    genome_seq = load_genome_sequence(genome_fasta, chrom)
    print(f"  Loaded {len(genome_seq)} bp")
    
    print("Loading CDS sequences...")
    cds_seqs = load_cds_sequences(cds_fasta)
    print(f"  Loaded {len(cds_seqs)} genes")
    
    print("Parsing GFF3...")
    genes = parse_gff3(gff3_file, chrom)
    print(f"  Found {len(genes)} genes")
    
    # Validate CDS sequences
    validation_results = validate_cds_sequences(genes, genome_seq, cds_seqs, max_genes_to_check=20)
    
    print("Annotating positions...")
    annotations = annotate_positions(chrom, genes, genome_seq, cds_seqs)
    print(f"  Annotated {len(annotations)} positions")
    
    # Write output
    output_file = f"{chrom}.genic_bases.annotated.txt"
    print(f"Writing to {output_file}...")
    
    with open(output_file, 'w') as out:
        out.write("Chr\tGene\tPosition\tBase\tCodon_Position\tDegeneracy\tRef_Codon\tAmino_Acid\n")
        for chr_name, gene_id, pos, base, codon_pos, degeneracy, ref_codon, amino_acid in sorted(annotations, key=lambda x: x[2]):
            out.write(f"{chr_name}\t{gene_id}\t{pos}\t{base}\t{codon_pos}\t{degeneracy}\t{ref_codon}\t{amino_acid}\n")
    
    print(f"Done! Annotated {len(annotations)} positions for {len(genes)} genes.")

if __name__ == "__main__":
    main()
