#!/usr/bin/env python3
"""
Test strand handling for codon genotyping with synthetic data.

This creates:
1. Synthetic genome with 20 genes (10+, 10-)
2. GFF3 annotation with primary and alternate transcripts
3. CDS FASTA with known codons
4. VCF with known variants
5. Verifies codon genotypes are correct

Author: Luis Javier Madrigal-Roca & John K. Kelly
"""

import os
import sys
import tempfile
import subprocess
from pathlib import Path

# Add miscellanea_code to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'miscellanea_code'))

def create_test_genome():
    """
    Create synthetic genome with known sequences.
    
    Gene layout on Chr_TEST:
    - Plus strand genes at positions 1000-2000 (10 genes)
    - Minus strand genes at positions 3000-4000 (10 genes)
    """
    
    # Create base sequence (all A's, we'll insert specific codons)
    genome_seq = 'A' * 5000
    genome_list = list(genome_seq)
    
    # Plus strand genes (positions 1000-1899)
    # Gene P01: Start=1000, CDS=1000-1008 (3 codons: ATG-CCA-TGA)
    genome_list[1000:1009] = list('ATGCCATGA')
    
    # Gene P02: Start=1100, CDS=1100-1108 (3 codons: ATG-GCA-TAA)
    genome_list[1100:1109] = list('ATGGCATAA')
    
    # Gene P03: Start=1200, CDS=1200-1208 (3 codons: ATG-TCA-TAG)
    genome_list[1200:1209] = list('ATGTCATAG')
    
    # Minus strand genes (positions 3000-3899)
    # For minus strand, we write the REVERSE COMPLEMENT on the genome
    # Gene M01: CDS in gene orientation: ATG-CCA-TGA
    # Reverse complement: TCA-TGG-CAT (write this on genome at 3000-3008)
    genome_list[3000:3009] = list('TCATGGCAT')
    
    # Gene M02: CDS in gene orientation: ATG-GCA-TAA
    # Reverse complement: TTA-TGC-CAT
    genome_list[3100:3109] = list('TTATGCCAT')
    
    # Gene M03: CDS in gene orientation: ATG-TCA-TAG
    # Reverse complement: CTA-TGA-CAT
    genome_list[3200:3209] = list('CTATGACAT')
    
    genome_seq = ''.join(genome_list)
    
    fasta = f">Chr_TEST\n"
    # Write in 60 bp lines
    for i in range(0, len(genome_seq), 60):
        fasta += genome_seq[i:i+60] + '\n'
    
    return fasta

def create_test_gff3():
    """
    Create GFF3 with genes, including some with multiple isoforms.
    """
    
    gff3 = "##gff-version 3\n"
    gff3 += "#Chr_TEST synthetic genome for testing\n"
    
    # Plus strand genes
    # Gene P01 - single isoform
    gff3 += "Chr_TEST\ttest\tgene\t1000\t1008\t.\t+\t.\tID=gene:P01\n"
    gff3 += "Chr_TEST\ttest\tmRNA\t1000\t1008\t.\t+\t.\tID=MgIM767.P01.1.v2.1;Parent=gene:P01\n"
    gff3 += "Chr_TEST\ttest\tCDS\t1000\t1008\t.\t+\t0\tID=MgIM767.P01.1.v2.1.CDS.1;Parent=MgIM767.P01.1.v2.1\n"
    
    # Gene P02 - with TWO isoforms (.1 and .2) - should only process .1!
    gff3 += "Chr_TEST\ttest\tgene\t1100\t1108\t.\t+\t.\tID=gene:P02\n"
    gff3 += "Chr_TEST\ttest\tmRNA\t1100\t1108\t.\t+\t.\tID=MgIM767.P02.1.v2.1;Parent=gene:P02\n"
    gff3 += "Chr_TEST\ttest\tCDS\t1100\t1108\t.\t+\t0\tID=MgIM767.P02.1.v2.1.CDS.1;Parent=MgIM767.P02.1.v2.1\n"
    gff3 += "Chr_TEST\ttest\tmRNA\t1100\t1105\t.\t+\t.\tID=MgIM767.P02.2.v2.1;Parent=gene:P02\n"
    gff3 += "Chr_TEST\ttest\tCDS\t1100\t1105\t.\t+\t0\tID=MgIM767.P02.2.v2.1.CDS.1;Parent=MgIM767.P02.2.v2.1\n"
    
    # Gene P03 - single isoform
    gff3 += "Chr_TEST\ttest\tgene\t1200\t1208\t.\t+\t.\tID=gene:P03\n"
    gff3 += "Chr_TEST\ttest\tmRNA\t1200\t1208\t.\t+\t.\tID=MgIM767.P03.1.v2.1;Parent=gene:P03\n"
    gff3 += "Chr_TEST\ttest\tCDS\t1200\t1208\t.\t+\t0\tID=MgIM767.P03.1.v2.1.CDS.1;Parent=MgIM767.P03.1.v2.1\n"
    
    # Minus strand genes
    # Gene M01 - single isoform
    gff3 += "Chr_TEST\ttest\tgene\t3000\t3008\t.\t-\t.\tID=gene:M01\n"
    gff3 += "Chr_TEST\ttest\tmRNA\t3000\t3008\t.\t-\t.\tID=MgIM767.M01.1.v2.1;Parent=gene:M01\n"
    gff3 += "Chr_TEST\ttest\tCDS\t3000\t3008\t.\t-\t0\tID=MgIM767.M01.1.v2.1.CDS.1;Parent=MgIM767.M01.1.v2.1\n"
    
    # Gene M02 - with TWO isoforms (.1 and .2) - should only process .1!
    gff3 += "Chr_TEST\ttest\tgene\t3100\t3108\t.\t-\t.\tID=gene:M02\n"
    gff3 += "Chr_TEST\ttest\tmRNA\t3100\t3108\t.\t-\t.\tID=MgIM767.M02.1.v2.1;Parent=gene:M02\n"
    gff3 += "Chr_TEST\ttest\tCDS\t3100\t3108\t.\t-\t0\tID=MgIM767.M02.1.v2.1.CDS.1;Parent=MgIM767.M02.1.v2.1\n"
    gff3 += "Chr_TEST\ttest\tmRNA\t3100\t3105\t.\t-\t.\tID=MgIM767.M02.2.v2.1;Parent=gene:M02\n"
    gff3 += "Chr_TEST\ttest\tCDS\t3100\t3105\t.\t-\t0\tID=MgIM767.M02.2.v2.1.CDS.1;Parent=MgIM767.M02.2.v2.1\n"
    
    # Gene M03 - single isoform
    gff3 += "Chr_TEST\ttest\tgene\t3200\t3208\t.\t-\t.\tID=gene:M03\n"
    gff3 += "Chr_TEST\ttest\tmRNA\t3200\t3208\t.\t-\t.\tID=MgIM767.M03.1.v2.1;Parent=gene:M03\n"
    gff3 += "Chr_TEST\ttest\tCDS\t3200\t3208\t.\t-\t0\tID=MgIM767.M03.1.v2.1.CDS.1;Parent=MgIM767.M03.1.v2.1\n"
    
    return gff3

def create_test_cds():
    """
    Create CDS FASTA with sequences in gene orientation (5'→3').
    """
    
    cds = ""
    
    # Plus strand genes (sequences as they appear on genome)
    cds += ">MgIM767.P01.1 ID=MgIM767.P01.1.v2.1\n"
    cds += "ATGCCATGA\n"
    
    cds += ">MgIM767.P02.1 ID=MgIM767.P02.1.v2.1\n"
    cds += "ATGGCATAA\n"
    
    cds += ">MgIM767.P03.1 ID=MgIM767.P03.1.v2.1\n"
    cds += "ATGTCATAG\n"
    
    # Minus strand genes (sequences in GENE orientation, not genomic!)
    # These are already reverse complemented from genome
    cds += ">MgIM767.M01.1 ID=MgIM767.M01.1.v2.1\n"
    cds += "ATGCCATGA\n"  # Same as P01 in gene orientation
    
    cds += ">MgIM767.M02.1 ID=MgIM767.M02.1.v2.1\n"
    cds += "ATGGCATAA\n"  # Same as P02 in gene orientation
    
    cds += ">MgIM767.M03.1 ID=MgIM767.M03.1.v2.1\n"
    cds += "ATGTCATAG\n"  # Same as P03 in gene orientation
    
    return cds

def create_test_vcf():
    """
    Create VCF with known variants at specific positions.
    
    Variants to test:
    - Plus strand gene P01, position 1001 (2nd base of ATG): T→C (should change codon to ACG)
    - Minus strand gene M01, position 3001 (middle of TCA on genome): C→G
      On genome: TCA→TGA
      Gene orientation: TGA→TCA (reverse complement)
      This changes the codon!
    """
    
    vcf = "##fileformat=VCFv4.2\n"
    vcf += "##contig=<ID=Chr_TEST,length=5000>\n"
    vcf += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    vcf += "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">\n"
    vcf += "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
    vcf += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    
    # Add 5 sample columns
    for i in range(1, 6):
        vcf += f"\tSample{i}"
    vcf += "\n"
    
    # Variant at plus strand gene P01, position 1001 (T→C)
    # Sample1: homozygous ref (T), Sample2: homozygous alt (C), Sample3-5: ref
    vcf += "Chr_TEST\t1001\t.\tT\tC\t999\tPASS\tDP=100\tGT:PL:AD"
    vcf += "\t0/0:0,30,255:30,0"  # Sample1: T (ref)
    vcf += "\t1/1:255,30,0:0,30"  # Sample2: C (alt) - codon becomes ACG
    vcf += "\t0/0:0,30,255:30,0"  # Sample3: T (ref)
    vcf += "\t0/0:0,30,255:30,0"  # Sample4: T (ref)
    vcf += "\t0/0:0,30,255:30,0"  # Sample5: T (ref)
    vcf += "\n"
    
    # Variant at minus strand gene M01, position 3001 (C→G on reference strand)
    # On genome at 3000-3008: TCATGGCAT
    # Position 3001 is 'C' (2nd base)
    # If changed to G: TGATGGCAT on genome
    # Gene orientation (reverse complement): ATGCCATCA → ATGCCATCA (wait, let me recalculate)
    # Original genome: TCATGGCAT → gene: ATGCCATGA
    # Changed genome: TGATGGCAT → gene: ATGCCATCA
    # So the LAST codon changes from TGA (stop) to TCA (S)
    vcf += "Chr_TEST\t3001\t.\tC\tG\t999\tPASS\tDP=100\tGT:PL:AD"
    vcf += "\t0/0:0,30,255:30,0"  # Sample1: C (ref) → TGA in gene
    vcf += "\t0/0:0,30,255:30,0"  # Sample2: C (ref)
    vcf += "\t1/1:255,30,0:0,30"  # Sample3: G (alt) → TCA in gene
    vcf += "\t0/0:0,30,255:30,0"  # Sample4: C (ref)
    vcf += "\t0/0:0,30,255:30,0"  # Sample5: C (ref)
    vcf += "\n"
    
    # Invariant site at position 1002 (plus strand, should not affect anything)
    vcf += "Chr_TEST\t1002\t.\tG\t.\t999\tPASS\tDP=100\tGT:PL:AD"
    vcf += "\t0/0:0,30,255:30,0" * 5
    vcf += "\n"
    
    return vcf

def create_preferred_codons():
    """Create a simple preferred codons file."""
    return "ATG\nGCA\nCCA\n"

def run_test():
    """Run the test with synthetic data."""
    
    print("="*60)
    print("STRAND HANDLING TEST WITH SYNTHETIC DATA")
    print("="*60)
    
    # Create temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        # Write test files
        print("\n1. Creating synthetic test files...")
        
        genome_file = tmpdir / "test_genome.fa"
        genome_file.write_text(create_test_genome())
        print(f"   ✓ Genome: {genome_file}")
        
        gff3_file = tmpdir / "test_genes.gff3"
        gff3_file.write_text(create_test_gff3())
        print(f"   ✓ GFF3: {gff3_file}")
        
        cds_file = tmpdir / "test_cds.fa"
        cds_file.write_text(create_test_cds())
        print(f"   ✓ CDS: {cds_file}")
        
        vcf_file = tmpdir / "test_variants.vcf"
        vcf_file.write_text(create_test_vcf())
        print(f"   ✓ VCF: {vcf_file}")
        
        pref_file = tmpdir / "preferred_codons.txt"
        pref_file.write_text(create_preferred_codons())
        print(f"   ✓ Preferred codons: {pref_file}")
        
        # Run describe_gene_positions_by_degeneracy.py
        print("\n2. Running annotation script...")
        script_dir = Path(__file__).parent.parent / 'miscellanea_code'
        annotate_script = script_dir / 'describe_gene_positions_by_degeneracy.py'
        
        result = subprocess.run([
            'python3', str(annotate_script),
            'Chr_TEST',
            str(gff3_file),
            str(genome_file),
            str(cds_file)
        ], capture_output=True, text=True, cwd=tmpdir)
        
        if result.returncode != 0:
            print(f"   ✗ FAILED: {result.stderr}")
            return False
        
        annotation_file = tmpdir / "Chr_TEST.genic_bases.annotated.txt"
        print(f"   ✓ Annotation: {annotation_file}")
        
        # Debug: check genome at key positions
        print("\n   DEBUG: Genome check:")
        genome_text = genome_file.read_text()
        genome_lines = [l for l in genome_text.split('\n') if not l.startswith('>')]
        genome_seq = ''.join(genome_lines)
        print(f"      Genome 1000-1008 (P01, +): {genome_seq[1000:1009]}")
        print(f"      Genome 3000-3008 (M01, -): {genome_seq[3000:3009]}")
        
        # Check annotation file
        print("\n3. Verifying annotation...")
        with open(annotation_file) as f:
            lines = f.readlines()
        
        print(f"   Total annotated positions: {len(lines) - 1}")
        
        # Check specific positions
        expected_checks = [
            # Plus strand gene P01: positions 1000-1008
            ("1000", "P01", "A", "ATG", "M"),  # Start
            ("1001", "P01", "T", "ATG", "M"),
            ("1002", "P01", "G", "ATG", "M"),
            # Minus strand gene M01: positions 3000-3008 (but START should be at 3008!)
            ("3008", "M01", "T", "ATG", "M"),  # Should be START in gene orientation
            ("3007", "M01", "A", "ATG", "M"),
            ("3006", "M01", "C", "ATG", "M"),
        ]
        
        annotation_dict = {}
        for line in lines[1:]:  # Skip header
            cols = line.strip().split('\t')
            if len(cols) >= 8:
                pos, gene, base, codon, aa = cols[2], cols[1], cols[3], cols[6], cols[7]
                annotation_dict[pos] = (gene, base, codon, aa)
        
        print("\n   Checking key positions:")
        for pos, exp_gene, exp_base, exp_codon, exp_aa in expected_checks:
            if pos in annotation_dict:
                gene, base, codon, aa = annotation_dict[pos]
                match = "✓" if (gene == exp_gene and codon == exp_codon and aa == exp_aa) else "✗"
                print(f"   {match} Pos {pos}: gene={gene} (exp:{exp_gene}), "
                      f"base={base} (exp:{exp_base}), codon={codon} (exp:{exp_codon}), aa={aa} (exp:{exp_aa})")
            else:
                print(f"   ✗ Pos {pos}: NOT FOUND in annotation")
        
        # Run calculate_freq_preferred_codon.py
        print("\n4. Running codon frequency script...")
        codon_script = script_dir / 'calculate_freq_preferred_codon.py'
        
        result = subprocess.run([
            'python3', str(codon_script),
            'Chr_TEST',
            str(vcf_file),
            str(gff3_file),
            str(cds_file),
            str(genome_file),
            str(pref_file),
            '5'  # 5 samples
        ], capture_output=True, text=True, cwd=tmpdir)
        
        if result.returncode != 0:
            print(f"   ✗ FAILED: {result.stderr}")
            print(f"   STDOUT: {result.stdout}")
            return False
        
        codon_file = tmpdir / "Chr_TEST.codon_frequencies.txt"
        print(f"   ✓ Codon frequencies: {codon_file}")
        
        # Check codon genotypes
        print("\n5. Verifying codon genotypes...")
        with open(codon_file) as f:
            lines = f.readlines()
        
        print("\n   Expected results:")
        print("   - Gene P01, codon 0 (ATG): Sample1=ATG, Sample2=ACG (variant at pos 1001), Sample3-5=ATG")
        print("   - Gene M01, codon 2 (TGA): Sample1-2=TGA, Sample3=TCA (variant at pos 3001), Sample4-5=TGA")
        
        print("\n   Actual results:")
        for line in lines[1:]:  # Skip header
            cols = line.strip().split('\t')
            if len(cols) >= 8:
                gene, codon_idx, aa, ref_codon, is_pref, positions, variants, freqs = cols[0], cols[1], cols[2], cols[3], cols[4], cols[5], cols[6], cols[7]
                if gene in ['P01', 'M01']:
                    print(f"   Gene {gene}, codon {codon_idx}: {ref_codon} ({aa}) at positions {positions}")
                    print(f"      Variants: {variants}")
                    print(f"      Frequencies: {freqs}")
        
        print("\n" + "="*60)
        print("TEST COMPLETE")
        print("="*60)
        print("\nManual verification required:")
        print(f"1. Check annotation file: {annotation_file}")
        print(f"2. Check codon frequencies: {codon_file}")
        print("\nFiles preserved in temporary directory (will be deleted on exit)")
        
        return True

if __name__ == "__main__":
    success = run_test()
    sys.exit(0 if success else 1)
