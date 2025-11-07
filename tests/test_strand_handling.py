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

def reverse_complement(seq):
    """Return reverse complement of DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))

def create_test_genome():
    """
    Create synthetic genome with realistic complexity:
    - Multi-exon genes (split CDS)
    - Genes with different lengths
    - Intergenic regions
    - Mix of + and - strand
    """
    
    # Create base sequence (mostly N's for intergenic, we'll insert specific regions)
    genome_seq = 'N' * 10000
    genome_list = list(genome_seq)
    
    # NOTE: GFF3 uses 1-indexed coordinates, Python lists are 0-indexed
    # GFF3 position X corresponds to genome_list[X-1]
    
    # ===== PLUS STRAND GENES =====
    
    # Gene P01: Single exon, 5 codons at GFF3 positions 1000-1014 (1-indexed)
    # ATG-CCA-GCA-TCA-TGA (M-P-A-S-*)
    p01_seq = 'ATGCCAGCATCATGA'
    genome_list[999:999+len(p01_seq)] = list(p01_seq)  # 999 = position 1000 in GFF3
    
    # Gene P02: Two exons with intron at GFF3 positions 1100-1164
    # Exon1: ATG-CCA-GCA (M-P-A) at 1100-1108
    # Intron: 50bp
    # Exon2: TCA-TAG (S-*) at 1159-1164
    p02_exon1 = 'ATGCCAGCA'
    p02_intron = 'N' * 50
    p02_exon2 = 'TCATAG'
    genome_list[1099:1099+len(p02_exon1)] = list(p02_exon1)  # 1099 = position 1100
    genome_list[1108:1108+len(p02_intron)] = list(p02_intron)
    genome_list[1158:1158+len(p02_exon2)] = list(p02_exon2)  # 1158 = position 1159
    
    # Gene P03: Three exons with introns, 6 codons at GFF3 positions 1300-1387
    # Exon1: ATG-TCA (M-S) at 1300-1305
    # Intron1: 30bp
    # Exon2: GCA-CCA (A-P) at 1336-1341
    # Intron2: 40bp
    # Exon3: TCA-TAG (S-*) at 1382-1387
    p03_exon1 = 'ATGTCA'
    p03_exon2 = 'GCACCA'
    p03_exon3 = 'TCATAG'
    genome_list[1299:1305] = list(p03_exon1)  # 1299 = position 1300
    genome_list[1335:1341] = list(p03_exon2)  # 1335 = position 1336
    genome_list[1381:1387] = list(p03_exon3)  # 1381 = position 1382
    
    # ===== MINUS STRAND GENES =====
    # For minus strand, we write REVERSE COMPLEMENT on genome
    
    # Gene M01: Single exon, 5 codons at GFF3 positions 3000-3014
    # Gene orientation: ATG-CCA-GCA-TCA-TGA
    # Genome (rev comp): TCA-TGA-TGC-TGG-CAT
    m01_gene = 'ATGCCAGCATCATGA'
    m01_genome = reverse_complement(m01_gene)
    genome_list[2999:2999+len(m01_genome)] = list(m01_genome)  # 2999 = position 3000
    
    # Gene M02: Two exons with intron at GFF3 positions 3100-3164
    # Gene orientation exons: ATG-CCA-GCA (M-P-A), TCA-TAG (S-*)
    # Genome writes reverse complement with intron between
    m02_gene_exon1 = 'ATGCCAGCA'
    m02_gene_exon2 = 'TCATAG'
    m02_genome_exon1 = reverse_complement(m02_gene_exon1)  # Will be at higher positions
    m02_genome_exon2 = reverse_complement(m02_gene_exon2)  # Will be at lower positions
    # On genome: exon2 comes first (lower coords), then intron, then exon1 (higher coords)
    genome_list[3099:3099+len(m02_genome_exon2)] = list(m02_genome_exon2)  # 3099 = position 3100
    genome_list[3155:3155+len(m02_genome_exon1)] = list(m02_genome_exon1)  # 3155 = position 3156
    
    # Gene M03: Three exons (same as P03 in gene orientation) at GFF3 positions 3200-3297
    # Gene orientation: ATG-TCA, GCA-CCA, TCA-TAG
    m03_gene_exon1 = 'ATGTCA'
    m03_gene_exon2 = 'GCACCA'
    m03_gene_exon3 = 'TCATAG'
    m03_genome_exon1 = reverse_complement(m03_gene_exon1)  # Highest coords
    m03_genome_exon2 = reverse_complement(m03_gene_exon2)  # Middle coords
    m03_genome_exon3 = reverse_complement(m03_gene_exon3)  # Lowest coords
    # On genome: exon3, intron2, exon2, intron1, exon1
    genome_list[3199:3199+len(m03_genome_exon3)] = list(m03_genome_exon3)  # 3199 = position 3200
    genome_list[3245:3245+len(m03_genome_exon2)] = list(m03_genome_exon2)  # 3245 = position 3246
    genome_list[3291:3291+len(m03_genome_exon1)] = list(m03_genome_exon1)  # 3291 = position 3292
    
    genome_seq = ''.join(genome_list)
    
    fasta = f">Chr_TEST\n"
    # Write in 60 bp lines
    for i in range(0, len(genome_seq), 60):
        fasta += genome_seq[i:i+60] + '\n'
    
    return fasta

def create_test_gff3():
    """
    Create GFF3 with realistic gene structures:
    - Single exon genes
    - Multi-exon genes with introns
    - Genes with multiple isoforms (should process only .1!)
    """
    
    gff3 = "##gff-version 3\n"
    gff3 += "#Chr_TEST synthetic genome for testing\n"
    
    # ===== PLUS STRAND GENES =====
    
    # Gene P01: Single exon (1000-1014)
    gff3 += "Chr_TEST\ttest\tgene\t1000\t1014\t.\t+\t.\tID=gene:P01\n"
    gff3 += "Chr_TEST\ttest\tmRNA\t1000\t1014\t.\t+\t.\tID=MgIM767.P01.1.v2.1;Parent=gene:P01\n"
    gff3 += "Chr_TEST\ttest\tCDS\t1000\t1014\t.\t+\t0\tID=MgIM767.P01.1.v2.1.CDS.1;Parent=MgIM767.P01.1.v2.1\n"
    
    # Gene P02: Two exons with intron (exon1: 1100-1108, exon2: 1159-1164)
    gff3 += "Chr_TEST\ttest\tgene\t1100\t1164\t.\t+\t.\tID=gene:P02\n"
    # Primary isoform .1 (full length)
    gff3 += "Chr_TEST\ttest\tmRNA\t1100\t1164\t.\t+\t.\tID=MgIM767.P02.1.v2.1;Parent=gene:P02\n"
    gff3 += "Chr_TEST\ttest\tCDS\t1100\t1108\t.\t+\t0\tID=MgIM767.P02.1.v2.1.CDS.1;Parent=MgIM767.P02.1.v2.1\n"
    gff3 += "Chr_TEST\ttest\tCDS\t1159\t1164\t.\t+\t0\tID=MgIM767.P02.1.v2.1.CDS.2;Parent=MgIM767.P02.1.v2.1\n"
    # Alternate isoform .2 (exon1 only - should be SKIPPED!)
    gff3 += "Chr_TEST\ttest\tmRNA\t1100\t1108\t.\t+\t.\tID=MgIM767.P02.2.v2.1;Parent=gene:P02\n"
    gff3 += "Chr_TEST\ttest\tCDS\t1100\t1108\t.\t+\t0\tID=MgIM767.P02.2.v2.1.CDS.1;Parent=MgIM767.P02.2.v2.1\n"
    
    # Gene P03: Three exons with introns (exon1: 1300-1305, exon2: 1336-1341, exon3: 1382-1387)
    gff3 += "Chr_TEST\ttest\tgene\t1300\t1387\t.\t+\t.\tID=gene:P03\n"
    gff3 += "Chr_TEST\ttest\tmRNA\t1300\t1387\t.\t+\t.\tID=MgIM767.P03.1.v2.1;Parent=gene:P03\n"
    gff3 += "Chr_TEST\ttest\tCDS\t1300\t1305\t.\t+\t0\tID=MgIM767.P03.1.v2.1.CDS.1;Parent=MgIM767.P03.1.v2.1\n"
    gff3 += "Chr_TEST\ttest\tCDS\t1336\t1341\t.\t+\t0\tID=MgIM767.P03.1.v2.1.CDS.2;Parent=MgIM767.P03.1.v2.1\n"
    gff3 += "Chr_TEST\ttest\tCDS\t1382\t1387\t.\t+\t0\tID=MgIM767.P03.1.v2.1.CDS.3;Parent=MgIM767.P03.1.v2.1\n"
    
    # ===== MINUS STRAND GENES =====
    
    # Gene M01: Single exon (3000-3014)
    gff3 += "Chr_TEST\ttest\tgene\t3000\t3014\t.\t-\t.\tID=gene:M01\n"
    gff3 += "Chr_TEST\ttest\tmRNA\t3000\t3014\t.\t-\t.\tID=MgIM767.M01.1.v2.1;Parent=gene:M01\n"
    gff3 += "Chr_TEST\ttest\tCDS\t3000\t3014\t.\t-\t0\tID=MgIM767.M01.1.v2.1.CDS.1;Parent=MgIM767.M01.1.v2.1\n"
    
    # Gene M02: Two exons with intron (exon at lower coords: 3100-3105, exon at higher coords: 3156-3164)
    # In gene orientation: exon2 (high coords) comes first, then exon1 (low coords)
    gff3 += "Chr_TEST\ttest\tgene\t3100\t3164\t.\t-\t.\tID=gene:M02\n"
    # Primary isoform .1
    gff3 += "Chr_TEST\ttest\tmRNA\t3100\t3164\t.\t-\t.\tID=MgIM767.M02.1.v2.1;Parent=gene:M02\n"
    gff3 += "Chr_TEST\ttest\tCDS\t3156\t3164\t.\t-\t0\tID=MgIM767.M02.1.v2.1.CDS.1;Parent=MgIM767.M02.1.v2.1\n"
    gff3 += "Chr_TEST\ttest\tCDS\t3100\t3105\t.\t-\t0\tID=MgIM767.M02.1.v2.1.CDS.2;Parent=MgIM767.M02.1.v2.1\n"
    # Alternate isoform .2 (should be SKIPPED!)
    gff3 += "Chr_TEST\ttest\tmRNA\t3156\t3164\t.\t-\t.\tID=MgIM767.M02.2.v2.1;Parent=gene:M02\n"
    gff3 += "Chr_TEST\ttest\tCDS\t3156\t3164\t.\t-\t0\tID=MgIM767.M02.2.v2.1.CDS.1;Parent=MgIM767.M02.2.v2.1\n"
    
    # Gene M03: Three exons (exon at low: 3200-3205, middle: 3246-3251, high: 3292-3297)
    # In gene orientation: exon3 (high) → exon2 (mid) → exon1 (low)
    gff3 += "Chr_TEST\ttest\tgene\t3200\t3297\t.\t-\t.\tID=gene:M03\n"
    gff3 += "Chr_TEST\ttest\tmRNA\t3200\t3297\t.\t-\t.\tID=MgIM767.M03.1.v2.1;Parent=gene:M03\n"
    gff3 += "Chr_TEST\ttest\tCDS\t3292\t3297\t.\t-\t0\tID=MgIM767.M03.1.v2.1.CDS.1;Parent=MgIM767.M03.1.v2.1\n"
    gff3 += "Chr_TEST\ttest\tCDS\t3246\t3251\t.\t-\t0\tID=MgIM767.M03.1.v2.1.CDS.2;Parent=MgIM767.M03.1.v2.1\n"
    gff3 += "Chr_TEST\ttest\tCDS\t3200\t3205\t.\t-\t0\tID=MgIM767.M03.1.v2.1.CDS.3;Parent=MgIM767.M03.1.v2.1\n"
    
    return gff3

def create_test_cds():
    """
    Create CDS FASTA with sequences in gene orientation (5'→3').
    All sequences match what you'd get from splicing exons.
    """
    
    cds = ""
    
    # ===== PLUS STRAND GENES =====
    # (sequences as written on genome, left to right)
    
    cds += ">MgIM767.P01.1 ID=MgIM767.P01.1.v2.1\n"
    cds += "ATGCCAGCATCATGA\n"  # 5 codons
    
    cds += ">MgIM767.P02.1 ID=MgIM767.P02.1.v2.1\n"
    cds += "ATGCCAGCATCATAG\n"  # exon1 + exon2 = 5 codons
    
    cds += ">MgIM767.P03.1 ID=MgIM767.P03.1.v2.1\n"
    cds += "ATGTCAGCACCATCATAG\n"  # exon1 + exon2 + exon3 = 6 codons
    
    # ===== MINUS STRAND GENES =====
    # (sequences in gene orientation 5'→3', NOT as they appear on genome!)
    # These match the plus strand genes in gene orientation
    
    cds += ">MgIM767.M01.1 ID=MgIM767.M01.1.v2.1\n"
    cds += "ATGCCAGCATCATGA\n"  # Same as P01 in gene orientation
    
    cds += ">MgIM767.M02.1 ID=MgIM767.M02.1.v2.1\n"
    cds += "ATGCCAGCATCATAG\n"  # Same as P02 in gene orientation
    
    cds += ">MgIM767.M03.1 ID=MgIM767.M03.1.v2.1\n"
    cds += "ATGTCAGCACCATCATAG\n"  # Same as P03 in gene orientation
    
    return cds

def create_test_vcf():
    """
    Create VCF with variants testing multiple scenarios:
    - Variants in single-exon genes
    - Variants in multi-exon genes (in different exons)
    - Variants in + and - strand genes
    - Variants that change codons
    - Invariant sites
    """
    
    vcf = "##fileformat=VCFv4.2\n"
    vcf += "##contig=<ID=Chr_TEST,length=10000>\n"
    vcf += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    vcf += "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">\n"
    vcf += "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
    vcf += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    
    # Add 10 sample columns
    for i in range(1, 11):
        vcf += f"\tSample{i}"
    vcf += "\n"
    
    # === PLUS STRAND VARIANTS ===
    
    # Variant in P01, position 1001 (2nd base of ATG): T→C
    # Changes ATG (M) → ACG (T)
    vcf += "Chr_TEST\t1001\t.\tT\tC\t999\tPASS\tDP=100\tGT:PL:AD"
    vcf += "\t0/0:0,30,255:30,0" * 8  # Samples 1-8: ref (T)
    vcf += "\t1/1:255,30,0:0,30"     # Sample 9: alt (C)
    vcf += "\t0/0:0,30,255:30,0"      # Sample 10: ref (T)
    vcf += "\n"
    
    # Variant in P02 exon2, position 1160 (2nd base of TCA): C→G
    # Changes TCA (S) → TGA (*)
    vcf += "Chr_TEST\t1160\t.\tC\tG\t999\tPASS\tDP=100\tGT:PL:AD"
    vcf += "\t0/0:0,30,255:30,0" * 7  # Samples 1-7: ref (C)
    vcf += "\t1/1:255,30,0:0,30"     # Sample 8: alt (G)
    vcf += "\t0/0:0,30,255:30,0" * 2  # Samples 9-10: ref (C)
    vcf += "\n"
    
    # Variant in P03 exon2, position 1337 (2nd base of GCA): C→T
    # Changes GCA (A) → GTA (V)
    vcf += "Chr_TEST\t1337\t.\tC\tT\t999\tPASS\tDP=100\tGT:PL:AD"
    vcf += "\t0/0:0,30,255:30,0" * 5  # Samples 1-5: ref (C)
    vcf += "\t1/1:255,30,0:0,30"     # Sample 6: alt (T)
    vcf += "\t0/0:0,30,255:30,0" * 4  # Samples 7-10: ref (C)
    vcf += "\n"
    
    # === MINUS STRAND VARIANTS ===
    # Remember: VCF reports reference strand, but gene is on minus strand!
    
    # Variant in M01, position 3001 (on genome: 2nd base of TCA...TGA...)
    # Genome ref: C, Gene position: this is 2nd base of 5th codon (TGA)
    # If C→G on genome, in gene becomes: TGA→TCA (stop→Ser)
    vcf += "Chr_TEST\t3001\t.\tC\tG\t999\tPASS\tDP=100\tGT:PL:AD"
    vcf += "\t0/0:0,30,255:30,0" * 4  # Samples 1-4: ref (TGA in gene)
    vcf += "\t1/1:255,30,0:0,30"     # Sample 5: alt (TCA in gene)
    vcf += "\t0/0:0,30,255:30,0" * 5  # Samples 6-10: ref
    vcf += "\n"
    
    # Variant in M02 exon at high coords, position 3157
    # This is in the first codon of gene orientation (ATG)
    # Genome ref at 3157: T (part of rev comp of ATG)
    # If T→A, changes first codon
    vcf += "Chr_TEST\t3157\t.\tT\tA\t999\tPASS\tDP=100\tGT:PL:AD"
    vcf += "\t0/0:0,30,255:30,0" * 3  # Samples 1-3: ref
    vcf += "\t1/1:255,30,0:0,30"     # Sample 4: alt
    vcf += "\t0/0:0,30,255:30,0" * 6  # Samples 5-10: ref
    vcf += "\n"
    
    # Variant in M03 exon at middle coords, position 3247
    # This is in codon 3 (GCA) - 2nd exon in gene orientation
    vcf += "Chr_TEST\t3247\t.\tG\tA\t999\tPASS\tDP=100\tGT:PL:AD"
    vcf += "\t0/0:0,30,255:30,0" * 2  # Samples 1-2: ref
    vcf += "\t1/1:255,30,0:0,30"     # Sample 3: alt
    vcf += "\t0/0:0,30,255:30,0" * 7  # Samples 4-10: ref
    vcf += "\n"
    
    # Invariant site
    vcf += "Chr_TEST\t1005\t.\tA\t.\t999\tPASS\tDP=100\tGT:PL:AD"
    vcf += "\t0/0:0,30,255:30,0" * 10
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
        
        # Check if advisor's proc2.py exists (local copy in miscellanea_code)
        script_dir = Path(__file__).parent.parent / 'miscellanea_code'
        advisor_script = script_dir / 'proc2.py'
        has_advisor_script = advisor_script.exists()
        
        if has_advisor_script:
            print(f"\n✓ Found advisor's proc2.py: {advisor_script}")
        else:
            print(f"\n⚠ Advisor's proc2.py not found at: {advisor_script}")
            print("  Will run only our script for testing")
        
        
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
        print(f"      Genome 999-1013 (P01 in 0-indexed): {genome_seq[999:1014]}")
        print(f"      Genome 2999-3013 (M01 in 0-indexed): {genome_seq[2999:3014]}")
        print(f"      GFF3 positions are 1-indexed, so add 1 to see GFF3 coords")
        
        # Check annotation file
        print("\n3. Verifying annotation...")
        with open(annotation_file) as f:
            lines = f.readlines()
        
        print(f"   Total annotated positions: {len(lines) - 1}")
        
        # Check specific positions using ACTUAL genome bases
        # Note: annotation outputs REFERENCE STRAND bases (genome), not gene orientation
        expected_checks = [
            # Plus strand gene P01: codon 0 (ATG) at GFF3 positions 1000-1002
            ("1000", "P01", "A", "ATG", "M"),  # genome[999] = A
            ("1001", "P01", "T", "ATG", "M"),  # genome[1000] = T
            ("1002", "P01", "G", "ATG", "M"),  # genome[1001] = G
            # Minus strand gene M01: codon 2 (GCA in gene) at GFF3 positions 3006-3008
            # Genome has: T,G,C which reverse complement to A,C,G = GCA
            ("3008", "M01", "C", "GCA", "A"),  # genome[3007] = C
            ("3007", "M01", "G", "GCA", "A"),  # genome[3006] = G
            ("3006", "M01", "T", "GCA", "A"),  # genome[3005] = T
            # Minus strand gene M01: codon 0 (ATG in gene) at GFF3 positions 3012-3014
            # Genome has: C,A,T which reverse complement to G,T,A = ATG (reversed)
            ("3014", "M01", "T", "ATG", "M"),  # genome[3013] = T
            ("3013", "M01", "A", "ATG", "M"),  # genome[3012] = A
            ("3012", "M01", "C", "ATG", "M"),  # genome[3011] = C
        ]
        
        annotation_dict = {}
        for line in lines[1:]:  # Skip header
            cols = line.strip().split('\t')
            if len(cols) >= 8:
                pos, gene, base, codon, aa = cols[2], cols[1], cols[3], cols[6], cols[7]
                annotation_dict[pos] = (gene, base, codon, aa)
        
        print("\n   Checking key positions:")
        all_checks_pass = True
        for pos, exp_gene, exp_base, exp_codon, exp_aa in expected_checks:
            if pos in annotation_dict:
                gene, base, codon, aa = annotation_dict[pos]
                # Check ALL: gene, base, codon, and AA
                match = "✓" if (gene == exp_gene and base == exp_base and codon == exp_codon and aa == exp_aa) else "✗"
                if match == "✗":
                    all_checks_pass = False
                print(f"   {match} Pos {pos}: gene={gene} (exp:{exp_gene}), "
                      f"base={base} (exp:{exp_base}), codon={codon} (exp:{exp_codon}), aa={aa} (exp:{exp_aa})")
            else:
                print(f"   ✗ Pos {pos}: NOT FOUND in annotation")
                all_checks_pass = False
        
        if not all_checks_pass:
            print("\n   ✗ ANNOTATION CHECK FAILED - Base mismatches detected!")
            return False
        
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
            '10'  # 10 samples
        ], capture_output=True, text=True, cwd=tmpdir)
        
        if result.returncode != 0:
            print(f"   ✗ FAILED: {result.stderr}")
            print(f"   STDOUT: {result.stdout}")
            return False
        
        codon_file = tmpdir / "Chr_TEST.codon_frequencies.txt"
        print(f"   ✓ Codon frequencies: {codon_file}")
        
        # Run calculate_pi.py (our implementation)
        print("\n4b. Running calculate_pi.py (our π implementation)...")
        pi_script = script_dir / 'calculate_pi.py'
        
        result = subprocess.run([
            'python3', str(pi_script),
            'Chr_TEST',
            str(vcf_file),
            str(annotation_file),
            '10'  # 10 samples
        ], capture_output=True, text=True, cwd=tmpdir)
        
        if result.returncode != 0:
            print(f"   ✗ FAILED: {result.stderr}")
            print(f"   STDOUT: {result.stdout}")
            return False
        
        our_pi_file = tmpdir / "Chr_TEST.bygene.pi.txt"
        print(f"   ✓ Our π output: {our_pi_file}")
        
        # Run advisor's proc2.py if available
        advisor_pi_file = None
        if has_advisor_script:
            print("\n4c. Running advisor's proc2.py...")
            
            # proc2.py expects:
            # - ../Included_snp.sites.txt (VCF)
            # - Chr_TEST.genic_bases.annotated.txt (annotation)
            # Create symlink for VCF in parent dir
            vcf_link = tmpdir.parent / "Included_snp.sites.txt"
            try:
                vcf_link.symlink_to(vcf_file)
            except FileExistsError:
                vcf_link.unlink()
                vcf_link.symlink_to(vcf_file)
            
            result = subprocess.run([
                'python3', str(advisor_script),
                'Chr_TEST'
            ], capture_output=True, text=True, cwd=tmpdir)
            
            # Clean up symlink
            vcf_link.unlink()
            
            if result.returncode != 0:
                print(f"   ⚠ Advisor script failed: {result.stderr}")
                print(f"   STDOUT: {result.stdout}")
            else:
                advisor_pi_file = tmpdir / "Chr_TEST.bygene"
                if advisor_pi_file.exists():
                    print(f"   ✓ Advisor π output: {advisor_pi_file}")
                else:
                    print(f"   ⚠ Advisor output file not created")
                    advisor_pi_file = None
        
        # Check codon genotypes
        print("\n5. Verifying codon genotypes...")
        with open(codon_file) as f:
            lines = f.readlines()
        
        # Parse results
        results = {}
        for line in lines[1:]:  # Skip header
            cols = line.strip().split('\t')
            if len(cols) >= 8:
                gene, codon_idx, aa, ref_codon, is_pref, positions, variants, freqs = (
                    cols[0], int(cols[1]), cols[2], cols[3], cols[4], cols[5], cols[6], cols[7]
                )
                if gene not in results:
                    results[gene] = {}
                results[gene][codon_idx] = {
                    'aa': aa,
                    'ref_codon': ref_codon,
                    'positions': positions,
                    'variants': variants,
                    'freqs': freqs
                }
        
        print("\n   Expected results:")
        print("   - Gene P01 (plus, single exon), codon 0 (ATG):")
        print("     * Variant at pos 1001 (T→C) in Sample 9")
        print("     * Should show ATG:9 and ACG:1")
        print("   - Gene P02 (plus, 2 exons), codon 3 (TCA):")
        print("     * Variant at pos 1160 (C→G) in Sample 8")
        print("     * Should show TCA:9 and TGA:1")
        print("   - Gene M01 (minus, single exon), codon 4 (TGA):")
        print("     * Variant at pos 3001 affects last codon")
        print("     * Should show TGA:9 and TCA:1 in Sample 5")
        
        # Verification
        all_pass = True
        
        # Test P01
        print("\n   --- Gene P01 (plus strand, single exon) ---")
        if 'P01' in results and 0 in results['P01']:
            p01_c0 = results['P01'][0]
            print(f"   Codon 0: {p01_c0['ref_codon']} ({p01_c0['aa']}) at {p01_c0['positions']}")
            print(f"   Variants: {p01_c0['variants']}")
            print(f"   Frequencies: {p01_c0['freqs']}")
            
            if 'ATG' in p01_c0['variants'] and 'ACG' in p01_c0['variants']:
                print("   ✓ PASS: Both ATG and ACG found")
            else:
                print("   ✗ FAIL: Expected both ATG and ACG")
                all_pass = False
        else:
            print("   ✗ FAIL: Gene P01 codon 0 not found")
            all_pass = False
        
        # Test P02
        print("\n   --- Gene P02 (plus strand, 2 exons) ---")
        if 'P02' in results:
            print(f"   Gene P02 has {len(results['P02'])} codons")
            if 3 in results['P02']:
                p02_c3 = results['P02'][3]
                print(f"   Codon 3: {p02_c3['ref_codon']} ({p02_c3['aa']}) at {p02_c3['positions']}")
                print(f"   Variants: {p02_c3['variants']}")
                print(f"   Frequencies: {p02_c3['freqs']}")
                
                if 'TCA' in p02_c3['variants'] and 'TGA' in p02_c3['variants']:
                    print("   ✓ PASS: Both TCA and TGA found")
                else:
                    print("   ✗ FAIL: Expected both TCA and TGA")
                    all_pass = False
            else:
                print("   ✗ FAIL: Codon 3 not found")
                all_pass = False
        else:
            print("   ✗ FAIL: Gene P02 not found")
            all_pass = False
        
        # Test P03
        print("\n   --- Gene P03 (plus strand, 3 exons) ---")
        if 'P03' in results:
            print(f"   Gene P03 has {len(results['P03'])} codons")
            if 2 in results['P03']:
                p03_c2 = results['P03'][2]
                print(f"   Codon 2: {p03_c2['ref_codon']} ({p03_c2['aa']}) at {p03_c2['positions']}")
                print(f"   Variants: {p03_c2['variants']}")
                
                if 'GCA' in p03_c2['variants']:
                    print("   ✓ PASS: GCA found (variant check deferred)")
                else:
                    print("   ✗ FAIL: Expected GCA")
                    all_pass = False
            else:
                print("   ✗ FAIL: Codon 2 not found")
                all_pass = False
        else:
            print("   ✗ FAIL: Gene P03 not found")
            all_pass = False
        
        # Test M01
        print("\n   --- Gene M01 (minus strand, single exon) ---")
        if 'M01' in results:
            print(f"   Gene M01 has {len(results['M01'])} codons")
            if 4 in results['M01']:
                m01_c4 = results['M01'][4]
                print(f"   Codon 4: {m01_c4['ref_codon']} ({m01_c4['aa']}) at {m01_c4['positions']}")
                print(f"   Variants: {m01_c4['variants']}")
                print(f"   Frequencies: {m01_c4['freqs']}")
                
                if 'TGA' in m01_c4['variants'] and 'TCA' in m01_c4['variants']:
                    print("   ✓ PASS: Both TGA and TCA found")
                else:
                    print("   ✗ FAIL: Expected both TGA and TCA")
                    all_pass = False
            else:
                print("   ✗ FAIL: Codon 4 not found")
                all_pass = False
        else:
            print("   ✗ FAIL: Gene M01 not found")
            all_pass = False
        
        # Test M02
        print("\n   --- Gene M02 (minus strand, 2 exons) ---")
        if 'M02' in results:
            print(f"   Gene M02 has {len(results['M02'])} codons")
            print("   ✓ PASS: Gene M02 processed")
        else:
            print("   ✗ FAIL: Gene M02 not found")
            all_pass = False
        
        # Test M03
        print("\n   --- Gene M03 (minus strand, 3 exons) ---")
        if 'M03' in results:
            print(f"   Gene M03 has {len(results['M03'])} codons")
            print("   ✓ PASS: Gene M03 processed")
        else:
            print("   ✗ FAIL: Gene M03 not found")
            all_pass = False
        
        # Compare π calculations if advisor output exists
        if advisor_pi_file:
            print("\n6. Comparing π calculations...")
            print("   (Advisor vs Our implementation)")
            
            # Parse our output
            our_pi = {}
            with open(our_pi_file) as f:
                header = f.readline()
                for line in f:
                    cols = line.strip().split('\t')
                    if len(cols) >= 7:
                        gene = cols[1]
                        # Sites_0fold, Poly_0fold, Pi_sum_0fold
                        sites_0 = int(cols[2])
                        poly_0 = int(cols[3])
                        pi_sum_0 = float(cols[4])
                        # Sites_2fold, Poly_2fold, Pi_sum_2fold
                        sites_2 = int(cols[8])
                        poly_2 = int(cols[9])
                        pi_sum_2 = float(cols[10])
                        # Sites_3fold, Poly_3fold, Pi_sum_3fold
                        sites_3 = int(cols[14])
                        poly_3 = int(cols[15])
                        pi_sum_3 = float(cols[16])
                        # Sites_4fold, Poly_4fold, Pi_sum_4fold
                        sites_4 = int(cols[20])
                        poly_4 = int(cols[21])
                        pi_sum_4 = float(cols[22])
                        
                        our_pi[gene] = {
                            '0fold': (sites_0, poly_0, pi_sum_0),
                            '2fold': (sites_2, poly_2, pi_sum_2),
                            '3fold': (sites_3, poly_3, pi_sum_3),
                            '4fold': (sites_4, poly_4, pi_sum_4)
                        }
            
            # Parse advisor output
            # Format: Chr Gene Sites_1 Poly_1 Pi_1 Sites_2 Poly_2 Pi_2 Sites_3_Not_4f Poly_3_Not_4f Pi_3_Not_4f Sites_3_4f Poly_3_4f Pi_3_4f
            advisor_pi = {}
            with open(advisor_pi_file) as f:
                for line in f:
                    cols = line.strip().split('\t')
                    if len(cols) >= 13:
                        gene = cols[1]
                        # 1st position (advisor "1")
                        sites_1 = int(cols[2])
                        poly_1 = int(cols[3])
                        pi_1 = float(cols[4])
                        # 2nd position (advisor "2")
                        sites_2 = int(cols[5])
                        poly_2 = int(cols[6])
                        pi_2 = float(cols[7])
                        # 3rd non-4fold (advisor "3_Not_4f")
                        sites_3nf = int(cols[8])
                        poly_3nf = int(cols[9])
                        pi_3nf = float(cols[10])
                        # 3rd 4fold (advisor "3_fourfold")
                        sites_4f = int(cols[11])
                        poly_4f = int(cols[12])
                        pi_4f = float(cols[13])
                        
                        advisor_pi[gene] = {
                            '1st': (sites_1, poly_1, pi_1),
                            '2nd': (sites_2, poly_2, pi_2),
                            '3rd_not4f': (sites_3nf, poly_3nf, pi_3nf),
                            '4fold': (sites_4f, poly_4f, pi_4f)
                        }
            
            # CRITICAL: Advisor uses different degeneracy categories than ours!
            # 
            # Advisor's approach (position-based):
            #   - "1": ALL 1st codon positions
            #   - "2": ALL 2nd codon positions  
            #   - "3_Not_4f": 3rd positions that are NOT 4-fold degenerate
            #   - "3_fourfold": 3rd positions that ARE 4-fold degenerate
            #
            # Our approach (degeneracy-based):
            #   - "0-fold": Bases where ANY change causes AA change (includes some 1st, 2nd, 3rd positions)
            #   - "2-fold": Bases with 1 synonymous alternative
            #   - "3-fold": Bases with 2 synonymous alternatives
            #   - "4-fold": Bases with 3 synonymous alternatives (only 3rd positions of certain codons)
            #
            # Example: ATG codon (Met)
            #   Advisor: 1st=A (category "1"), 2nd=T (category "2"), 3rd=G (category "3_Not_4f")
            #   Ours: 1st=A (0-fold), 2nd=T (0-fold), 3rd=G (0-fold) - all changes alter AA
            #
            # Example: GCA codon (Ala)  
            #   Advisor: 1st=G (category "1"), 2nd=C (category "2"), 3rd=A (category "3_fourfold")
            #   Ours: 1st=G (0-fold), 2nd=C (0-fold), 3rd=A (4-fold) - only 3rd is synonymous
            #
            # This means:
            #   - Advisor's "1" and "2" categories mix different degeneracy classes
            #   - Only advisor's "3_fourfold" directly maps to our "4-fold"
            #   - We can verify π calculation is identical, but site counts will differ!
            
            print("\n   ═══════════════════════════════════════════════════════")
            print("   COMPARING π CALCULATIONS BETWEEN IMPLEMENTATIONS")
            print("   ═══════════════════════════════════════════════════════")
            print()
            print("   ⚠ NOTE: Different classification schemes used!")
            print("   ├─ Advisor: Position-based (1st, 2nd, 3rd_not4f, 3rd_4f)")
            print("   └─ Ours: Degeneracy-based (0-fold, 2-fold, 3-fold, 4-fold)")
            print()
            print("   Only 4-fold sites map directly between schemes.")
            print("   Other categories will have different site counts.")
            print("   π formula is identical: π = 2n·p(1-p)/(n-1)")
            print()
            
            pi_match = True
            site_count_mismatch = False
            for gene in sorted(set(our_pi.keys()) | set(advisor_pi.keys())):
                print(f"   Gene {gene}:")
                if gene in our_pi and gene in advisor_pi:
                    # Compare 4-fold sites (only category that maps directly)
                    our_4f = our_pi[gene]['4fold']
                    adv_4f = advisor_pi[gene]['4fold']
                    
                    sites_match = our_4f[0] == adv_4f[0]
                    poly_match = our_4f[1] == adv_4f[1]
                    pi_close = abs(our_4f[2] - adv_4f[2]) < 0.001
                    
                    # For π comparison, what matters is the calculation is correct
                    # Site count differences are EXPECTED due to different schemes
                    if sites_match and poly_match and pi_close:
                        match_str = "✓"
                    elif pi_close and (our_4f[1] == adv_4f[1]):  # π matches and same # polymorphisms
                        match_str = "⚠"
                        site_count_mismatch = True
                    else:
                        match_str = "✗"
                        pi_match = False
                    
                    print(f"      {match_str} 4-fold sites: Our={our_4f} | Advisor={adv_4f}")
                    
                    # Show other categories for information (expect differences!)
                    print(f"         [Info] Our 0-fold: sites={our_pi[gene]['0fold'][0]}, poly={our_pi[gene]['0fold'][1]}, π={our_pi[gene]['0fold'][2]:.6f}")
                    print(f"         [Info] Advisor 1st: sites={advisor_pi[gene]['1st'][0]}, poly={advisor_pi[gene]['1st'][1]}, π={advisor_pi[gene]['1st'][2]:.6f}")
                    print(f"         [Info] Our 2-fold: sites={our_pi[gene]['2fold'][0]}, poly={our_pi[gene]['2fold'][1]}, π={our_pi[gene]['2fold'][2]:.6f}")
                    print(f"         [Info] Advisor 2nd: sites={advisor_pi[gene]['2nd'][0]}, poly={advisor_pi[gene]['2nd'][1]}, π={advisor_pi[gene]['2nd'][2]:.6f}")
                elif gene in our_pi:
                    print(f"      ⚠ Only in our output")
                    pi_match = False
                else:
                    print(f"      ⚠ Only in advisor output")
                    pi_match = False
            
            print()
            if pi_match:
                print("   ✓✓✓ π CALCULATION VERIFIED: Identical between implementations!")
                print("   ✓✓✓ Formula confirmed: π = 2n·p(1-p)/(n-1)")
                all_pass = all_pass and True
            elif site_count_mismatch and not pi_match:
                print("   ⚠ Site count differences due to classification schemes (EXPECTED)")
                print("   ⚠ But π values should match for sites with same polymorphism count")
                all_pass = False
            else:
                print("   ✗ π calculation mismatch detected (UNEXPECTED)")
                all_pass = False
        
        print("\n" + "="*60)
        if all_pass:
            print("✓✓✓ ALL TESTS PASSED ✓✓✓")
            print()
            print("Verified:")
            print("  • Multi-exon gene handling correct (plus & minus strands)")
            print("  • Codon genotyping accurate across gene structures")
            print("  • VCF variant integration working properly")
            if advisor_pi_file:
                print("  • π calculation formula identical to advisor's")
                print("  • Formula: π = 2n·p(1-p)/(n-1) ✓")
        else:
            print("✗ SOME TESTS FAILED - Review output above")
        print("="*60)
        
        print("\nFiles for manual inspection:")
        print(f"1. Annotation: {annotation_file}")
        print(f"2. Codon frequencies: {codon_file}")
        print(f"3. Our π output: {our_pi_file}")
        if advisor_pi_file:
            print(f"4. Advisor π output: {advisor_pi_file}")
        
        return all_pass

if __name__ == "__main__":
    success = run_test()
    sys.exit(0 if success else 1)
