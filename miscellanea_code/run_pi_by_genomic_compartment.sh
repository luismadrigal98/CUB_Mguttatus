#!/bin/bash
#SBATCH --job-name=pi_compartments
#SBATCH --output=logs/pi_compartments_%j.out
#SBATCH --error=logs/pi_compartments_%j.err
#SBATCH --time=10-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=kucg,eeb,kelly
#
# run_pi_by_genomic_compartment.sh
#
# Calculate nucleotide diversity (π) across genomic compartments:
#   1. Intergenic windows (50kb) - neutral baseline
#   2. Introns (50bp trimmed) - neutral within genes
#   3. First exons - 4-fold degenerate sites only
#   4. Non-first exons - 4-fold degenerate sites only
#
# Usage:
#   sbatch run_pi_by_genomic_compartment.sh
#   OR
#   bash run_pi_by_genomic_compartment.sh
#
# Author: Luis Javier Madrigal-Roca & GitHub Copilot
# Date: 2026-01-14

set -e
set -u

echo "========================================="
echo "π by Genomic Compartment Analysis"
echo "========================================="
echo "Started: $(date)"
echo ""

# ============================================================================
# Configuration
# ============================================================================

# Input files
GFF3="/home/l338m483/scratch/IMLines_to_767/data/Mguttatusvar_IM767_887_v2.1.gene.gff3"
GENOME_FA="/home/l338m483/scratch/IMLines_to_767/data/Mguttatusvar_IM767_887_v2.0.fa"
CDS_FA="/home/l338m483/scratch/IMLines_to_767/data/Mguttatusvar_IM767_887_v2.1.cds_primaryTranscriptOnly.fa"
VCF="/home/l338m483/scratch/IMLines_to_767/Included_snp.sites.txt"  # Your all-sites VCF

# Python scripts
SCRIPT_DIR="/home/l338m483/scratch/IMLines_to_767/miscellanea_code"
MISCELLANEA_DIR="/home/l338m483/scratch/IMLines_to_767/miscellanea_code"

# Output
OUTPUT_DIR="results/pi_compartment_analysis"
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# Chromosomes
CHROMOSOMES=(Chr_01 Chr_02 Chr_03 Chr_04 Chr_05 Chr_06 Chr_07 Chr_08 Chr_09 Chr_10 Chr_11 Chr_12 Chr_13 Chr_14)

# ============================================================================
# Step 1: Generate degeneracy annotations (if not already done)
# ============================================================================

echo "Step 1: Checking degeneracy annotation files..."

DEGENERACY_FILES=""
MISSING_ANNOT=0

for CHR in "${CHROMOSOMES[@]}"; do
    ANNOT_FILE="$OUTPUT_DIR/${CHR}.genic_bases.annotated.txt"
    
    if [ ! -f "$ANNOT_FILE" ]; then
        echo "  Generating annotation for $CHR..."
        
        python3 "${MISCELLANEA_DIR}/describe_gene_positions_by_degeneracy.py" \
            "$CHR" \
            "$GFF3" \
            "$GENOME_FA" \
            "$CDS_FA"
        
        # Move to output directory
        if [ -f "${CHR}.genic_bases.annotated.txt" ]; then
            mv "${CHR}.genic_bases.annotated.txt" "$ANNOT_FILE"
        fi
    else
        echo "  ✓ $CHR annotation exists"
    fi
    
    if [ -f "$ANNOT_FILE" ]; then
        DEGENERACY_FILES="$DEGENERACY_FILES $ANNOT_FILE"
    else
        echo "  ✗ WARNING: Missing annotation for $CHR"
        MISSING_ANNOT=1
    fi
done

if [ $MISSING_ANNOT -eq 1 ]; then
    echo ""
    echo "ERROR: Some degeneracy annotation files are missing."
    echo "Please run describe_gene_positions_by_degeneracy.py for all chromosomes first."
    exit 1
fi

echo ""
echo "All degeneracy annotations ready."
echo ""

# ============================================================================
# Step 2: Calculate π by compartment
# ============================================================================

echo "Step 2: Calculating π by genomic compartment..."
echo ""

OUTPUT_FILE="$OUTPUT_DIR/pi_by_compartment.txt"

# Check VCF format
if [ -f "$VCF" ]; then
    echo "Reading VCF from file: $VCF"
    
    if [[ "$VCF" == *.gz ]]; then
        zcat "$VCF" | python3 "${SCRIPT_DIR}/calculate_pi_by_genomic_compartment.py" \
            --stream \
            --gff "$GFF3" \
            --degeneracy $DEGENERACY_FILES \
            --output "$OUTPUT_FILE"
    else
        python3 "${SCRIPT_DIR}/calculate_pi_by_genomic_compartment.py" \
            --vcf "$VCF" \
            --gff "$GFF3" \
            --degeneracy $DEGENERACY_FILES \
            --output "$OUTPUT_FILE"
    fi
else
    echo "ERROR: VCF file not found: $VCF"
    exit 1
fi

echo ""

# ============================================================================
# Step 3: Generate summary report
# ============================================================================

echo "Step 3: Generating summary report..."

REPORT_FILE="$OUTPUT_DIR/pi_compartment_report.txt"

cat > "$REPORT_FILE" << EOF
================================================================================
NUCLEOTIDE DIVERSITY (π) BY GENOMIC COMPARTMENT
Mimulus guttatus IM767 Reference Genome
================================================================================

Analysis Date: $(date)
GFF3: $GFF3
VCF: $VCF

CONFIGURATION:
- Intergenic window size: 50 kb
- Intron trimming: 50 bp from splice sites
- Exon sites: 4-fold degenerate only
- Minimum sample size: 10

RESULTS:
EOF

# Append the results
if [ -f "$OUTPUT_FILE" ]; then
    cat "$OUTPUT_FILE" >> "$REPORT_FILE"
fi

cat >> "$REPORT_FILE" << EOF

--------------------------------------------------------------------------------
INTERPRETATION:

1. INTERGENIC: Represents the neutral baseline (no coding constraints,
   minimal selection). π should be determined primarily by mutation rate
   and demographic history.

2. INTRON: Neutral sites within genes. May be slightly lower than intergenic
   due to linked selection (background selection, hitchhiking from coding
   sites), but should be close to neutral.

3. FIRST EXON (4-fold): 4-fold degenerate sites in first exons. These are
   synonymous positions. May show different patterns if first exons have
   distinct characteristics (e.g., GC content, expression-linked selection).

4. NON-FIRST EXON (4-fold): 4-fold degenerate sites in exons 2+. These
   represent the bulk of synonymous sites.

EXPECTED PATTERN UNDER SELECTION + GC-BIASED GENE CONVERSION:
- If C/G-biased forces (selection + gBGC) increase π at synonymous sites,
  expect: π(4-fold exon) > π(intron) ≈ π(intergenic)
  
- The "hump" effect from weak selection would manifest as elevated π
  at sites under moderate translational selection.
================================================================================
EOF

echo "Report saved to: $REPORT_FILE"
echo ""

# ============================================================================
# Done
# ============================================================================

echo "========================================="
echo "Analysis Complete"
echo "========================================="
echo "Finished: $(date)"
echo ""
echo "Output files:"
echo "  - $OUTPUT_FILE"
echo "  - $REPORT_FILE"
echo ""
