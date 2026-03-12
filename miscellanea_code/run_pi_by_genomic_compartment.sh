#!/bin/bash
#SBATCH --job-name=pi_compartments
#SBATCH --output=logs/pi_compartments_%A_%a.out
#SBATCH --error=logs/pi_compartments_%A_%a.err
#SBATCH --array=1-14
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
#   sbatch run_pi_by_genomic_compartment.sh        # Parallel (HPC)
#   bash run_pi_by_genomic_compartment.sh          # Sequential (local)
#   bash run_pi_by_genomic_compartment.sh concat   # Concatenate after parallel jobs
#
# Author: Luis Javier Madrigal-Roca & GitHub Copilot
# Date: 2026-01-14

set -e
set -u

# ============================================================================
# Configuration
# ============================================================================

# Input files
GFF3="/home/l338m483/scratch/IMLines_to_767/data/Mguttatusvar_IM767_887_v2.1.gene.gff3"
GENOME_FA="/home/l338m483/scratch/IMLines_to_767/data/Mguttatusvar_IM767_887_v2.0.fa"
CDS_FA="/home/l338m483/scratch/IMLines_to_767/data/Mguttatusvar_IM767_887_v2.1.cds_primaryTranscriptOnly.fa"
VCF="/home/l338m483/scratch/IMLines_to_767/Included_snp.sites.txt"  # Your all-sites VCF

# Python scripts
SCRIPT_DIR="/home/l338m483/scratch/IMLines_to_767/python_scripts"
MISCELLANEA_DIR="/home/l338m483/scratch/IMLines_to_767/miscellanea_code"

# Output
OUTPUT_DIR="results/pi_compartment_analysis"
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# Chromosomes
CHROMOSOMES_ARRAY=(Chr_01 Chr_02 Chr_03 Chr_04 Chr_05 Chr_06 Chr_07 Chr_08 Chr_09 Chr_10 Chr_11 Chr_12 Chr_13 Chr_14)

# Pre-filtered VCF directory (optional, for speed) - USE ABSOLUTE PATH
# These should be named like: Chr_01.vcf, Chr_02.vcf, etc. OR Included_snp.sites.Chr_01.vcf
PREFILTERED_VCF_DIR="/home/l338m483/scratch/IMLines_to_767/vcf_by_chromosome"

# VCF base name (without path) for pattern matching
VCF_BASENAME=$(basename "$VCF")

# ============================================================================
# Detect execution mode
# ============================================================================

# Check for concatenation-only mode
if [ "${1:-}" = "concat" ]; then
    CONCAT_ONLY=1
else
    CONCAT_ONLY=0
fi

# Detect if running in SLURM environment (parallel) or standalone (sequential)
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    PARALLEL_MODE=1
    echo "Running in PARALLEL mode (SLURM array job)"
else
    PARALLEL_MODE=0
    if [ $CONCAT_ONLY -eq 1 ]; then
        echo "Running in CONCATENATION-ONLY mode"
    else
        echo "Running in SEQUENTIAL mode"
    fi
fi

echo "========================================="
echo "π by Genomic Compartment Analysis"
echo "========================================="
echo "Started: $(date)"
echo ""

# ============================================================================
# Determine chromosomes to process
# ============================================================================

if [ $PARALLEL_MODE -eq 1 ]; then
    # Parallel mode: process single chromosome based on array task ID
    CHR=${CHROMOSOMES_ARRAY[$SLURM_ARRAY_TASK_ID-1]}
    CHROMOSOMES_TO_PROCESS=($CHR)
    
    echo "PARALLEL MODE - Processing: $CHR"
    echo "Job ID: $SLURM_JOB_ID"
    echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
    echo "Node: ${SLURM_NODELIST:-localhost}"
    echo ""
elif [ $CONCAT_ONLY -eq 1 ]; then
    CHROMOSOMES_TO_PROCESS=()
    echo "Skipping processing, will concatenate only."
    echo ""
else
    # Sequential mode: process all chromosomes
    CHROMOSOMES_TO_PROCESS=("${CHROMOSOMES_ARRAY[@]}")
    echo "SEQUENTIAL MODE - Processing all chromosomes"
    echo ""
fi

# ============================================================================
# Step 1: Generate degeneracy annotations (if not already done)
# ============================================================================

if [ ${#CHROMOSOMES_TO_PROCESS[@]} -gt 0 ]; then
    echo "Step 1: Checking degeneracy annotation files..."
    
    DEGENERACY_FILES=""
    MISSING_ANNOT=0
    
    for CHR in "${CHROMOSOMES_TO_PROCESS[@]}"; do
        ANNOT_FILE="$OUTPUT_DIR/${CHR}.genic_bases.annotated.txt"
        
        if [ -f "$ANNOT_FILE" ]; then
            echo "  ✓ $CHR annotation exists"
        else
            echo "  Generating annotation for $CHR..."
            
            python3 "${MISCELLANEA_DIR}/describe_gene_positions_by_degeneracy.py" \
                "$CHR" \
                "$GFF3" \
                "$GENOME_FA" \
                "$CDS_FA" \
                > "${OUTPUT_DIR}/${CHR}.step1.log" 2>&1
            
            # Move to output directory
            if [ -f "${CHR}.genic_bases.annotated.txt" ]; then
                mv "${CHR}.genic_bases.annotated.txt" "$ANNOT_FILE"
                echo "  ✓ Created $ANNOT_FILE"
            fi
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
fi

# ============================================================================
# Step 1b: Pre-filter VCF by chromosome (if not already done)
# ============================================================================

if [ ${#CHROMOSOMES_TO_PROCESS[@]} -gt 0 ]; then
    echo "Step 1b: Checking for pre-filtered VCF files..."
    
    # Create directory if needed
    mkdir -p "$PREFILTERED_VCF_DIR"
    
    for CHR in "${CHROMOSOMES_TO_PROCESS[@]}"; do
        CHR_VCF="${PREFILTERED_VCF_DIR}/${CHR}.vcf"
        
        if [ -f "$CHR_VCF" ]; then
            echo "  ✓ $CHR pre-filtered VCF exists"
        else
            echo "  Extracting $CHR from VCF (this may take a while first time)..."
            # Use awk for speed - match chromosome in first column
            awk -v chr="$CHR" '$1 == chr' "$VCF" > "$CHR_VCF"
            
            if [ -s "$CHR_VCF" ]; then
                LINE_COUNT=$(wc -l < "$CHR_VCF")
                echo "  ✓ Created $CHR_VCF ($LINE_COUNT lines)"
            else
                echo "  ⚠ WARNING: No lines extracted for $CHR"
                rm -f "$CHR_VCF"
            fi
        fi
    done
    
    echo ""
fi

# ============================================================================
# Step 2: Calculate π by compartment (per chromosome)
# ============================================================================

if [ ${#CHROMOSOMES_TO_PROCESS[@]} -gt 0 ]; then
    echo "Step 2: Calculating π by genomic compartment..."
    echo ""
    
    for CHR in "${CHROMOSOMES_TO_PROCESS[@]}"; do
        echo "-----------------------------------------"
        echo "Processing $CHR"
        echo "-----------------------------------------"
        
        CHR_START=$(date +%s)
        
        OUTPUT_FILE="$OUTPUT_DIR/${CHR}.pi_by_compartment.txt"
        ANNOT_FILE="$OUTPUT_DIR/${CHR}.genic_bases.annotated.txt"
        
        # Check if output already exists
        if [ -f "$OUTPUT_FILE" ]; then
            echo "  ⚠ Output file already exists, skipping: $OUTPUT_FILE"
            continue
        fi
        
        # Check for pre-filtered VCF (much faster!)
        # Try multiple naming patterns
        VCF_INPUT=""
        
        # Pattern 1: {VCF_BASENAME}.{CHR}.vcf (e.g., Included_snp.sites.txt.Chr_03.vcf)
        if [ -z "$VCF_INPUT" ] && [ -f "${PREFILTERED_VCF_DIR}/${VCF_BASENAME}.${CHR}.vcf" ]; then
            VCF_INPUT="${PREFILTERED_VCF_DIR}/${VCF_BASENAME}.${CHR}.vcf"
        fi
        
        # Pattern 2: {CHR}.vcf (e.g., Chr_03.vcf)
        if [ -z "$VCF_INPUT" ] && [ -f "${PREFILTERED_VCF_DIR}/${CHR}.vcf" ]; then
            VCF_INPUT="${PREFILTERED_VCF_DIR}/${CHR}.vcf"
        fi
        
        # Pattern 3: {VCF_BASENAME_NO_EXT}.{CHR}.vcf (e.g., Included_snp.sites.Chr_03.vcf)
        VCF_NO_EXT="${VCF_BASENAME%.*}"
        if [ -z "$VCF_INPUT" ] && [ -f "${PREFILTERED_VCF_DIR}/${VCF_NO_EXT}.${CHR}.vcf" ]; then
            VCF_INPUT="${PREFILTERED_VCF_DIR}/${VCF_NO_EXT}.${CHR}.vcf"
        fi
        
        # Pattern 4: Check .vcf.gz variants
        if [ -z "$VCF_INPUT" ] && [ -f "${PREFILTERED_VCF_DIR}/${CHR}.vcf.gz" ]; then
            VCF_INPUT="${PREFILTERED_VCF_DIR}/${CHR}.vcf.gz"
        fi
        
        # Fallback to full VCF
        if [ -z "$VCF_INPUT" ]; then
            VCF_INPUT="$VCF"
            echo "  Using full VCF: $VCF_INPUT"
            echo "  ⚠ Pre-filtered VCF not found in: $PREFILTERED_VCF_DIR"
            echo "  ⚠ Tried patterns: ${VCF_BASENAME}.${CHR}.vcf, ${CHR}.vcf, ${VCF_NO_EXT}.${CHR}.vcf"
            echo "  ⚠ Consider pre-filtering VCF by chromosome for faster processing:"
            echo "      grep -E '^#|^${CHR}\\s' $VCF > ${PREFILTERED_VCF_DIR}/${CHR}.vcf"
        else
            echo "  Using pre-filtered VCF: $VCF_INPUT"
            echo "  File size: $(du -h "$VCF_INPUT" | cut -f1)"
        fi
        
        # Run the compartment analysis
        PERGENE_FILE="$OUTPUT_DIR/${CHR}.pi_per_gene_feature.txt"
        if [[ "$VCF_INPUT" == *.gz ]]; then
            zcat "$VCF_INPUT" | python3 "${SCRIPT_DIR}/calculate_pi_by_genomic_compartment.py" \
                --stream \
                --gff "$GFF3" \
                --degeneracy "$ANNOT_FILE" \
                --chromosome "$CHR" \
                --output "$OUTPUT_FILE" \
                --per-gene-output "$PERGENE_FILE" \
                > "${OUTPUT_DIR}/${CHR}.step2.log" 2>&1
        else
            python3 "${SCRIPT_DIR}/calculate_pi_by_genomic_compartment.py" \
                --vcf "$VCF_INPUT" \
                --gff "$GFF3" \
                --degeneracy "$ANNOT_FILE" \
                --chromosome "$CHR" \
                --output "$OUTPUT_FILE" \
                --per-gene-output "$PERGENE_FILE" \
                > "${OUTPUT_DIR}/${CHR}.step2.log" 2>&1
        fi
        
        CHR_END=$(date +%s)
        CHR_ELAPSED=$((CHR_END - CHR_START))
        
        if [ -f "$OUTPUT_FILE" ]; then
            echo "  ✓ Created $OUTPUT_FILE ($(($CHR_ELAPSED / 60)) min $(($CHR_ELAPSED % 60)) sec)"
        else
            echo "  ✗ ERROR: Failed to create output for $CHR"
        fi
        echo ""
    done
fi

# ============================================================================
# Exit early if parallel mode (concatenation done separately)
# ============================================================================

if [ $PARALLEL_MODE -eq 1 ]; then
    echo ""
    echo "========================================="
    echo "Completed: $CHR"
    echo "Finished: $(date)"
    echo "========================================="
    echo ""
    echo "After all array jobs complete, run concatenation step:"
    echo "  bash ${0} concat"
    exit 0
fi

# ============================================================================
# Step 3: Concatenate results across chromosomes
# ============================================================================

echo "========================================="
echo "Step 3: Concatenating results across chromosomes"
echo "========================================="

COMBINED_FILE="$OUTPUT_DIR/all_chromosomes.pi_by_compartment.txt"

FIRST=1
for CHR in "${CHROMOSOMES_ARRAY[@]}"; do
    CHR_FILE="$OUTPUT_DIR/${CHR}.pi_by_compartment.txt"
    
    if [ ! -f "$CHR_FILE" ]; then
        echo "  ⚠ Warning: $CHR_FILE not found, skipping..."
        continue
    fi
    
    if [ $FIRST -eq 1 ]; then
        # Include header from first file
        cat "$CHR_FILE" > "$COMBINED_FILE"
        FIRST=0
    else
        # Skip header for subsequent files
        tail -n +2 "$CHR_FILE" >> "$COMBINED_FILE"
    fi
    
    echo "  ✓ Added $CHR"
done

if [ -f "$COMBINED_FILE" ]; then
    echo ""
    echo "Created: $COMBINED_FILE"
else
    echo ""
    echo "ERROR: No chromosome files found to concatenate."
    exit 1
fi

# Concatenate per-gene feature files
COMBINED_PERGENE="$OUTPUT_DIR/all_chromosomes.pi_per_gene_feature.txt"

FIRST=1
for CHR in "${CHROMOSOMES_ARRAY[@]}"; do
    CHR_PERGENE="$OUTPUT_DIR/${CHR}.pi_per_gene_feature.txt"
    
    if [ ! -f "$CHR_PERGENE" ]; then
        continue
    fi
    
    if [ $FIRST -eq 1 ]; then
        cat "$CHR_PERGENE" > "$COMBINED_PERGENE"
        FIRST=0
    else
        tail -n +2 "$CHR_PERGENE" >> "$COMBINED_PERGENE"
    fi
    
    echo "  ✓ Added $CHR (per-gene feature)"
done

if [ -f "$COMBINED_PERGENE" ]; then
    echo ""
    echo "Created: $COMBINED_PERGENE"
fi

# ============================================================================
# Step 4: Generate summary report
# ============================================================================

echo ""
echo "Step 4: Generating summary report..."

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

# Append the combined results
if [ -f "$COMBINED_FILE" ]; then
    echo "" >> "$REPORT_FILE"
    cat "$COMBINED_FILE" >> "$REPORT_FILE"
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
echo "  - $COMBINED_FILE"
echo "  - $REPORT_FILE"
echo ""
