#!/bin/bash
#SBATCH --job-name=codon_polymorphism
#SBATCH --output=logs/codon_poly_%A_%a.out
#SBATCH --error=logs/codon_poly_%A_%a.err
#SBATCH --array=1-14
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=kucg,eeb,kelly
#
# run_polymorphism_analysis.sh
#
# Wrapper script to run codon-aware polymorphism analysis on all chromosomes.
# 
# Usage:
#   Sequential (local):  bash run_polymorphism_analysis.sh
#   Parallel (HPC):      sbatch run_polymorphism_analysis.sh
#
# Requirements:
#   - Python 3.x
#   - Input files in expected locations
#   - Preferred codons list (create from R analysis first)
#
# Author: Luis Javier Madrigal-Roca & John K. Kelly

set -e  # Exit on error
set -u  # Exit on undefined variable

# ============================================================================
# Configuration
# ============================================================================

# Detect if running in SLURM environment (parallel) or standalone (sequential)
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    PARALLEL_MODE=1
    echo "Running in PARALLEL mode (SLURM array job)"
else
    PARALLEL_MODE=0
    echo "Running in SEQUENTIAL mode"
fi

# Input files
GFF3="data/Mguttatusvar_IM767_887_v2.1.gene.gff3"
GENOME_FA="data/Mguttatusvar_IM767_887_v2.0.fa"
CDS_FA="data/Mguttatusvar_IM767_887_v2.1.cds_primaryTranscriptOnly.fa"
VCF="Included_snp.sites.txt"
PREFERRED_CODONS="preferred_codons.txt"

# Chromosomes to process
CHROMOSOMES_ARRAY=(Chr_01 Chr_02 Chr_03 Chr_04 Chr_05 Chr_06 Chr_07 Chr_08 Chr_09 Chr_10 Chr_11 Chr_12 Chr_13 Chr_14)

# Python scripts directory
SCRIPT_DIR="miscellanea_code"

# Output directory
OUTPUT_DIR="results/polymorphism_analysis"

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# ============================================================================
# Determine chromosomes to process
# ============================================================================

if [ $PARALLEL_MODE -eq 1 ]; then
    # Parallel mode: process single chromosome based on array task ID
    CHR=${CHROMOSOMES_ARRAY[$SLURM_ARRAY_TASK_ID-1]}
    CHROMOSOMES_TO_PROCESS=($CHR)
    
    echo "========================================="
    echo "Codon-Aware Polymorphism Analysis"
    echo "PARALLEL MODE - Processing: $CHR"
    echo "========================================="
    echo "Job ID: $SLURM_JOB_ID"
    echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
    echo "Node: ${SLURM_NODELIST:-localhost}"
    echo "Started: $(date)"
    echo ""
else
    # Sequential mode: process all chromosomes
    CHROMOSOMES_TO_PROCESS=("${CHROMOSOMES_ARRAY[@]}")
    
    echo "========================================="
    echo "Codon-Aware Polymorphism Analysis"
    echo "SEQUENTIAL MODE - Processing all chromosomes"
    echo "========================================="
    echo "Started: $(date)"
    echo ""
fi

# ============================================================================
# Pre-flight checks
# ============================================================================

# Check if required files exist
echo "Checking input files..."

for FILE in "$GFF3" "$GENOME_FA" "$CDS_FA" "$VCF"; do
    if [ ! -f "$FILE" ]; then
        echo "ERROR: Required file not found: $FILE"
        exit 1
    fi
    echo "  ✓ $FILE"
done

# Check for preferred codons file
if [ ! -f "$PREFERRED_CODONS" ]; then
    echo ""
    echo "WARNING: Preferred codons file not found: $PREFERRED_CODONS"
    echo ""
    if [ $PARALLEL_MODE -eq 0 ]; then
        echo "Please create this file first from your R analysis:"
        echo ""
        echo "  preferred_codons <- cai_results\$w_table %>%"
        echo "    filter(relative_adaptiveness == 1.0) %>%"
        echo "    pull(codon)"
        echo "  writeLines(preferred_codons, \"$PREFERRED_CODONS\")"
        echo ""
    fi
    echo "Continuing without steps 3, 3b, 3c (codon frequency analysis)..."
    SKIP_STEP3=1
else
    echo "  ✓ $PREFERRED_CODONS"
    SKIP_STEP3=0
fi

# Check Python scripts
echo ""
echo "Checking Python scripts..."

SCRIPT1="$SCRIPT_DIR/describe_gene_positions_by_degeneracy.py"
SCRIPT2="$SCRIPT_DIR/calculate_pi.py"
SCRIPT3="$SCRIPT_DIR/calculate_freq_preferred_codon.py"
SCRIPT4="$SCRIPT_DIR/calculate_site_freq_by_preference.py"
SCRIPT5="$SCRIPT_DIR/reformat_codon_freq_by_preference.py"

for SCRIPT in "$SCRIPT1" "$SCRIPT2" "$SCRIPT3" "$SCRIPT4" "$SCRIPT5"; do
    if [ ! -f "$SCRIPT" ]; then
        echo "ERROR: Python script not found: $SCRIPT"
        exit 1
    fi
    echo "  ✓ $SCRIPT"
done

echo ""
echo "All checks passed! Starting analysis..."
echo ""

# ============================================================================
# Detect sample size from VCF (only once for speed)
# ============================================================================

echo "Detecting sample size from VCF..."
N_SAMPLES=$(grep -m1 "^#CHROM" "$VCF" | awk '{print NF-9}')

if [ -z "$N_SAMPLES" ] || [ "$N_SAMPLES" -lt 1 ]; then
    echo "ERROR: Could not determine sample size from VCF"
    exit 1
fi

echo "  ✓ Detected $N_SAMPLES samples in VCF"
echo ""

# ============================================================================
# Process each chromosome
# ============================================================================

START_TIME=$(date +%s)

for CHR in "${CHROMOSOMES_TO_PROCESS[@]}"; do
    echo "========================================="
    echo "Processing $CHR"
    echo "========================================="
    
    CHR_START=$(date +%s)
    
    # ------------------------------------------------------------------------
    # Step 1: Annotate positions by degeneracy
    # ------------------------------------------------------------------------
    
    echo ""
    echo "Step 1: Annotating genomic positions..."
    
    ANNOTATION_FILE="${OUTPUT_DIR}/${CHR}.genic_bases.annotated.txt"
    
    if [ -f "$ANNOTATION_FILE" ]; then
        echo "  ⚠ Annotation file already exists, skipping..."
    else
        python3 "$SCRIPT1" \
            "$CHR" \
            "$GFF3" \
            "$GENOME_FA" \
            "$CDS_FA" \
            > "${OUTPUT_DIR}/${CHR}.step1.log" 2>&1
        
        # Move output to results directory
        if [ -f "${CHR}.genic_bases.annotated.txt" ]; then
            mv "${CHR}.genic_bases.annotated.txt" "$ANNOTATION_FILE"
        fi
        
        echo "  ✓ Created $ANNOTATION_FILE"
    fi
    
    # ------------------------------------------------------------------------
    # Step 2: Calculate π and diversity statistics
    # ------------------------------------------------------------------------
    
    echo ""
    echo "Step 2: Calculating nucleotide diversity..."
    
    PI_FILE="${OUTPUT_DIR}/${CHR}.bygene.pi.txt"
    
    if [ -f "$PI_FILE" ]; then
        echo "  ⚠ Pi file already exists, skipping..."
    else
        python3 "$SCRIPT2" \
            "$CHR" \
            "$VCF" \
            "$ANNOTATION_FILE" \
            "$N_SAMPLES" \
            > "${OUTPUT_DIR}/${CHR}.step2.log" 2>&1
        
        # Move output to results directory
        if [ -f "${CHR}.bygene.pi.txt" ]; then
            mv "${CHR}.bygene.pi.txt" "$PI_FILE"
        fi
        
        echo "  ✓ Created $PI_FILE"
    fi
    
    # ------------------------------------------------------------------------
    # Step 3: Analyze preferred codon frequencies (if preferred codons exist)
    # ------------------------------------------------------------------------
    
    if [ $SKIP_STEP3 -eq 0 ]; then
        echo ""
        echo "Step 3: Analyzing codon-level frequencies..."
        
        CODON_FILE="${OUTPUT_DIR}/${CHR}.codon_frequencies.txt"
        
        if [ -f "$CODON_FILE" ]; then
            echo "  ⚠ Codon frequency file already exists, skipping..."
        else
            python3 "$SCRIPT3" \
                "$CHR" \
                "$VCF" \
                "$GFF3" \
                "$CDS_FA" \
                "$GENOME_FA" \
                "$PREFERRED_CODONS" \
                "$N_SAMPLES" \
                > "${OUTPUT_DIR}/${CHR}.step3.log" 2>&1
            
            # Move output to results directory
            if [ -f "${CHR}.codon_frequencies.txt" ]; then
                mv "${CHR}.codon_frequencies.txt" "$CODON_FILE"
            fi
            
            echo "  ✓ Created $CODON_FILE"
        fi
        
        # Reformat to preferred-centric view
        echo ""
        echo "Step 3b: Reformatting codons by preference..."
        
        CODON_PREF_FILE="${OUTPUT_DIR}/${CHR}.codon_frequencies_preferred.txt"
        
        if [ -f "$CODON_PREF_FILE" ]; then
            echo "  ⚠ Preferred-centric file already exists, skipping..."
        else
            python3 "$SCRIPT5" \
                "$CODON_FILE" \
                "$PREFERRED_CODONS" \
                "$CODON_PREF_FILE" \
                > "${OUTPUT_DIR}/${CHR}.step3b.log" 2>&1
            
            echo "  ✓ Created $CODON_PREF_FILE"
        fi
        
        # Calculate site-level frequencies by preference and degeneracy
        echo ""
        echo "Step 3c: Calculating site-level frequencies by preference..."
        
        SITE_FREQ_FILE="${OUTPUT_DIR}/${CHR}.site_freq_by_preference.txt"
        
        if [ -f "$SITE_FREQ_FILE" ]; then
            echo "  ⚠ Site frequency file already exists, skipping..."
        else
            python3 "$SCRIPT4" \
                "$CHR" \
                "$VCF" \
                "$ANNOTATION_FILE" \
                "$PREFERRED_CODONS" \
                > "${OUTPUT_DIR}/${CHR}.step3c.log" 2>&1
            
            # Move output to results directory
            if [ -f "${CHR}.site_freq_by_preference.txt" ]; then
                mv "${CHR}.site_freq_by_preference.txt" "$SITE_FREQ_FILE"
            fi
            
            echo "  ✓ Created $SITE_FREQ_FILE"
        fi
    fi
    
    CHR_END=$(date +%s)
    CHR_ELAPSED=$((CHR_END - CHR_START))
    
    echo ""
    echo "✓ $CHR complete! ($(($CHR_ELAPSED / 60)) min $(($CHR_ELAPSED % 60)) sec)"
    echo ""
done

# ============================================================================
# Concatenate results across chromosomes (sequential mode only)
# ============================================================================

if [ $PARALLEL_MODE -eq 1 ]; then
    echo ""
    echo "========================================="
    echo "Completed: $CHR"
    echo "Finished: $(date)"
    echo "========================================="
    echo ""
    echo "After all array jobs complete, run concatenation step:"
    echo "  bash ${SCRIPT_DIR}/concatenate_results.sh"
    exit 0
fi

echo "========================================="
echo "Concatenating results across chromosomes"
echo "========================================="

# Concatenate π statistics
echo ""
echo "Creating all_chromosomes.bygene.pi.txt..."

FIRST=1
for CHR in "${CHROMOSOMES_ARRAY[@]}"; do
    PI_FILE="${OUTPUT_DIR}/${CHR}.bygene.pi.txt"
    
    if [ ! -f "$PI_FILE" ]; then
        echo "  ⚠ Warning: $PI_FILE not found, skipping..."
        continue
    fi
    
    if [ $FIRST -eq 1 ]; then
        # Include header from first file
        cat "$PI_FILE" > "${OUTPUT_DIR}/all_chromosomes.bygene.pi.txt"
        FIRST=0
    else
        # Skip header for subsequent files
        tail -n +2 "$PI_FILE" >> "${OUTPUT_DIR}/all_chromosomes.bygene.pi.txt"
    fi
done

echo "  ✓ Created ${OUTPUT_DIR}/all_chromosomes.bygene.pi.txt"

# Concatenate codon frequencies
if [ $SKIP_STEP3 -eq 0 ]; then
    echo ""
    echo "Creating all_chromosomes.codon_frequencies.txt..."
    
    FIRST=1
    for CHR in "${CHROMOSOMES_ARRAY[@]}"; do
        CODON_FILE="${OUTPUT_DIR}/${CHR}.codon_frequencies.txt"
        
        if [ ! -f "$CODON_FILE" ]; then
            echo "  ⚠ Warning: $CODON_FILE not found, skipping..."
            continue
        fi
        
        if [ $FIRST -eq 1 ]; then
            cat "$CODON_FILE" > "${OUTPUT_DIR}/all_chromosomes.codon_frequencies.txt"
            FIRST=0
        else
            tail -n +2 "$CODON_FILE" >> "${OUTPUT_DIR}/all_chromosomes.codon_frequencies.txt"
        fi
    done
    
    echo "  ✓ Created ${OUTPUT_DIR}/all_chromosomes.codon_frequencies.txt"
    
    # Concatenate preferred-centric codon frequencies
    echo ""
    echo "Creating all_chromosomes.codon_frequencies_preferred.txt..."
    
    FIRST=1
    for CHR in "${CHROMOSOMES_ARRAY[@]}"; do
        CODON_PREF_FILE="${OUTPUT_DIR}/${CHR}.codon_frequencies_preferred.txt"
        
        if [ ! -f "$CODON_PREF_FILE" ]; then
            echo "  ⚠ Warning: $CODON_PREF_FILE not found, skipping..."
            continue
        fi
        
        if [ $FIRST -eq 1 ]; then
            cat "$CODON_PREF_FILE" > "${OUTPUT_DIR}/all_chromosomes.codon_frequencies_preferred.txt"
            FIRST=0
        else
            tail -n +2 "$CODON_PREF_FILE" >> "${OUTPUT_DIR}/all_chromosomes.codon_frequencies_preferred.txt"
        fi
    done
    
    echo "  ✓ Created ${OUTPUT_DIR}/all_chromosomes.codon_frequencies_preferred.txt"
    
    # Concatenate site-level frequencies by preference
    echo ""
    echo "Creating all_chromosomes.site_freq_by_preference.txt..."
    
    FIRST=1
    for CHR in "${CHROMOSOMES_ARRAY[@]}"; do
        SITE_FREQ_FILE="${OUTPUT_DIR}/${CHR}.site_freq_by_preference.txt"
        
        if [ ! -f "$SITE_FREQ_FILE" ]; then
            echo "  ⚠ Warning: $SITE_FREQ_FILE not found, skipping..."
            continue
        fi
        
        if [ $FIRST -eq 1 ]; then
            cat "$SITE_FREQ_FILE" > "${OUTPUT_DIR}/all_chromosomes.site_freq_by_preference.txt"
            FIRST=0
        else
            tail -n +2 "$SITE_FREQ_FILE" >> "${OUTPUT_DIR}/all_chromosomes.site_freq_by_preference.txt"
        fi
    done
    
    echo "  ✓ Created ${OUTPUT_DIR}/all_chromosomes.site_freq_by_preference.txt"
fi

# ============================================================================
# Generate summary statistics
# ============================================================================

echo ""
echo "========================================="
echo "Summary Statistics"
echo "========================================="

PI_FILE="${OUTPUT_DIR}/all_chromosomes.bygene.pi.txt"

if [ -f "$PI_FILE" ]; then
    N_GENES=$(tail -n +2 "$PI_FILE" | wc -l)
    echo ""
    echo "Total genes analyzed: $N_GENES"
    
    echo ""
    echo "Average diversity per degeneracy class:"
    echo "  (averaged across all genes)"
    
    # Calculate column averages using awk
    tail -n +2 "$PI_FILE" | awk -F'\t' '
    BEGIN {
        sum1=0; sum2=0; sum3nf=0; sum3f=0; n=0
    }
    {
        if ($3 > 0) sum1 += $5/$3
        if ($6 > 0) sum2 += $8/$6
        if ($9 > 0) sum3nf += $11/$9
        if ($12 > 0) sum3f += $14/$12
        n++
    }
    END {
        printf "  π (1st position):     %.6f\n", sum1/n
        printf "  π (2nd position):     %.6f\n", sum2/n
        printf "  π (3rd non-4f):       %.6f\n", sum3nf/n
        printf "  π (3rd 4-fold):       %.6f\n", sum3f/n
    }'
fi

# ============================================================================
# Completion
# ============================================================================

END_TIME=$(date +%s)
TOTAL_ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "========================================="
echo "Analysis Complete!"
echo "========================================="
echo ""
echo "Total time: $(($TOTAL_ELAPSED / 60)) minutes"
echo ""
echo "Output files in: $OUTPUT_DIR/"
echo ""
echo "Next steps:"
echo "  1. Load results into R:"
echo "     polymorphism_data <- read.table('$OUTPUT_DIR/all_chromosomes.bygene.pi.txt', header=TRUE)"
echo ""
echo "  2. Merge with CUB data:"
echo "     cub_polymorphism <- left_join(exp_enc_data_cai, polymorphism_data, by=c('Gene_name'='Gene'))"
echo ""
echo "  3. Test hypotheses:"
echo "     cor.test(cub_polymorphism\$CAI, cub_polymorphism\$Pi_3_4f)"
echo ""
