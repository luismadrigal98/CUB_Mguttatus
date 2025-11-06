#!/bin/bash
#SBATCH --job-name=codon_polymorphism
#SBATCH --output=logs/codon_poly_%A_%a.out
#SBATCH --error=logs/codon_poly_%A_%a.err
#SBATCH --array=1-14
#SBATCH --time=12-00:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=kucg,eeb,kelly

################################################################################
# HPC-Optimized Parallel Processing for Codon-Aware Polymorphism Analysis
#
# This script processes all 14 chromosomes in parallel using SLURM job arrays.
# Each chromosome runs as an independent job, maximizing throughput.
#
# Usage:
#   sbatch run_parallel_hpc.sh
#
# Requirements:
#   - Python 3.7+
#   - Input files in data/ directory
#   - Output directory: results/polymorphism/
#
# Author: Luis Javier Madrigal-Roca & John K. Kelly
################################################################################

# Create output directories
mkdir -p logs
mkdir -p results/polymorphism

# Define chromosomes
CHROMOSOMES=(Chr_01 Chr_02 Chr_03 Chr_04 Chr_05 Chr_06 Chr_07 Chr_08 Chr_09 Chr_10 Chr_11 Chr_12 Chr_13 Chr_14)

# Get chromosome for this array task
CHR=${CHROMOSOMES[$SLURM_ARRAY_TASK_ID-1]}

echo "========================================="
echo "Processing: $CHR"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURM_NODELIST"
echo "Started: $(date)"
echo "========================================="

# Set paths
WORK_DIR="/home/l338m483/scratch/IMLines_to_767"
SCRIPT_DIR="miscellanea_code"
DATA_DIR="data"
RESULTS_DIR="results/polymorphism"
VCF_FILE="Included_snp.sites.txt"

# Change to working directory
cd $WORK_DIR

GFF3="${DATA_DIR}/Mguttatusvar_IM767_887_v2.1.gene.gff3"
GENOME="${DATA_DIR}/Mguttatusvar_IM767_887_v2.0.fa"
CDS="${DATA_DIR}/Mguttatusvar_IM767_887_v2.1.cds_primaryTranscriptOnly.fa"
PREFERRED_CODONS="preferred_codons.txt"
PREFILTERED_VCF_DIR="vcf_by_chromosome"

# Check if VCF is already filtered for this chromosome (recommended for speed)
if [ -f "${PREFILTERED_VCF_DIR}/${VCF_FILE}.${CHR}.vcf" ]; then
    VCF_INPUT="${PREFILTERED_VCF_DIR}/${VCF_FILE}.${CHR}.vcf"
    echo "Using pre-filtered VCF: $VCF_INPUT"
    echo "  File size: $(du -h $VCF_INPUT | cut -f1)"
else
    VCF_INPUT="$VCF_FILE"
    echo "Using full VCF: $VCF_INPUT (consider pre-filtering for speed)"
fi

# Step 1: Annotate genomic positions by degeneracy
echo ""
echo "Step 1: Annotating positions for $CHR..."
python3 ${SCRIPT_DIR}/describe_gene_positions_by_degeneracy.py \
    $CHR \
    $GFF3 \
    $GENOME \
    $CDS \
    > ${RESULTS_DIR}/${CHR}.annotation.log 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: Position annotation failed for $CHR"
    exit 1
fi

ANNOTATION_FILE="${CHR}.genic_bases.annotated.txt"
echo "✓ Created: $ANNOTATION_FILE"

# Move annotation file to results directory
mv $ANNOTATION_FILE ${RESULTS_DIR}/

# Step 2: Calculate nucleotide diversity (π)
echo ""
echo "Step 2: Calculating diversity metrics for $CHR..."
python3 ${SCRIPT_DIR}/calculate_pi.py \
    $CHR \
    $VCF_INPUT \
    ${RESULTS_DIR}/$ANNOTATION_FILE \
    187 \
    > ${RESULTS_DIR}/${CHR}.pi_calculation.log 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: Pi calculation failed for $CHR"
    exit 1
fi

PI_FILE="${CHR}.bygene.pi.txt"
echo "✓ Created: $PI_FILE"

# Move pi file to results directory
mv $PI_FILE ${RESULTS_DIR}/

# Step 3: Analyze preferred codon frequencies
echo ""
echo "Step 3: Analyzing codon frequencies for $CHR..."
python3 ${SCRIPT_DIR}/calculate_freq_preferred_codon.py \
    $CHR \
    $VCF_INPUT \
    $GFF3 \
    $CDS \
    $GENOME \
    $PREFERRED_CODONS \
    187 \
    > ${RESULTS_DIR}/${CHR}.codon_freq.log 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: Codon frequency analysis failed for $CHR"
    exit 1
fi

CODON_FILE="${CHR}.codon_frequencies.txt"
echo "✓ Created: $CODON_FILE"

# Move codon file to results directory
mv $CODON_FILE ${RESULTS_DIR}/

echo ""
echo "========================================="
echo "Completed: $CHR"
echo "Finished: $(date)"
echo "========================================="
