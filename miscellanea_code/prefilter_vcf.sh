#!/bin/bash
################################################################################
# PARALLEL Pre-filter VCF by Chromosome for Faster Processing
#
# This script splits the massive 407 GB VCF into chromosome-specific files.
# Uses SLURM array to process all chromosomes in PARALLEL.
#
# Usage:
#   sbatch prefilter_vcf.sh
#
# Estimated time: 2-4 hours (all chromosomes at once!)
# Memory required: 8 GB per job
# Total jobs: 14 (one per chromosome)
#
# Author: Luis Javier Madrigal-Roca & John K. Kelly
################################################################################

#SBATCH --job-name=vcf_prefilter
#SBATCH --output=logs/vcf_prefilter_%A_%a.out
#SBATCH --error=logs/vcf_prefilter_%A_%a.err
#SBATCH --array=1-14
#SBATCH --time=06:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=eeb,kelly,kucg

WORK_DIR="/home/l338m483/scratch/IMLines_to_767"

cd $WORK_DIR

mkdir -p logs
mkdir -p vcf_by_chromosome

VCF_FILE="Included_snp.sites.txt"

# Define chromosomes
CHROMOSOMES=(Chr_01 Chr_02 Chr_03 Chr_04 Chr_05 Chr_06 Chr_07 Chr_08 Chr_09 Chr_10 Chr_11 Chr_12 Chr_13 Chr_14)

# Get chromosome for this array task
CHR=${CHROMOSOMES[$SLURM_ARRAY_TASK_ID-1]}

echo "========================================="
echo "Pre-filtering VCF for $CHR"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Started: $(date)"
echo "========================================="

# Check if VCF exists and create header file (only first job)
if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
    if [ ! -f "$VCF_FILE" ]; then
        echo "ERROR: VCF file not found: $VCF_FILE"
        exit 1
    fi
    
    echo "Checking for VCF header..."
    HEADER_LINES=$(head -100 $VCF_FILE | grep -c "^#" || true)
    echo "Header lines: $HEADER_LINES"
    
    if [ $HEADER_LINES -gt 0 ]; then
        grep "^#" $VCF_FILE > vcf_by_chromosome/header.txt
    else
        touch vcf_by_chromosome/header.txt
    fi
else
    # Wait for first job to create header file
    WAIT_COUNT=0
    while [ ! -f "vcf_by_chromosome/header.txt" ] && [ $WAIT_COUNT -lt 60 ]; do
        sleep 2
        WAIT_COUNT=$((WAIT_COUNT + 1))
    done
fi

echo "Processing $CHR..."
OUTPUT="vcf_by_chromosome/${VCF_FILE}.${CHR}.vcf"

# Copy header and extract chromosome
cp vcf_by_chromosome/header.txt $OUTPUT
grep "^${CHR}" $VCF_FILE >> $OUTPUT

# Get stats
LINES=$(wc -l < $OUTPUT)
SIZE=$(du -h $OUTPUT | cut -f1)

echo "✓ $CHR: $LINES lines, $SIZE"
echo "Finished: $(date)"
