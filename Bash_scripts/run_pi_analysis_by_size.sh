#!/bin/bash
# SLURM wrapper to compute mean nucleotide diversity (π) for exons by size bins
# Usage: sbatch Bash_scripts/run_pi_analysis_by_size.sh <vcf.gz> <gff3> <results_dir>

#SBATCH --job-name=exonic_pi_size
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=exonic_pi.%j.out
#SBATCH --error=exonic_pi.%j.err

set -euo pipefail

VCF_GZ=${1:-}
GFF3=${2:-}
RESULTS_DIR=${3:-.}

if [[ -z "$VCF_GZ" || -z "$GFF3" ]]; then
    echo "Usage: $0 <vcf.gz> <gff3> <results_dir>"
    exit 1
fi

# Ensure results dir exists
mkdir -p "$RESULTS_DIR"

# Output files
EXONIC_PI_CSV="$RESULTS_DIR/exonic_pi.csv"
SIZE_BIN_CSV="$RESULTS_DIR/exonic_pi_by_size_bin.csv"

echo "========================================"
echo "Exonic π Analysis by Size Bin"
echo "========================================"
echo "Job: ${SLURM_JOB_ID:-local}"
echo "VCF: $VCF_GZ"
echo "GFF3: $GFF3"
echo "Results: $RESULTS_DIR"
echo "Start: $(date)"

# Step 1: Compute per-exon π
echo ""
echo "Step 1: Computing π for each exon..."
zcat "$VCF_GZ" | python3 miscellanea_code/calculate_exonic_pi.py /dev/stdin "$GFF3" "$EXONIC_PI_CSV" \
    --buffer-mb 64

echo "Per-exon π saved to: $EXONIC_PI_CSV"

# Step 2: Aggregate by size bins
echo ""
echo "Step 2: Aggregating π by size bins [1-1000), [1000-2000), ..., [10000+]..."
Rscript src/aggregate_exonic_pi_by_size_bin.R "$EXONIC_PI_CSV" "$SIZE_BIN_CSV"

echo "Aggregated results saved to: $SIZE_BIN_CSV"

# Step 3: Summary
echo ""
echo "========================================"
echo "Summary"
echo "========================================"
echo "Per-exon CSV: $EXONIC_PI_CSV"
wc -l "$EXONIC_PI_CSV"
echo ""
echo "By-size-bin CSV: $SIZE_BIN_CSV"
wc -l "$SIZE_BIN_CSV"
echo ""
head -20 "$SIZE_BIN_CSV"

echo ""
echo "Finished: $(date)"
echo "========================================"

exit 0
