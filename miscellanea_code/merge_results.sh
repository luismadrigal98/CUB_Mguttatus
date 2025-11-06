#!/bin/bash
#SBATCH --job-name=merger_codon_analysis
#SBATCH --output=logs/merger_codon_analysis.out
#SBATCH --error=logs/merger_codon_analysis.err
#SBATCH --time=12-00:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=kucg,eeb,kelly

################################################################################
# Merge Results from Parallel Chromosome Processing
#
# Run this after all chromosome jobs complete to concatenate results.
#
# Usage:
#   bash merge_results.sh
#
# Author: Luis Javier Madrigal-Roca & John K. Kelly
################################################################################

RESULTS_DIR="results/polymorphism"

echo "Merging results from all chromosomes..."

# Check if all chromosome files exist
MISSING_FILES=0
for CHR in Chr_{01..14}; do
    if [ ! -f "${RESULTS_DIR}/${CHR}.bygene.pi.txt" ]; then
        echo "WARNING: Missing ${CHR}.bygene.pi.txt"
        MISSING_FILES=$((MISSING_FILES + 1))
    fi
done

if [ $MISSING_FILES -gt 0 ]; then
    echo "ERROR: $MISSING_FILES chromosome files are missing!"
    echo "Check job logs in logs/ directory"
    exit 1
fi

echo "All chromosome files present. Merging..."

# Merge diversity metrics (π, Tajima's D, etc.)
echo "Merging diversity metrics..."
cat ${RESULTS_DIR}/Chr_01.bygene.pi.txt | head -1 > ${RESULTS_DIR}/all_chromosomes.bygene.pi.txt
for CHR in Chr_{01..14}; do
    grep -v "^Chr" ${RESULTS_DIR}/${CHR}.bygene.pi.txt >> ${RESULTS_DIR}/all_chromosomes.bygene.pi.txt
done
echo "✓ Created: ${RESULTS_DIR}/all_chromosomes.bygene.pi.txt"

# Merge codon frequency data
echo "Merging codon frequencies..."
cat ${RESULTS_DIR}/Chr_01.codon_frequencies.txt | head -1 > ${RESULTS_DIR}/all_chromosomes.codon_frequencies.txt
for CHR in Chr_{01..14}; do
    grep -v "^Gene" ${RESULTS_DIR}/${CHR}.codon_frequencies.txt >> ${RESULTS_DIR}/all_chromosomes.codon_frequencies.txt
done
echo "✓ Created: ${RESULTS_DIR}/all_chromosomes.codon_frequencies.txt"

# Generate summary statistics
echo ""
echo "Summary Statistics:"
echo "==================="

TOTAL_GENES=$(grep -v "^Chr" ${RESULTS_DIR}/all_chromosomes.bygene.pi.txt | wc -l)
echo "Total genes analyzed: $TOTAL_GENES"

TOTAL_CODONS=$(grep -v "^Gene" ${RESULTS_DIR}/all_chromosomes.codon_frequencies.txt | wc -l)
echo "Total codon positions: $TOTAL_CODONS"

echo ""
echo "Diversity by chromosome:"
echo "------------------------"
for CHR in Chr_{01..14}; do
    GENES=$(grep -v "^Chr" ${RESULTS_DIR}/${CHR}.bygene.pi.txt | wc -l)
    printf "%s: %6d genes\n" "$CHR" "$GENES"
done

echo ""
echo "Files ready for R analysis:"
echo "  - ${RESULTS_DIR}/all_chromosomes.bygene.pi.txt"
echo "  - ${RESULTS_DIR}/all_chromosomes.codon_frequencies.txt"
echo ""
echo "Next step: Run integrate_polymorphism_analysis.R"
