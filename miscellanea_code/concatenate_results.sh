#!/bin/bash
#
# concatenate_results.sh
#
# Concatenate per-chromosome results into genome-wide files.
# Run this after all parallel jobs complete.
#
# Usage:
#   bash concatenate_results.sh
#
# Author: Luis Javier Madrigal-Roca & John K. Kelly

set -e  # Exit on error
set -u  # Exit on undefined variable

echo "========================================="
echo "Concatenating Polymorphism Results"
echo "========================================="
echo ""

# Configuration
OUTPUT_DIR="results/polymorphism_analysis"
CHROMOSOMES=(Chr_01 Chr_02 Chr_03 Chr_04 Chr_05 Chr_06 Chr_07 Chr_08 Chr_09 Chr_10 Chr_11 Chr_12 Chr_13 Chr_14)

# Check if output directory exists
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "ERROR: Output directory not found: $OUTPUT_DIR"
    echo "Have you run the chromosome processing jobs yet?"
    exit 1
fi

# ============================================================================
# Concatenate π statistics
# ============================================================================

echo "Creating all_chromosomes.bygene.pi.txt..."

FIRST=1
MISSING_COUNT=0

for CHR in "${CHROMOSOMES[@]}"; do
    PI_FILE="${OUTPUT_DIR}/${CHR}.bygene.pi.txt"
    
    if [ ! -f "$PI_FILE" ]; then
        echo "  ⚠ Warning: $PI_FILE not found, skipping..."
        MISSING_COUNT=$((MISSING_COUNT + 1))
        continue
    fi
    
    if [ $FIRST -eq 1 ]; then
        # Include header from first file
        cat "$PI_FILE" > "${OUTPUT_DIR}/all_chromosomes.bygene.pi.txt"
        FIRST=0
        echo "  ✓ Added $CHR (with header)"
    else
        # Skip header for subsequent files
        tail -n +2 "$PI_FILE" >> "${OUTPUT_DIR}/all_chromosomes.bygene.pi.txt"
        echo "  ✓ Added $CHR"
    fi
done

if [ $FIRST -eq 1 ]; then
    echo "  ✗ ERROR: No π files found!"
    exit 1
elif [ $MISSING_COUNT -gt 0 ]; then
    echo "  ⚠ Warning: $MISSING_COUNT chromosome(s) missing"
else
    echo "  ✓ All chromosomes concatenated successfully"
fi

# ============================================================================
# Concatenate codon frequencies
# ============================================================================

echo ""
echo "Checking for codon frequency files..."

# Check if any codon frequency files exist
CODON_FILES_EXIST=0
for CHR in "${CHROMOSOMES[@]}"; do
    if [ -f "${OUTPUT_DIR}/${CHR}.codon_frequencies.txt" ]; then
        CODON_FILES_EXIST=1
        break
    fi
done

if [ $CODON_FILES_EXIST -eq 0 ]; then
    echo "  ⚠ No codon frequency files found (preferred_codons.txt was missing?)"
    echo "  Skipping codon-level analyses..."
else
    echo "  ✓ Codon frequency files found"
    
    # Concatenate standard codon frequencies
    echo ""
    echo "Creating all_chromosomes.codon_frequencies.txt..."
    
    FIRST=1
    MISSING_COUNT=0
    
    for CHR in "${CHROMOSOMES[@]}"; do
        CODON_FILE="${OUTPUT_DIR}/${CHR}.codon_frequencies.txt"
        
        if [ ! -f "$CODON_FILE" ]; then
            echo "  ⚠ Warning: $CODON_FILE not found, skipping..."
            MISSING_COUNT=$((MISSING_COUNT + 1))
            continue
        fi
        
        if [ $FIRST -eq 1 ]; then
            cat "$CODON_FILE" > "${OUTPUT_DIR}/all_chromosomes.codon_frequencies.txt"
            FIRST=0
            echo "  ✓ Added $CHR (with header)"
        else
            tail -n +2 "$CODON_FILE" >> "${OUTPUT_DIR}/all_chromosomes.codon_frequencies.txt"
            echo "  ✓ Added $CHR"
        fi
    done
    
    if [ $FIRST -eq 1 ]; then
        echo "  ✗ ERROR: No codon frequency files found!"
    elif [ $MISSING_COUNT -gt 0 ]; then
        echo "  ⚠ Warning: $MISSING_COUNT chromosome(s) missing"
    else
        echo "  ✓ All chromosomes concatenated successfully"
    fi
    
    # Concatenate preferred-centric codon frequencies
    echo ""
    echo "Creating all_chromosomes.codon_frequencies_preferred.txt..."
    
    FIRST=1
    MISSING_COUNT=0
    
    for CHR in "${CHROMOSOMES[@]}"; do
        CODON_PREF_FILE="${OUTPUT_DIR}/${CHR}.codon_frequencies_preferred.txt"
        
        if [ ! -f "$CODON_PREF_FILE" ]; then
            echo "  ⚠ Warning: $CODON_PREF_FILE not found, skipping..."
            MISSING_COUNT=$((MISSING_COUNT + 1))
            continue
        fi
        
        if [ $FIRST -eq 1 ]; then
            cat "$CODON_PREF_FILE" > "${OUTPUT_DIR}/all_chromosomes.codon_frequencies_preferred.txt"
            FIRST=0
            echo "  ✓ Added $CHR (with header)"
        else
            tail -n +2 "$CODON_PREF_FILE" >> "${OUTPUT_DIR}/all_chromosomes.codon_frequencies_preferred.txt"
            echo "  ✓ Added $CHR"
        fi
    done
    
    if [ $FIRST -eq 1 ]; then
        echo "  ✗ ERROR: No preferred-centric codon files found!"
    elif [ $MISSING_COUNT -gt 0 ]; then
        echo "  ⚠ Warning: $MISSING_COUNT chromosome(s) missing"
    else
        echo "  ✓ All chromosomes concatenated successfully"
    fi
    
    # Concatenate site-level frequencies by preference
    echo ""
    echo "Creating all_chromosomes.site_freq_by_preference.txt..."
    
    FIRST=1
    MISSING_COUNT=0
    
    for CHR in "${CHROMOSOMES[@]}"; do
        SITE_FREQ_FILE="${OUTPUT_DIR}/${CHR}.site_freq_by_preference.txt"
        
        if [ ! -f "$SITE_FREQ_FILE" ]; then
            echo "  ⚠ Warning: $SITE_FREQ_FILE not found, skipping..."
            MISSING_COUNT=$((MISSING_COUNT + 1))
            continue
        fi
        
        if [ $FIRST -eq 1 ]; then
            cat "$SITE_FREQ_FILE" > "${OUTPUT_DIR}/all_chromosomes.site_freq_by_preference.txt"
            FIRST=0
            echo "  ✓ Added $CHR (with header)"
        else
            tail -n +2 "$SITE_FREQ_FILE" >> "${OUTPUT_DIR}/all_chromosomes.site_freq_by_preference.txt"
            echo "  ✓ Added $CHR"
        fi
    done
    
    if [ $FIRST -eq 1 ]; then
        echo "  ✗ ERROR: No site frequency files found!"
    elif [ $MISSING_COUNT -gt 0 ]; then
        echo "  ⚠ Warning: $MISSING_COUNT chromosome(s) missing"
    else
        echo "  ✓ All chromosomes concatenated successfully"
    fi
fi

# ============================================================================
# Generate summary statistics
# ============================================================================

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

echo ""
echo "========================================="
echo "Concatenation Complete!"
echo "========================================="
echo ""
echo "Output files in: $OUTPUT_DIR/"
echo ""
echo "Genome-wide files created:"
echo "  - all_chromosomes.bygene.pi.txt"
if [ $CODON_FILES_EXIST -eq 1 ]; then
    echo "  - all_chromosomes.codon_frequencies.txt"
    echo "  - all_chromosomes.codon_frequencies_preferred.txt"
    echo "  - all_chromosomes.site_freq_by_preference.txt"
fi
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
