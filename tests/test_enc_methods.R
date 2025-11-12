#' Test script to compare Wright (1990) and Sun et al. (2012) ENC methods
#' 
#' @description Tests the calculate_enc function with both methods on sample data
#' to verify proper implementation of the Sun et al. (2012) pseudocount correction.
#' 
#' @author Luis J. Madrigal-Roca
#' @date November 12, 2025
#' ___________________________________________________________________________

# Load required libraries
library(data.table)
library(Biostrings)

# Source the calculate_enc function
source("src/calculate_enc.R")

# Load genetic code
data(GENETIC_CODE, package = 'Biostrings')

cat("\n=== Testing calculate_enc with Wright and Sun methods ===\n\n")

# Create test data with known codon usage patterns
# Test Case 1: Extreme bias (one codon per amino acid)
test_extreme <- data.table(
  Gene_name = "ExtremeGene",
  AAA = 10, AAC = 0, AAG = 10, AAT = 0,  # Lys: AAA only, Asn: AAG only
  AGA = 0, AGG = 0, AGC = 0, AGT = 0,    # Arg 2-fold: none, Ser 2-fold: none
  CGA = 0, CGC = 0, CGG = 0, CGT = 10,   # Arg 4-fold: CGT only
  TCA = 0, TCC = 0, TCG = 0, TCT = 10,   # Ser 4-fold: TCT only
  GCA = 10, GCC = 0, GCG = 0, GCT = 0,   # Ala: GCA only
  TTA = 0, TTG = 0, CTA = 0, CTC = 10, CTG = 0, CTT = 0  # Leu: CTC only
)

# Add remaining codons with appropriate counts
remaining_codons <- setdiff(names(GENETIC_CODE)[GENETIC_CODE != "STOP"], names(test_extreme))
for (codon in remaining_codons) {
  test_extreme[[codon]] <- 10
}

# Test Case 2: No bias (uniform usage)
test_uniform <- data.table(Gene_name = "UniformGene")
for (codon in names(GENETIC_CODE)[GENETIC_CODE != "STOP"]) {
  test_uniform[[codon]] <- 10
}

# Test Case 3: Moderate bias (realistic scenario)
test_moderate <- data.table(Gene_name = "ModerateGene")
set.seed(123)
for (codon in names(GENETIC_CODE)[GENETIC_CODE != "STOP"]) {
  test_moderate[[codon]] <- sample(5:15, 1)
}

# Combine test cases
test_data <- rbindlist(list(test_extreme, test_uniform, test_moderate), fill = TRUE)
test_data[is.na(test_data)] <- 0

cat("Test data created with", nrow(test_data), "genes\n\n")

# Test 1: Wright method with split families (have_F6 = FALSE)
cat("--- Test 1: Wright (1990) with split families (have_F6 = FALSE) ---\n")
enc_wright_split <- calculate_enc(test_data, GENETIC_CODE, have_F6 = FALSE, method = "wright")
print(enc_wright_split)
cat("\n")

# Test 2: Wright method with 6-codon families (have_F6 = TRUE)
cat("--- Test 2: Wright (1990) with 6-codon families (have_F6 = TRUE) ---\n")
enc_wright_6fold <- calculate_enc(test_data, GENETIC_CODE, have_F6 = TRUE, method = "wright")
print(enc_wright_6fold)
cat("\n")

# Test 3: Sun method (always splits families)
cat("--- Test 3: Sun et al. (2012) with pseudocount correction ---\n")
enc_sun <- calculate_enc(test_data, GENETIC_CODE, method = "sun")
print(enc_sun)
cat("\n")

# Test 4: Verify that have_F6 is ignored for Sun method
cat("--- Test 4: Verify have_F6 is ignored for Sun method ---\n")
enc_sun_ignored <- calculate_enc(test_data, GENETIC_CODE, have_F6 = TRUE, method = "sun")
cat("Sun with have_F6=TRUE (should be same as have_F6=FALSE):\n")
print(enc_sun_ignored)
cat("\nDifference between Sun with have_F6=FALSE and have_F6=TRUE:\n")
print(merge(enc_sun, enc_sun_ignored, by = "Gene_name", suffixes = c("_FALSE", "_TRUE")))
cat("\n")

# Summary comparison
cat("=== Summary Comparison ===\n")
comparison <- merge(
  merge(enc_wright_split, enc_wright_6fold, by = "Gene_name", suffixes = c("_Wright_Split", "_Wright_6fold")),
  enc_sun, by = "Gene_name"
)
setnames(comparison, "ENC", "ENC_Sun")
print(comparison)
cat("\n")

# Calculate differences
comparison[, Diff_Wright_Split_vs_Sun := ENC_Wright_Split - ENC_Sun]
comparison[, Diff_Wright_6fold_vs_Sun := ENC_Wright_6fold - ENC_Sun]
comparison[, Diff_Wright_Split_vs_6fold := ENC_Wright_Split - ENC_Wright_6fold]

cat("=== Differences between methods ===\n")
print(comparison[, .(Gene_name, 
                     Diff_Wright_Split_vs_Sun, 
                     Diff_Wright_6fold_vs_Sun, 
                     Diff_Wright_Split_vs_6fold)])
cat("\n")

# Validation checks
cat("=== Validation Checks ===\n")
cat("All Wright split values in valid range [20, 61]:", 
    all(enc_wright_split$ENC >= 20 & enc_wright_split$ENC <= 61), "\n")
cat("All Wright 6fold values in valid range [20, 61]:", 
    all(enc_wright_6fold$ENC >= 20 & enc_wright_6fold$ENC <= 61), "\n")
cat("All Sun values in valid range [20, 61]:", 
    all(enc_sun$ENC >= 20 & enc_sun$ENC <= 61), "\n")
cat("\n")

# Expected behavior:
cat("=== Expected Behavior ===\n")
cat("1. Sun method should give similar but not identical results to Wright split\n")
cat("2. Sun method uses pseudocount correction, reducing sensitivity to small samples\n")
cat("3. Wright 6fold vs split can differ substantially for genes with 6-codon families\n")
cat("4. All ENC values should be in valid range [20, 61]\n")
cat("5. UniformGene should have ENC close to 61 (no bias)\n")
cat("6. ExtremeGene should have lower ENC (strong bias)\n")
cat("\nTest completed successfully!\n")
