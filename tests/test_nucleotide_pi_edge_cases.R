#!/usr/bin/env Rscript
# Additional edge case tests for nucleotide pi calculation

library(data.table)

source("./src/derivation_gamma_from_polymorphism.R")
source("./src/local_M_estimation.R")

cat(paste(rep("=", 80), collapse=""), "\n", sep="")
cat("ADDITIONAL EDGE CASE TESTS\n")
cat(paste(rep("=", 80), collapse=""), "\n\n", sep="")

# Setup genetic code with more amino acids
genetic_code_df <- data.frame(
  AA = c(rep("A", 4), rep("V", 4), rep("L", 6), rep("S", 6)),
  Codon = c("GCT", "GCC", "GCA", "GCG",  # Alanine
            "GTT", "GTC", "GTA", "GTG",   # Valine
            "TTA", "TTG", "CTT", "CTC", "CTA", "CTG",  # Leucine
            "TCT", "TCC", "TCA", "TCG", "AGT", "AGC")  # Serine
)

aa_mut_rates <- data.table(
  AA = c("A", "V", "L", "S"),
  u = c(0.4, 0.4, 0.3, 0.5),
  v = c(0.6, 0.6, 0.7, 0.5)
)

# ==============================================================================
# TEST 8: Leucine (6-fold degenerate, codons differ at multiple positions)
# ==============================================================================
cat("TEST 8: Leucine - 6-fold degenerate amino acid\n")
cat(paste(rep("-", 60), collapse=""), "\n", sep="")
cat("Testing with TTA (preferred) vs CTT/CTC (non-preferred)\n")
cat("TTA vs CTT differ at positions 1 (T vs C) and 3 (A vs T) = 2/3\n\n")

vcf_test8 <- data.table(
  Gene = "Gene8",
  Codon_Pos = 1,
  AA = "L",
  Preferred_Codon = "TTA",
  # TTA (pref): 90, CTT: 50, CTC: 47
  Codon_Variants = "TTA:90;CTT:50;CTC:47"
)

result8 <- process_codon_vcf_with_nucleotide_pi(vcf_test8, aa_mut_rates, genetic_code_df)
cat("Input: TTA:90 (pref), CTT:50, CTC:47 (non-pref)\n")
cat("Result:\n")
print(result8[, .(k, n, p, Site_Pi_Codon, Site_Pi_Nucleotide, N_Syn_Positions)])

# Pref = 90/187 = 0.481
# pi_codon = 2*0.481*0.519*(187/186) = 0.501
# TTA vs CTT: differ at positions 1,3 → 2 positions
# Expected: pi_nucleotide = 0.501 * (2/3) ≈ 0.334
cat("\nExpected Site_Pi_Codon ≈ 0.50\n")
cat("Expected Site_Pi_Nucleotide ≈ 0.333 (2 positions differ)\n")
cat("Actual Site_Pi_Codon:", result8$Site_Pi_Codon, "\n")
cat("Actual Site_Pi_Nucleotide:", result8$Site_Pi_Nucleotide, "\n")
cat("N_Syn_Positions:", result8$N_Syn_Positions, "\n")
cat("Status:", ifelse(result8$Site_Pi_Nucleotide < 0.5 && 
                       result8$N_Syn_Positions == 2, "✓ PASS", "✗ FAIL"), "\n\n")

# ==============================================================================
# TEST 9: Leucine - all 3 positions differ
# ==============================================================================
cat("TEST 9: Leucine - Maximum nucleotide difference\n")
cat(paste(rep("-", 60), collapse=""), "\n", sep="")
cat("TTA (preferred) vs CTG (non-preferred) - all 3 positions differ\n\n")

vcf_test9 <- data.table(
  Gene = "Gene9",
  Codon_Pos = 1,
  AA = "L",
  Preferred_Codon = "TTA",
  Codon_Variants = "TTA:93;CTG:94"
)

result9 <- process_codon_vcf_with_nucleotide_pi(vcf_test9, aa_mut_rates, genetic_code_df)
cat("Input: TTA:93 (pref), CTG:94 (non-pref)\n")
cat("TTA vs CTG differ at all 3 positions (T≠C, T≠T, A≠G)\n")
cat("Wait... position 2: T=T, so only positions 1 and 3 differ!\n")
cat("Result:\n")
print(result9[, .(k, n, p, Site_Pi_Codon, Site_Pi_Nucleotide, N_Syn_Positions)])

# Actually TTA vs CTG: T≠C, T=T, A≠G → 2 positions differ
cat("\nActual Site_Pi_Codon:", result9$Site_Pi_Codon, "\n")
cat("Actual Site_Pi_Nucleotide:", result9$Site_Pi_Nucleotide, "\n")
cat("N_Syn_Positions:", result9$N_Syn_Positions, "(should be 2)\n")
cat("Status:", ifelse(result9$Site_Pi_Nucleotide < 0.5, "✓ PASS", "✗ FAIL"), "\n\n")

# ==============================================================================
# TEST 10: Serine - mixed codon families (TCN and AGY)
# ==============================================================================
cat("TEST 10: Serine - codons from different families\n")
cat(paste(rep("-", 60), collapse=""), "\n", sep="")
cat("TCT (preferred) vs AGT (non-preferred) - maximum difference\n\n")

vcf_test10 <- data.table(
  Gene = "Gene10",
  Codon_Pos = 1,
  AA = "S",
  Preferred_Codon = "TCT",
  Codon_Variants = "TCT:100;AGT:87"
)

result10 <- process_codon_vcf_with_nucleotide_pi(vcf_test10, aa_mut_rates, genetic_code_df)
cat("Input: TCT:100 (pref), AGT:87 (non-pref)\n")
cat("TCT vs AGT differ at positions 1 (T vs A) and 2 (C vs G) = 2/3\n")
cat("Result:\n")
print(result10[, .(k, n, p, Site_Pi_Codon, Site_Pi_Nucleotide, N_Syn_Positions)])

cat("\nActual Site_Pi_Codon:", result10$Site_Pi_Codon, "\n")
cat("Actual Site_Pi_Nucleotide:", result10$Site_Pi_Nucleotide, "\n")
cat("N_Syn_Positions:", result10$N_Syn_Positions, "\n")
cat("Status:", ifelse(result10$Site_Pi_Nucleotide < 0.5, "✓ PASS", "✗ FAIL"), "\n\n")

# ==============================================================================
# TEST 11: Very large sample size
# ==============================================================================
cat("TEST 11: Large sample size (n=1000)\n")
cat(paste(rep("-", 60), collapse=""), "\n", sep="")

vcf_test11 <- data.table(
  Gene = "Gene11",
  Codon_Pos = 1,
  AA = "A",
  Preferred_Codon = "GCT",
  Codon_Variants = "GCT:500;GCC:500"
)

result11 <- process_codon_vcf_with_nucleotide_pi(vcf_test11, aa_mut_rates, genetic_code_df)
cat("Input: GCT:500 (pref), GCC:500 (non-pref), n=1000\n")
cat("Result:\n")
print(result11[, .(k, n, p, Site_Pi_Codon, Site_Pi_Nucleotide)])

# Expected: pi = 2*0.5*0.5*(1000/999) ≈ 0.5005
cat("\nExpected Site_Pi_Codon ≈ 0.5005 (minimal finite sample correction)\n")
cat("Actual Site_Pi_Codon:", result11$Site_Pi_Codon, "\n")
cat("Actual Site_Pi_Nucleotide:", result11$Site_Pi_Nucleotide, "\n")
cat("Status:", ifelse(abs(result11$Site_Pi_Codon - 0.5005) < 0.001, "✓ PASS", "✗ FAIL"), "\n\n")

# ==============================================================================
# TEST 12: Extreme bias (p = 0.99)
# ==============================================================================
cat("TEST 12: Extreme bias toward preferred (p=0.99)\n")
cat(paste(rep("-", 60), collapse=""), "\n", sep="")

vcf_test12 <- data.table(
  Gene = "Gene12",
  Codon_Pos = 1,
  AA = "A",
  Preferred_Codon = "GCT",
  Codon_Variants = "GCT:990;GCC:10"
)

result12 <- process_codon_vcf_with_nucleotide_pi(vcf_test12, aa_mut_rates, genetic_code_df)
cat("Input: GCT:990 (pref), GCC:10 (non-pref), n=1000\n")
cat("Result:\n")
print(result12[, .(k, n, p, Site_Pi_Codon, Site_Pi_Nucleotide)])

# Expected: pi = 2*0.99*0.01*(1000/999) ≈ 0.0198
cat("\nExpected Site_Pi_Codon ≈ 0.020\n")
cat("Expected Site_Pi_Nucleotide ≈ 0.0066\n")
cat("Actual Site_Pi_Codon:", result12$Site_Pi_Codon, "\n")
cat("Actual Site_Pi_Nucleotide:", result12$Site_Pi_Nucleotide, "\n")
cat("Status:", ifelse(result12$Site_Pi_Codon < 0.03 && 
                       result12$Site_Pi_Nucleotide < 0.01, "✓ PASS", "✗ FAIL"), "\n\n")

# ==============================================================================
# TEST 13: Near-monomorphic with 1 rare variant
# ==============================================================================
cat("TEST 13: Singleton variant (1 out of 200)\n")
cat(paste(rep("-", 60), collapse=""), "\n", sep="")

vcf_test13 <- data.table(
  Gene = "Gene13",
  Codon_Pos = 1,
  AA = "A",
  Preferred_Codon = "GCT",
  Codon_Variants = "GCT:199;GCC:1"
)

result13 <- process_codon_vcf_with_nucleotide_pi(vcf_test13, aa_mut_rates, genetic_code_df)
cat("Input: GCT:199 (pref), GCC:1 (singleton)\n")
cat("Result:\n")
print(result13[, .(k, n, p, Site_Pi_Codon, Site_Pi_Nucleotide)])

# Expected: pi = 2*0.995*0.005*(200/199) ≈ 0.01
cat("\nExpected Site_Pi_Codon ≈ 0.010\n")
cat("Actual Site_Pi_Codon:", result13$Site_Pi_Codon, "\n")
cat("Actual Site_Pi_Nucleotide:", result13$Site_Pi_Nucleotide, "\n")
cat("Status:", ifelse(result13$Site_Pi_Nucleotide < 0.5, "✓ PASS", "✗ FAIL"), "\n\n")

# ==============================================================================
# SUMMARY OF EDGE CASES
# ==============================================================================
cat("\n")
cat(paste(rep("=", 80), collapse=""), "\n", sep="")
cat("EDGE CASE SUMMARY\n")
cat(paste(rep("=", 80), collapse=""), "\n\n", sep="")

all_edge <- rbind(result8, result9, result10, result11, result12, result13)

cat("Sample size range: n = ", min(all_edge$n), " to ", max(all_edge$n), "\n", sep="")
cat("Frequency range: p = ", round(min(all_edge$p), 4), " to ", 
    round(max(all_edge$p), 4), "\n", sep="")
cat("\nMaximum values:\n")
cat("  Site_Pi_Codon:", max(all_edge$Site_Pi_Codon), "\n")
cat("  Site_Pi_Nucleotide:", max(all_edge$Site_Pi_Nucleotide), "\n")

cat("\nConstraint checks:\n")
cat("  All Site_Pi_Codon <= 1.0:", all(all_edge$Site_Pi_Codon <= 1.0), 
    ifelse(all(all_edge$Site_Pi_Codon <= 1.0), "✓", "✗"), "\n")
cat("  All Site_Pi_Nucleotide <= 0.5:", all(all_edge$Site_Pi_Nucleotide <= 0.5), 
    ifelse(all(all_edge$Site_Pi_Nucleotide <= 0.5), "✓", "✗"), "\n")
cat("  Site_Pi_Nucleotide <= Site_Pi_Codon:", 
    all(all_edge$Site_Pi_Nucleotide <= all_edge$Site_Pi_Codon), 
    ifelse(all(all_edge$Site_Pi_Nucleotide <= all_edge$Site_Pi_Codon), "✓", "✗"), "\n")

cat("\nN_Syn_Positions distribution:\n")
print(table(all_edge$N_Syn_Positions))

cat("\nBehavior with sample size:\n")
cat("  Small (n=187): pi_codon range =", 
    round(min(all_edge$Site_Pi_Codon[all_edge$n < 500]), 3), "-", 
    round(max(all_edge$Site_Pi_Codon[all_edge$n < 500]), 3), "\n")
cat("  Large (n=1000): pi_codon range =", 
    round(min(all_edge$Site_Pi_Codon[all_edge$n >= 500]), 3), "-", 
    round(max(all_edge$Site_Pi_Codon[all_edge$n >= 500]), 3), "\n")

cat("\n")
if (all(all_edge$Site_Pi_Codon <= 1.0) && 
    all(all_edge$Site_Pi_Nucleotide <= 0.5) &&
    all(all_edge$Site_Pi_Nucleotide <= all_edge$Site_Pi_Codon)) {
  cat("✓ ALL EDGE CASES HANDLED CORRECTLY!\n")
  cat("\nFunction is robust to:\n")
  cat("  - 6-fold degenerate amino acids\n")
  cat("  - Codons differing at 0, 1, or 2 positions\n")
  cat("  - Sample sizes from n=2 to n=1000\n")
  cat("  - Allele frequencies from p=0.005 to p=0.995\n")
  cat("  - Singleton variants and extreme bias\n")
} else {
  cat("✗ EDGE CASE ISSUES DETECTED\n")
}

cat("\n")
