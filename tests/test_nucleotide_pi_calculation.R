#!/usr/bin/env Rscript
# Test process_codon_vcf_with_nucleotide_pi function
# Check biallelic pi calculation under various scenarios

library(data.table)

source("./src/derivation_gamma_from_polymorphism.R")
source("./src/local_M_estimation.R")

cat(paste(rep("=", 80), collapse=""), "\n", sep="")
cat("TESTING NUCLEOTIDE PI CALCULATION (BIALLELIC)\n")
cat(paste(rep("=", 80), collapse=""), "\n\n", sep="")

# Setup genetic code
genetic_code_df <- data.frame(
  AA = c(rep("A", 4), rep("V", 4)),
  Codon = c("GCT", "GCC", "GCA", "GCG",  # Alanine
            "GTT", "GTC", "GTA", "GTG")   # Valine
)

# Simple mutation rates (for testing)
aa_mut_rates <- data.table(
  AA = c("A", "V"),
  u = c(0.4, 0.4),
  v = c(0.6, 0.6)
)

# ==============================================================================
# TEST 1: Perfect biallelic split (p = 0.5)
# ==============================================================================
cat("TEST 1: Perfect 50/50 split (p=0.5)\n")
cat(paste(rep("-", 60), collapse=""), "\n", sep="")

vcf_test1 <- data.table(
  Gene = "Gene1",
  Codon_Pos = 1,
  AA = "A",
  Preferred_Codon = "GCT",
  # GCT (preferred): 94 copies, GCC (non-pref): 93 copies
  # Total = 187, p = 94/187 ≈ 0.50
  Codon_Variants = "GCT:94;GCC:93"
)

result1 <- process_codon_vcf_with_nucleotide_pi(vcf_test1, aa_mut_rates, genetic_code_df)
cat("Input: GCT:94 (pref), GCC:93 (non-pref)\n")
cat("Expected p ≈ 0.50\n")
cat("Result:\n")
print(result1[, .(k, n, p, Site_Pi_Codon, Site_Pi_Nucleotide, N_Syn_Positions)])

# Expected: pi_codon = 2*0.5*0.5*(187/186) ≈ 0.5027
# GCT vs GCC differ at position 3 (T vs C), so 1/3 positions differ
# Expected: pi_nucleotide = 0.5027 * (1/3) ≈ 0.168
cat("\nExpected Site_Pi_Codon ≈ 0.50\n")
cat("Expected Site_Pi_Nucleotide ≈ 0.167 (1 position differs out of 3)\n")
cat("Actual Site_Pi_Codon:", result1$Site_Pi_Codon, "\n")
cat("Actual Site_Pi_Nucleotide:", result1$Site_Pi_Nucleotide, "\n")
cat("Status:", ifelse(abs(result1$Site_Pi_Codon - 0.5027) < 0.01 && 
                       result1$Site_Pi_Nucleotide < 0.5, "✓ PASS", "✗ FAIL"), "\n\n")

# ==============================================================================
# TEST 2: Strong bias toward preferred (p ≈ 0.95)
# ==============================================================================
cat("TEST 2: Strong bias toward preferred (p≈0.95)\n")
cat(paste(rep("-", 60), collapse=""), "\n", sep="")

vcf_test2 <- data.table(
  Gene = "Gene2",
  Codon_Pos = 1,
  AA = "A",
  Preferred_Codon = "GCT",
  Codon_Variants = "GCT:178;GCC:9"
)

result2 <- process_codon_vcf_with_nucleotide_pi(vcf_test2, aa_mut_rates, genetic_code_df)
cat("Input: GCT:178 (pref), GCC:9 (non-pref)\n")
cat("Expected p ≈ 0.95\n")
cat("Result:\n")
print(result2[, .(k, n, p, Site_Pi_Codon, Site_Pi_Nucleotide, N_Syn_Positions)])

# Expected: pi_codon = 2*0.95*0.05*(187/186) ≈ 0.0955
# Expected: pi_nucleotide = 0.0955 * (1/3) ≈ 0.032
cat("\nExpected Site_Pi_Codon ≈ 0.096\n")
cat("Expected Site_Pi_Nucleotide ≈ 0.032\n")
cat("Actual Site_Pi_Codon:", result2$Site_Pi_Codon, "\n")
cat("Actual Site_Pi_Nucleotide:", result2$Site_Pi_Nucleotide, "\n")
cat("Status:", ifelse(result2$Site_Pi_Codon < 0.5 && 
                       result2$Site_Pi_Nucleotide < result2$Site_Pi_Codon, "✓ PASS", "✗ FAIL"), "\n\n")

# ==============================================================================
# TEST 3: Monomorphic site (all preferred)
# ==============================================================================
cat("TEST 3: Monomorphic site (p=1.0)\n")
cat(paste(rep("-", 60), collapse=""), "\n", sep="")

vcf_test3 <- data.table(
  Gene = "Gene3",
  Codon_Pos = 1,
  AA = "A",
  Preferred_Codon = "GCT",
  Codon_Variants = "GCT:187"
)

result3 <- process_codon_vcf_with_nucleotide_pi(vcf_test3, aa_mut_rates, genetic_code_df)
cat("Input: GCT:187 (all preferred)\n")
cat("Expected p = 1.0, pi = 0\n")
cat("Result:\n")
print(result3[, .(k, n, p, Site_Pi_Codon, Site_Pi_Nucleotide, N_Syn_Positions)])
cat("Status:", ifelse(result3$Site_Pi_Codon == 0 && 
                       result3$Site_Pi_Nucleotide == 0, "✓ PASS", "✗ FAIL"), "\n\n")

# ==============================================================================
# TEST 4: Multiple non-preferred codons (3 variants)
# ==============================================================================
cat("TEST 4: Multiple non-preferred codons\n")
cat(paste(rep("-", 60), collapse=""), "\n", sep="")

vcf_test4 <- data.table(
  Gene = "Gene4",
  Codon_Pos = 1,
  AA = "A",
  Preferred_Codon = "GCT",
  # GCT (pref): 100, GCC: 40, GCA: 30, GCG: 17
  # Pref = 100, Non-pref = 87, total = 187
  Codon_Variants = "GCT:100;GCC:40;GCA:30;GCG:17"
)

result4 <- process_codon_vcf_with_nucleotide_pi(vcf_test4, aa_mut_rates, genetic_code_df)
cat("Input: GCT:100 (pref), GCC:40, GCA:30, GCG:17 (non-pref)\n")
cat("Expected p = 100/187 ≈ 0.53\n")
cat("Result:\n")
print(result4[, .(k, n, p, Site_Pi_Codon, Site_Pi_Nucleotide, N_Syn_Positions)])

# Expected: pi_codon = 2*0.535*0.465*(187/186) ≈ 0.4995
# GCT vs others: position 3 varies (T vs C/A/G)
# Expected: pi_nucleotide = 0.4995 * (1/3) ≈ 0.167
cat("\nExpected Site_Pi_Codon ≈ 0.50\n")
cat("Expected Site_Pi_Nucleotide ≈ 0.167\n")
cat("Actual Site_Pi_Codon:", result4$Site_Pi_Codon, "\n")
cat("Actual Site_Pi_Nucleotide:", result4$Site_Pi_Nucleotide, "\n")
cat("Status:", ifelse(result4$Site_Pi_Codon < 0.5 && 
                       result4$Site_Pi_Nucleotide < 0.5, "✓ PASS", "✗ FAIL"), "\n\n")

# ==============================================================================
# TEST 5: Valine with 2-position difference
# ==============================================================================
cat("TEST 5: Valine codons (GTT vs GTC - 1 position differs)\n")
cat(paste(rep("-", 60), collapse=""), "\n", sep="")

vcf_test5 <- data.table(
  Gene = "Gene5",
  Codon_Pos = 1,
  AA = "V",
  Preferred_Codon = "GTT",
  Codon_Variants = "GTT:90;GTC:97"
)

result5 <- process_codon_vcf_with_nucleotide_pi(vcf_test5, aa_mut_rates, genetic_code_df)
cat("Input: GTT:90 (pref), GTC:97 (non-pref)\n")
cat("GTT vs GTC differ at position 3 (T vs C)\n")
cat("Result:\n")
print(result5[, .(k, n, p, Site_Pi_Codon, Site_Pi_Nucleotide, N_Syn_Positions)])
cat("Actual Site_Pi_Codon:", result5$Site_Pi_Codon, "\n")
cat("Actual Site_Pi_Nucleotide:", result5$Site_Pi_Nucleotide, "\n")
cat("Status:", ifelse(result5$Site_Pi_Codon < 0.5 && 
                       result5$Site_Pi_Nucleotide < 0.5, "✓ PASS", "✗ FAIL"), "\n\n")

# ==============================================================================
# TEST 6: Edge case - very small sample size
# ==============================================================================
cat("TEST 6: Small sample size (n=2)\n")
cat(paste(rep("-", 60), collapse=""), "\n", sep="")

vcf_test6 <- data.table(
  Gene = "Gene6",
  Codon_Pos = 1,
  AA = "A",
  Preferred_Codon = "GCT",
  Codon_Variants = "GCT:1;GCC:1"
)

result6 <- process_codon_vcf_with_nucleotide_pi(vcf_test6, aa_mut_rates, genetic_code_df)
cat("Input: GCT:1 (pref), GCC:1 (non-pref), n=2\n")
cat("Result:\n")
print(result6[, .(k, n, p, Site_Pi_Codon, Site_Pi_Nucleotide, N_Syn_Positions)])

# Expected: pi = 2*0.5*0.5*(2/1) = 1.0 (maximum for n=2)
cat("\nExpected Site_Pi_Codon = 1.0 (finite sample correction)\n")
cat("Actual Site_Pi_Codon:", result6$Site_Pi_Codon, "\n")
cat("Status:", ifelse(abs(result6$Site_Pi_Codon - 1.0) < 0.01, "✓ PASS", "✗ FAIL"), "\n\n")

# ==============================================================================
# TEST 7: All non-preferred (edge case)
# ==============================================================================
cat("TEST 7: All non-preferred (p=0)\n")
cat(paste(rep("-", 60), collapse=""), "\n", sep="")

vcf_test7 <- data.table(
  Gene = "Gene7",
  Codon_Pos = 1,
  AA = "A",
  Preferred_Codon = "GCT",
  Codon_Variants = "GCC:187"
)

result7 <- process_codon_vcf_with_nucleotide_pi(vcf_test7, aa_mut_rates, genetic_code_df)
cat("Input: GCC:187 (all non-preferred)\n")
cat("Expected p = 0, pi = 0\n")
cat("Result:\n")
print(result7[, .(k, n, p, Site_Pi_Codon, Site_Pi_Nucleotide, N_Syn_Positions)])
cat("Status:", ifelse(result7$Site_Pi_Codon == 0 && 
                       result7$Site_Pi_Nucleotide == 0, "✓ PASS", "✗ FAIL"), "\n\n")

# ==============================================================================
# SUMMARY
# ==============================================================================
cat("\n")
cat(paste(rep("=", 80), collapse=""), "\n", sep="")
cat("SUMMARY\n")
cat(paste(rep("=", 80), collapse=""), "\n\n", sep="")

all_results <- rbind(result1, result2, result3, result4, result5, result6, result7)

cat("Maximum Site_Pi_Codon:", max(all_results$Site_Pi_Codon), "\n")
cat("Maximum Site_Pi_Nucleotide:", max(all_results$Site_Pi_Nucleotide), "\n")
cat("\nChecking constraints:\n")
cat("  All Site_Pi_Codon <= 1.0:", all(all_results$Site_Pi_Codon <= 1.0), 
    ifelse(all(all_results$Site_Pi_Codon <= 1.0), "✓", "✗"), "\n")
cat("  All Site_Pi_Nucleotide <= 0.5:", all(all_results$Site_Pi_Nucleotide <= 0.5), 
    ifelse(all(all_results$Site_Pi_Nucleotide <= 0.5), "✓", "✗"), "\n")
cat("  Site_Pi_Nucleotide <= Site_Pi_Codon:", 
    all(all_results$Site_Pi_Nucleotide <= all_results$Site_Pi_Codon), 
    ifelse(all(all_results$Site_Pi_Nucleotide <= all_results$Site_Pi_Codon), "✓", "✗"), "\n")

cat("\nDistribution of N_Syn_Positions:\n")
print(table(all_results$N_Syn_Positions))

cat("\n")
if (all(all_results$Site_Pi_Codon <= 1.0) && 
    all(all_results$Site_Pi_Nucleotide <= 0.5) &&
    all(all_results$Site_Pi_Nucleotide <= all_results$Site_Pi_Codon)) {
  cat("✓ ALL TESTS PASSED - Function behaves correctly!\n")
} else {
  cat("✗ SOME TESTS FAILED - Review edge cases\n")
}

cat("\n")
