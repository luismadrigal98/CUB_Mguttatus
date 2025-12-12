#!/usr/bin/env Rscript
# Quick diagnostic to check if your data has reasonable pi distribution

library(data.table)

source("./src/derivation_gamma_from_polymorphism.R")
source("./src/local_M_estimation.R")

cat(paste(rep("=", 80), collapse=""), "\n")
cat("DIAGNOSTIC: Check your data's pi distribution\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

# Setup genetic code
genetic_code_df <- data.frame(
  AA = c(rep("A", 4), rep("V", 4)),
  Codon = c("GCT", "GCC", "GCA", "GCG",  # Alanine
            "GTT", "GTC", "GTA", "GTG")   # Valine
)

aa_mut_rates <- data.table(
  AA = c("A", "V"),
  u = c(0.4, 0.4),
  v = c(0.6, 0.6)
)

# Simulate realistic data distribution matching theta = 0.0312
set.seed(123)
n_sites <- 1000
theta <- 0.0312

cat("Simulating", n_sites, "sites with theta =", theta, "\n\n")

# Generate sites with realistic SFS
sim_data <- lapply(1:n_sites, function(i) {
  # Most sites monomorphic
  if (runif(1) > theta * 5) {  # ~85% monomorphic
    data.table(
      Gene = paste0("Gene", (i-1) %/% 10 + 1),
      Codon_Pos = (i-1) %% 10 + 1,
      AA = "A",
      Preferred_Codon = "GCT",
      Codon_Variants = "GCT:187"
    )
  } else {
    # Polymorphic site - use U-shaped SFS
    # Most variants are rare or near-fixation
    p_pref <- rbeta(1, 0.5, 0.5)  # U-shaped
    n_pref <- round(p_pref * 187)
    n_nonpref <- 187 - n_pref
    
    data.table(
      Gene = paste0("Gene", (i-1) %/% 10 + 1),
      Codon_Pos = (i-1) %% 10 + 1,
      AA = "A",
      Preferred_Codon = "GCT",
      Codon_Variants = sprintf("GCT:%d;GCC:%d", n_pref, n_nonpref)
    )
  }
})

sim_vcf <- rbindlist(sim_data)

# Process with the function
results <- process_codon_vcf_with_nucleotide_pi(sim_vcf, aa_mut_rates, genetic_code_df)

cat("SUMMARY STATISTICS\n")
cat(paste(rep("-", 60), collapse=""), "\n")
cat("Total sites:", nrow(results), "\n")
cat("Monomorphic sites:", sum(results$Site_Pi_Nucleotide == 0), 
    sprintf("(%.1f%%)\n", 100 * mean(results$Site_Pi_Nucleotide == 0)))
cat("Polymorphic sites:", sum(results$Site_Pi_Nucleotide > 0), 
    sprintf("(%.1f%%)\n", 100 * mean(results$Site_Pi_Nucleotide > 0)))

cat("\nNUCLEOTIDE PI DISTRIBUTION\n")
cat(paste(rep("-", 60), collapse=""), "\n")
cat("Mean Site_Pi_Nucleotide:", round(mean(results$Site_Pi_Nucleotide), 4), "\n")
cat("Median Site_Pi_Nucleotide:", round(median(results$Site_Pi_Nucleotide), 4), "\n")
cat("Max Site_Pi_Nucleotide:", round(max(results$Site_Pi_Nucleotide), 4), "\n")
cat("SD Site_Pi_Nucleotide:", round(sd(results$Site_Pi_Nucleotide), 4), "\n")

cat("\nQUANTILES\n")
cat(paste(rep("-", 60), collapse=""), "\n")
quantiles <- quantile(results$Site_Pi_Nucleotide, probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1))
print(round(quantiles, 4))

cat("\nDISTRIBUTION BY RANGE\n")
cat(paste(rep("-", 60), collapse=""), "\n")
ranges <- cut(results$Site_Pi_Nucleotide, 
              breaks = c(-0.01, 0, 0.01, 0.05, 0.1, 0.2, 0.5, 1),
              labels = c("0 (mono)", "(0, 0.01]", "(0.01, 0.05]", 
                        "(0.05, 0.1]", "(0.1, 0.2]", "(0.2, 0.5]", ">0.5"))
table_ranges <- table(ranges)
print(table_ranges)
cat("\nProportions:\n")
print(round(prop.table(table_ranges), 3))

cat("\nSITES WITH HIGH PI (> 0.05)\n")
cat(paste(rep("-", 60), collapse=""), "\n")
high_pi <- results[Site_Pi_Nucleotide > 0.05]
cat("Count:", nrow(high_pi), sprintf("(%.1f%%)\n", 100 * nrow(high_pi) / nrow(results)))
if (nrow(high_pi) > 0) {
  cat("Range:", round(min(high_pi$Site_Pi_Nucleotide), 4), "to", 
      round(max(high_pi$Site_Pi_Nucleotide), 4), "\n")
  cat("Mean frequency (p):", round(mean(high_pi$p), 3), "\n")
  cat("\nFirst few examples:\n")
  print(head(high_pi[, .(Gene, Codon_Pos, k, n, p, Site_Pi_Codon, Site_Pi_Nucleotide)], 5))
}

cat("\n")
cat(paste(rep("=", 80), collapse=""), "\n")
cat("INTERPRETATION\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

mean_pi <- mean(results$Site_Pi_Nucleotide)
pct_high_pi <- 100 * mean(results$Site_Pi_Nucleotide > 0.05)
max_pi <- max(results$Site_Pi_Nucleotide)

if (max_pi > 0.5) {
  cat("✗ ERROR: Maximum pi > 0.5 detected! This violates biallelic constraint.\n")
} else if (mean_pi > 0.1) {
  cat("⚠ WARNING: Mean pi very high (", round(mean_pi, 4), "). Check for:\n")
  cat("  - Data quality issues\n")
  cat("  - Non-synonymous variants in data\n")
  cat("  - Incorrect theta estimate\n")
} else if (pct_high_pi > 20) {
  cat("⚠ NOTICE: ", round(pct_high_pi, 1), "% of sites have pi > 0.05\n")
  cat("  This is higher than expected with theta = 0.0312\n")
  cat("  But individual high-pi sites are biologically possible\n")
  cat("  Consider: Are these real balanced polymorphisms?\n")
} else {
  cat("✓ HEALTHY DISTRIBUTION\n")
  cat("  - Mean pi =", round(mean_pi, 4), "(expected ~0.01-0.03)\n")
  cat("  - Max pi =", round(max_pi, 4), "(<0.5 biallelic limit)\n")
  cat("  - ", round(pct_high_pi, 1), "% sites with pi > 0.05 (acceptable)\n")
  cat("\n  High-pi sites (pi > 0.05) represent balanced polymorphisms.\n")
  cat("  These are biologically real and expected in ~1-10% of sites.\n")
}

cat("\n")
