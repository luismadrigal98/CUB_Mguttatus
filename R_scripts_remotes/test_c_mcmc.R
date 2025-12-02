#!/usr/bin/env Rscript
#
# Quick test of C-accelerated MCMC
#

cat("============================================\n")
cat("  MCMC Speed Test: Pure R vs C-accelerated\n")
cat("============================================\n\n")

script_dir <- "/mnt/1692B2EF92B2D28B/Ongoing_projects/Codon_bias_analysis/R_scripts_remotes"

# Source the main file (which loads C if available)
old_args <- commandArgs
commandArgs <- function(...) c()
source(file.path(script_dir, "ROC_model_pure_R.R"))
commandArgs <- old_args

# Generate small test dataset
set.seed(42)
n_genes <- 500
aa_to_codons <- get_genetic_code()
syn_aas <- get_synonymous_aas()
all_codons <- unlist(aa_to_codons)

codon_counts <- matrix(0L, nrow = n_genes, ncol = length(all_codons))
colnames(codon_counts) <- all_codons
rownames(codon_counts) <- paste0("Gene_", seq_len(n_genes))

for (g in seq_len(n_genes)) {
  for (aa in syn_aas) {
    codons <- aa_to_codons[[aa]]
    n_aa <- rpois(1, 10)
    if (n_aa > 0) {
      counts <- rmultinom(1, n_aa, rep(1, length(codons)))
      codon_counts[g, codons] <- as.integer(counts)
    }
  }
}

cat(sprintf("Test dataset: %d genes\n\n", n_genes))

# Test 1: With C acceleration (if available)
cat("--- Test with C acceleration ---\n")
start <- Sys.time()
result1 <- run_roc_mcmc(
  codon_counts = codon_counts,
  n_samples = 50,
  thin = 2,
  verbose = TRUE
)
time1 <- as.numeric(Sys.time() - start)
cat(sprintf("Time: %.2f seconds\n\n", time1))

# Test 2: Force pure R by temporarily disabling C
cat("--- Test with pure R (C disabled) ---\n")
old_use_c <- .USE_C_CODE
.USE_C_CODE <<- FALSE

start <- Sys.time()
result2 <- run_roc_mcmc(
  codon_counts = codon_counts,
  n_samples = 50,
  thin = 2,
  n_cores = 4,  # Use parallelization
  verbose = TRUE
)
time2 <- as.numeric(Sys.time() - start)
cat(sprintf("Time: %.2f seconds\n\n", time2))

# Restore
.USE_C_CODE <<- old_use_c

# Summary
cat("============================================\n")
cat("  Summary\n")
cat("============================================\n")
cat(sprintf("C-accelerated: %.2f sec\n", time1))
cat(sprintf("Pure R (4 cores): %.2f sec\n", time2))
cat(sprintf("Speedup: %.1fx\n", time2 / time1))

# Verify results are similar
phi_cor <- cor(result1$phi_mean, result2$phi_mean)
cat(sprintf("\nPhi correlation between methods: %.4f\n", phi_cor))
if (phi_cor > 0.99) {
  cat("✓ Results are consistent\n")
} else {
  cat("⚠ Results differ more than expected\n")
}

cat("============================================\n")
