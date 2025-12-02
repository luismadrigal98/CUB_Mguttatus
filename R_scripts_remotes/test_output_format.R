#!/usr/bin/env Rscript
#
# Test the output format of ROC_model_pure_R.R
#

cat("============================================\n")
cat("  Testing Output Format\n")
cat("============================================\n\n")

script_dir <- "/mnt/1692B2EF92B2D28B/Ongoing_projects/Codon_bias_analysis/R_scripts_remotes"
output_dir <- tempfile("roc_test_output_")

# Source the main file
old_args <- commandArgs
commandArgs <- function(...) c()
source(file.path(script_dir, "ROC_model_pure_R.R"))
commandArgs <- old_args

# Generate small test dataset
set.seed(42)
n_genes <- 100
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

# Create mock multi-tissue expression
obs_phi_matrix <- matrix(rlnorm(n_genes * 2, -0.5, 1), ncol = 2)
colnames(obs_phi_matrix) <- c("Leaf", "Bud")
rownames(obs_phi_matrix) <- rownames(codon_counts)
# Normalize
obs_phi_matrix[, 1] <- obs_phi_matrix[, 1] / mean(obs_phi_matrix[, 1])
obs_phi_matrix[, 2] <- obs_phi_matrix[, 2] / mean(obs_phi_matrix[, 2])

cat(sprintf("Test dataset: %d genes, 2 tissues\n", n_genes))
cat(sprintf("Output directory: %s\n\n", output_dir))

# Run short MCMC
cat("Running MCMC...\n")
result <- run_roc_mcmc(
  codon_counts = codon_counts,
  n_samples = 30,
  thin = 2,
  obs_phi_matrix = obs_phi_matrix,
  with_phi = TRUE,
  verbose = TRUE
)

# Save results
cat("\nSaving results...\n")
save_results(result, output_dir, create_plots = TRUE)

# Check output files
cat("\n============================================\n")
cat("  Output Files Created\n")
cat("============================================\n")

list_files_recursive <- function(path, prefix = "") {
  files <- list.files(path, full.names = TRUE)
  for (f in files) {
    if (file.info(f)$isdir) {
      cat(sprintf("%s%s/\n", prefix, basename(f)))
      list_files_recursive(f, paste0(prefix, "  "))
    } else {
      size <- file.info(f)$size
      cat(sprintf("%s%s (%.1f KB)\n", prefix, basename(f), size / 1024))
    }
  }
}

list_files_recursive(output_dir)

# Show sample content
cat("\n============================================\n")
cat("  Sample Output Content\n")
cat("============================================\n")

cat("\n--- gene_expression.csv (first 5 rows) ---\n")
print(head(read.csv(file.path(output_dir, "Parameter_est", "gene_expression.csv")), 5))

cat("\n--- selection_coefficients.csv (first 10 rows) ---\n")
print(head(read.csv(file.path(output_dir, "Parameter_est", "selection_coefficients.csv")), 10))

cat("\n--- hyperparameters.csv ---\n")
print(read.csv(file.path(output_dir, "Parameter_est", "hyperparameters.csv")))

cat("\n--- convergence_diagnostics.csv ---\n")
print(read.csv(file.path(output_dir, "Parameter_est", "convergence_diagnostics.csv")))

cat("\n--- summary.txt ---\n")
cat(readLines(file.path(output_dir, "summary.txt")), sep = "\n")

# Cleanup
cat("\n\nCleaning up temporary directory...\n")
unlink(output_dir, recursive = TRUE)

cat("\n✓ Test complete!\n")
