#!/usr/bin/env Rscript
#
# Benchmark Pure R vs C-accelerated ROC likelihood
#
# This script compares the performance of the pure R implementation
# against the C-accelerated version.

cat("============================================\n")
cat("  ROC Likelihood: R vs C Benchmark\n")
cat("============================================\n\n")

# Get script directory
script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg)))
  } else {
    "/mnt/1692B2EF92B2D28B/Ongoing_projects/Codon_bias_analysis/R_scripts_remotes"
  }
}, error = function(e) {
  "/mnt/1692B2EF92B2D28B/Ongoing_projects/Codon_bias_analysis/R_scripts_remotes"
})

# Source the pure R implementation (without CLI)
cat("Loading pure R implementation...\n")
old_args <- commandArgs
commandArgs <- function(...) c()  # Suppress argparse

source(file.path(script_dir, "ROC_model_pure_R.R"))

commandArgs <- old_args  # Restore

# Source C wrapper
cat("Loading C-accelerated functions...\n")
source(file.path(script_dir, "roc_likelihood_c.R"))

cat(sprintf("\nC acceleration available: %s\n\n", is_c_available()))

# ==============================================================================
# Generate test data
# ==============================================================================

cat("Generating test data...\n")

# Parameters
n_genes <- 1000
n_codons_per_gene <- 200  # Average gene length in codons

# Get genetic code
aa_to_codons <- get_genetic_code()
syn_aas <- get_synonymous_aas()
all_codons <- unlist(aa_to_codons)

# Generate random codon counts
set.seed(42)
codon_counts <- matrix(0L, nrow = n_genes, ncol = length(all_codons))
colnames(codon_counts) <- all_codons
rownames(codon_counts) <- paste0("Gene_", seq_len(n_genes))

for (g in seq_len(n_genes)) {
  # Random amino acid composition
  for (aa in syn_aas) {
    codons <- aa_to_codons[[aa]]
    # Random total for this AA
    n_aa <- rpois(1, n_codons_per_gene / length(syn_aas))
    if (n_aa > 0) {
      # Distribute among codons
      counts <- rmultinom(1, n_aa, rep(1, length(codons)))
      codon_counts[g, codons] <- as.integer(counts)
    }
  }
}

# Generate random parameters
dM_list <- list()
dEta_list <- list()
for (aa in syn_aas) {
  n_params <- length(aa_to_codons[[aa]]) - 1
  dM_list[[aa]] <- rnorm(n_params, 0, 0.3)
  dEta_list[[aa]] <- rnorm(n_params, 0, 0.3)
}

# Generate random phi
phi <- rlnorm(n_genes, -0.5, 1)
phi <- phi / mean(phi)

cat(sprintf("  Genes: %d\n", n_genes))
cat(sprintf("  Total codons: %d\n", sum(codon_counts)))
cat(sprintf("  Amino acids with synonymous codons: %d\n", length(syn_aas)))

# ==============================================================================
# Benchmark: Total log likelihood
# ==============================================================================

cat("\n--- Benchmark: Total Log Likelihood ---\n")

n_reps <- 10

# Pure R (serial)
cat("Running pure R (1 core)...\n")
r_times <- numeric(n_reps)
for (i in seq_len(n_reps)) {
  start <- Sys.time()
  ll_r <- calc_total_log_likelihood(codon_counts, dM_list, dEta_list, 
                                     phi, aa_to_codons, n_cores = 1)
  r_times[i] <- as.numeric(Sys.time() - start)
}
cat(sprintf("  Mean time: %.4f sec (SD: %.4f)\n", mean(r_times), sd(r_times)))
cat(sprintf("  Log likelihood: %.2f\n", ll_r))

# Pure R (parallel)
n_cores <- parallel::detectCores() - 1
cat(sprintf("Running pure R (%d cores)...\n", n_cores))
r_par_times <- numeric(n_reps)
for (i in seq_len(n_reps)) {
  start <- Sys.time()
  ll_r_par <- calc_total_log_likelihood(codon_counts, dM_list, dEta_list,
                                         phi, aa_to_codons, n_cores = n_cores)
  r_par_times[i] <- as.numeric(Sys.time() - start)
}
cat(sprintf("  Mean time: %.4f sec (SD: %.4f)\n", mean(r_par_times), sd(r_par_times)))

# C version
if (is_c_available()) {
  cat("Running C version...\n")
  c_times <- numeric(n_reps)
  for (i in seq_len(n_reps)) {
    start <- Sys.time()
    ll_c <- calc_total_log_likelihood_c(codon_counts, dM_list, dEta_list,
                                         phi, aa_to_codons)
    c_times[i] <- as.numeric(Sys.time() - start)
  }
  cat(sprintf("  Mean time: %.4f sec (SD: %.4f)\n", mean(c_times), sd(c_times)))
  cat(sprintf("  Log likelihood: %.2f\n", ll_c))
  
  # Verify correctness
  if (abs(ll_r - ll_c) < 1e-6) {
    cat("  ✓ Results match pure R\n")
  } else {
    cat(sprintf("  ✗ MISMATCH: R=%.6f, C=%.6f, diff=%.6f\n", ll_r, ll_c, ll_r - ll_c))
  }
  
  cat(sprintf("\n  SPEEDUP: %.1fx vs R(serial), %.1fx vs R(parallel)\n",
              mean(r_times) / mean(c_times),
              mean(r_par_times) / mean(c_times)))
}

# ==============================================================================
# Benchmark: Single gene likelihood
# ==============================================================================

cat("\n--- Benchmark: Single Gene Likelihood (x1000 calls) ---\n")

n_gene_reps <- 1000

# Pure R
cat("Running pure R...\n")
start <- Sys.time()
for (i in seq_len(n_gene_reps)) {
  g <- ((i - 1) %% n_genes) + 1
  ll <- calc_log_likelihood_gene(codon_counts[g, ], dM_list, dEta_list,
                                  phi[g], aa_to_codons)
}
r_gene_time <- as.numeric(Sys.time() - start)
cat(sprintf("  Time: %.4f sec (%.4f ms/call)\n", r_gene_time, r_gene_time * 1000 / n_gene_reps))

# C version
if (is_c_available()) {
  cat("Running C version...\n")
  start <- Sys.time()
  for (i in seq_len(n_gene_reps)) {
    g <- ((i - 1) %% n_genes) + 1
    ll <- calc_log_likelihood_gene_c(codon_counts[g, ], dM_list, dEta_list,
                                      phi[g], aa_to_codons)
  }
  c_gene_time <- as.numeric(Sys.time() - start)
  cat(sprintf("  Time: %.4f sec (%.4f ms/call)\n", c_gene_time, c_gene_time * 1000 / n_gene_reps))
  cat(sprintf("\n  SPEEDUP: %.1fx\n", r_gene_time / c_gene_time))
}

# ==============================================================================
# Benchmark: Batch phi update (most critical for MCMC)
# ==============================================================================

cat("\n--- Benchmark: Batch Phi Update (100 genes) ---\n")

n_batch_reps <- 10
batch_size <- 100
gene_indices <- seq_len(batch_size)
prop_sd <- rep(0.3, n_genes)

# Pure R fallback
cat("Running pure R...\n")
r_batch_times <- numeric(n_batch_reps)
for (i in seq_len(n_batch_reps)) {
  start <- Sys.time()
  result_r <- .batch_update_phi_r(codon_counts, dM_list, dEta_list, phi,
                                   NULL, 0.5, 1.0, prop_sd,
                                   aa_to_codons, FALSE, gene_indices)
  r_batch_times[i] <- as.numeric(Sys.time() - start)
}
cat(sprintf("  Mean time: %.4f sec\n", mean(r_batch_times)))

# C version
if (is_c_available()) {
  cat("Running C version...\n")
  c_batch_times <- numeric(n_batch_reps)
  for (i in seq_len(n_batch_reps)) {
    start <- Sys.time()
    result_c <- batch_update_phi_c(codon_counts, dM_list, dEta_list, phi,
                                    NULL, 0.5, 1.0, prop_sd,
                                    aa_to_codons, FALSE, gene_indices)
    c_batch_times[i] <- as.numeric(Sys.time() - start)
  }
  cat(sprintf("  Mean time: %.4f sec\n", mean(c_batch_times)))
  cat(sprintf("\n  SPEEDUP: %.1fx\n", mean(r_batch_times) / mean(c_batch_times)))
}

# ==============================================================================
# Summary
# ==============================================================================

cat("\n============================================\n")
cat("  Summary\n")
cat("============================================\n")

if (is_c_available()) {
  total_speedup <- mean(r_times) / mean(c_times)
  cat(sprintf("Total likelihood speedup: %.1fx\n", total_speedup))
  cat(sprintf("Single gene speedup: %.1fx\n", r_gene_time / c_gene_time))
  cat(sprintf("Batch update speedup: %.1fx\n", mean(r_batch_times) / mean(c_batch_times)))
  
  # Estimate MCMC speedup
  # In MCMC, most time is in:
  # - phi updates (n_genes per iteration): ~70% of time
  # - dM/dEta updates (recalc all genes): ~30% of time
  est_mcmc_speedup <- 0.7 * (mean(r_batch_times) / mean(c_batch_times)) + 
                       0.3 * total_speedup
  cat(sprintf("\nEstimated MCMC speedup: ~%.1fx\n", est_mcmc_speedup))
} else {
  cat("C acceleration not available.\n")
  cat("Run compile_roc_c() to compile and enable C acceleration.\n")
}

cat("============================================\n")
