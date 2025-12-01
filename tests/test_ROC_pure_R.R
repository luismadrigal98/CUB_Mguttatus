#!/usr/bin/env Rscript
#
# Test script for the Pure R ROC Model implementation
# Uses the actual Mguttatus data files
#
# Usage:
#   Rscript test_ROC_pure_R.R
#
# Or in an interactive R session:
#   source("test_ROC_pure_R.R")
#

# Source the main implementation
# Get script directory robustly (works with Rscript and source())
get_script_dir <- function() {
  # Try commandArgs first (Rscript)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  # Try sys.frame (source())
  for (i in seq_len(sys.nframe())) {
    ofile <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(ofile)) {
      return(dirname(normalizePath(ofile)))
    }
  }
  # Fallback to working directory
  return(getwd())
}

script_dir <- get_script_dir()
message(sprintf("Script directory: %s", script_dir))

# ROC_model_pure_R.R is in R_scripts_remotes/, test is in tests/
roc_script <- file.path(dirname(script_dir), "R_scripts_remotes", "ROC_model_pure_R.R")
message(sprintf("Sourcing: %s", roc_script))
source(roc_script)

# ==============================================================================
# Configuration
# ==============================================================================

# Paths to data files (adjust if needed)
DATA_DIR <- file.path(dirname(script_dir), "data")
RESULTS_DIR <- file.path(dirname(script_dir), "results", "pure_R_ROC_test")

# Input files
FASTA_FILE <- file.path(DATA_DIR, "IM767_887_v2.1.cds_primaryTranscriptOnlyCleanFiltered.fa")
EXPRESSION_FILE <- file.path(DATA_DIR, "observed_expression_multitissue.csv")
DM_FILE <- file.path(DATA_DIR, "Mguttatus_intron_derived_dM.csv")

# Test parameters (small run for testing)
N_SAMPLES <- 100      # Increase for real analysis (e.g., 10000)
THIN <- 5             # Increase for real analysis (e.g., 10-100)
N_CORES <- 4          # Set to 0 for auto-detect

# ==============================================================================
# Check files exist
# ==============================================================================

message("============================================")
message("  Pure R ROC Model - Test Run")
message("============================================")
message("")

check_file <- function(path, desc) {
  if (file.exists(path)) {
    message(sprintf("✓ Found %s: %s", desc, basename(path)))
    return(TRUE)
  } else {
    message(sprintf("✗ Missing %s: %s", desc, path))
    return(FALSE)
  }
}

files_ok <- all(c(

check_file(FASTA_FILE, "FASTA file"),
check_file(EXPRESSION_FILE, "Expression file"),
check_file(DM_FILE, "dM file")
))

if (!files_ok) {
  stop("Some required files are missing. Please check paths.")
}

message("")

# ==============================================================================
# Load data
# ==============================================================================

message("Loading FASTA file...")
fasta_data <- read_fasta_codon_counts(FASTA_FILE)
message(sprintf("  Loaded %d genes", fasta_data$n_genes))
message(sprintf("  Example gene IDs: %s", paste(head(fasta_data$gene_ids, 3), collapse = ", ")))

message("")
message("Loading expression data...")
expr_df <- read.csv(EXPRESSION_FILE)
message(sprintf("  Columns: %s", paste(colnames(expr_df), collapse = ", ")))
message(sprintf("  Rows: %d", nrow(expr_df)))

# Use leaf expression (column 2) as phi
# Match gene IDs
gene_col <- colnames(expr_df)[1]
phi_col <- "Exp_leaf"  # or "Exp_bud"

matched_idx <- match(fasta_data$gene_ids, expr_df[[gene_col]])
n_matched <- sum(!is.na(matched_idx))
message(sprintf("  Matched %d / %d genes", n_matched, fasta_data$n_genes))

obs_phi <- expr_df[[phi_col]][matched_idx]
obs_phi[is.na(obs_phi)] <- 1  # Default for unmatched genes

# Normalize phi to have mean = 1 (standard for ROC model)
obs_phi <- obs_phi / mean(obs_phi, na.rm = TRUE)
message(sprintf("  Phi range: %.4f - %.4f (mean = %.4f)", 
                min(obs_phi), max(obs_phi), mean(obs_phi)))

message("")
message("Loading dM values...")
dM_df <- read.csv(DM_FILE)
message(sprintf("  Columns: %s", paste(colnames(dM_df), collapse = ", ")))
message(sprintf("  Rows: %d codons", nrow(dM_df)))

# Parse dM into list format
aa_to_codons <- get_genetic_code()
syn_aas <- get_synonymous_aas()

init_dM <- list()
for (aa in syn_aas) {
  codons <- aa_to_codons[[aa]]
  non_ref_codons <- codons[-length(codons)]  # Exclude reference (last alphabetically)
  
  values <- numeric(length(non_ref_codons))
  for (i in seq_along(non_ref_codons)) {
    row_idx <- which(dM_df$Codon == non_ref_codons[i])
    if (length(row_idx) > 0) {
      values[i] <- dM_df$dM[row_idx[1]]
    }
  }
  init_dM[[aa]] <- values
}

message(sprintf("  Parsed dM for %d amino acids", length(init_dM)))
message(sprintf("  Example (Phe): %s", paste(round(init_dM[["F"]], 3), collapse = ", ")))

# ==============================================================================
# Test 1: Quick run without phi (baseline)
# ==============================================================================

message("")
message("============================================")
message("TEST 1: Quick run without phi (baseline)")
message("============================================")

# Use only first 100 genes for quick test
test_genes <- min(100, fasta_data$n_genes)
test_counts <- fasta_data$codon_counts[1:test_genes, ]

results_baseline <- run_roc_mcmc(
  codon_counts = test_counts,
  n_samples = 50,
  thin = 2,
  obs_phi = NULL,
  fix_phi = FALSE,
  fix_dM = FALSE,
  fix_dEta = FALSE,
  n_cores = N_CORES,
  verbose = TRUE
)

message(sprintf("  Final sphi: %.3f", results_baseline$sphi_mean))
message("  Test 1 PASSED")

# ==============================================================================
# Test 2: Run with fixed dM (the scenario that crashes AnaCoDa)
# ==============================================================================

message("")
message("============================================")
message("TEST 2: Run with fixed dM and observed phi")
message("  (This is the scenario that crashes AnaCoDa)")
message("============================================")

results_fixed_dM <- run_roc_mcmc(
  codon_counts = test_counts,
  n_samples = 50,
  thin = 2,
  obs_phi = obs_phi[1:test_genes],
  fix_phi = FALSE,
  fix_dM = TRUE,
  fix_dEta = FALSE,
  init_dM = init_dM,
  n_cores = N_CORES,
  verbose = TRUE
)

message(sprintf("  Final sphi: %.3f", results_fixed_dM$sphi_mean))
message("  Test 2 PASSED - No crash!")

# ==============================================================================
# Test 3: Full run with fixed dM and fixed phi
# ==============================================================================

message("")
message("============================================")
message("TEST 3: Run with both dM and phi fixed")
message("============================================")

results_both_fixed <- run_roc_mcmc(
  codon_counts = test_counts,
  n_samples = 50,
  thin = 2,
  obs_phi = obs_phi[1:test_genes],
  fix_phi = TRUE,
  fix_dM = TRUE,
  fix_dEta = FALSE,
  init_dM = init_dM,
  n_cores = N_CORES,
  verbose = TRUE
)

message(sprintf("  Final sphi: %.3f", results_both_fixed$sphi_mean))
message("  Test 3 PASSED")

# ==============================================================================
# Save test results
# ==============================================================================

message("")
message("============================================")
message("Saving test results...")
message("============================================")

if (!dir.exists(RESULTS_DIR)) {
  dir.create(RESULTS_DIR, recursive = TRUE)
}

save_results(results_fixed_dM, RESULTS_DIR)

# Also save a summary
summary_df <- data.frame(
  Test = c("Baseline", "Fixed_dM", "Both_Fixed"),
  sphi_mean = c(results_baseline$sphi_mean, 
                results_fixed_dM$sphi_mean, 
                results_both_fixed$sphi_mean),
  final_log_posterior = c(tail(results_baseline$log_posterior, 1),
                          tail(results_fixed_dM$log_posterior, 1),
                          tail(results_both_fixed$log_posterior, 1))
)
write.csv(summary_df, file.path(RESULTS_DIR, "test_summary.csv"), row.names = FALSE)

message(sprintf("Results saved to: %s", RESULTS_DIR))

# ==============================================================================
# Plot trace (if in interactive mode)
# ==============================================================================

if (interactive()) {
  message("")
  message("Plotting trace...")
  
  par(mfrow = c(2, 2))
  
  # Log posterior trace
  plot(results_fixed_dM$log_posterior, type = "l", 
       main = "Log Posterior Trace", xlab = "Sample", ylab = "Log Posterior")
  
  # sphi trace
  plot(results_fixed_dM$sphi_samples, type = "l",
       main = "Sphi Trace", xlab = "Sample", ylab = "Sphi")
  
  # Example phi trace (first gene)
  plot(results_fixed_dM$phi_samples[, 1], type = "l",
       main = "Phi[1] Trace", xlab = "Sample", ylab = "Phi")
  
  # Phi: observed vs estimated
  plot(obs_phi[1:test_genes], results_fixed_dM$phi_mean,
       main = "Observed vs Estimated Phi", 
       xlab = "Observed Phi", ylab = "Estimated Phi")
  abline(0, 1, col = "red", lty = 2)
  
  par(mfrow = c(1, 1))
}

message("")
message("============================================")
message("ALL TESTS PASSED!")
message("============================================")
message("")
message("The pure R implementation successfully handles the scenario")
message("that causes AnaCoDa to crash (fix_dM + with.phi).")
message("")
message("For production runs, increase n_samples and thin parameters.")
