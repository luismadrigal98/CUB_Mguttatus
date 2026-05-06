#!/usr/bin/env Rscript
#' ============================================================================
#' ROC-SEMPPR vs Polymorphism Consistency Tests
#' ============================================================================
#' 
#' This script tests whether the ROC-SEMPPR model results are consistent with
#' Site Frequency Spectrum (SFS) based selection estimates (γ = 4Nes).
#' 
#' KEY QUESTION: Why does ROC show expression-dependent selection intensity
#' (ROC_eff varies 100-1000× with expression) while SFS shows uniform γ (~1.4)?
#' 
#' CRITICAL FIX: The original γ estimation used global α,β from introns.
#' For proper per-amino-acid analysis, we need AA-specific parameters.
#' 
#' @author Luis Javier Madrigal-Roca
#' @date 2025-02-02
#' ============================================================================

# ==============================================================================
# SETUP
# ==============================================================================

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════════════╗\n")
cat("║     ROC-SEMPPR vs Polymorphism Consistency Analysis                   ║\n")
cat("║     Testing whether fixed and segregating site data are consistent    ║\n")
cat("╚═══════════════════════════════════════════════════════════════════════╝\n\n")

# Set working directory - use HPC path
work_dir <- "/home/l338m483/scratch/CUB/CUB_Mguttatus"
setwd(work_dir)
cat(sprintf("Working directory: %s\n\n", work_dir))

# Load required libraries
required_packages <- c(
  "data.table", "dplyr", "ggplot2", "cowplot",
  "mgcv", "gsl", "Biostrings"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("WARNING: Package '%s' not available, some analyses may fail\n", pkg))
  }
}

# Source required functions
source("./src/derivation_gamma_from_polymorphism.R")
source("./src/integrate_intronic_polymorphism.R")
source("./src/roc_polymorphism_consistency.R")
source("./src/roc_model_validation.R")
if (file.exists("./src/theme_custom.R")) {
  source("./src/theme_custom.R")
}

# Create output directories
output_dir <- "./results/consistency"
figures_dir <- file.path(output_dir, "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Output directory: %s\n", output_dir))
cat(sprintf("Figures directory: %s\n\n", figures_dir))

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading data...\n\n")

# 1. Load existing gamma vs ROC_eff data (already computed!)
gamma_sroc_file <- "./results/gamma_vs_Sroc_scaling_data.csv"

if (file.exists(gamma_sroc_file)) {
  gamma_sroc_data <- fread(gamma_sroc_file)
  cat(sprintf("✓ Loaded gamma vs ROC_eff data: %d expression bins\n", nrow(gamma_sroc_data)))
  has_gamma_data <- TRUE
} else {
  cat("⚠ gamma_vs_Sroc_scaling_data.csv not found\n")
  has_gamma_data <- FALSE
}

# 2. Load neutral parameters (already computed!)
neutral_params_file <- "./results/neutral_mutation_parameters.csv"

if (file.exists(neutral_params_file)) {
  np_df <- fread(neutral_params_file)
  neutral_params <- list(
    alpha_G = np_df[Parameter == "alpha_G", Value],
    beta_G = np_df[Parameter == "beta_G", Value],
    alpha_C = np_df[Parameter == "alpha_C", Value],
    beta_C = np_df[Parameter == "beta_C", Value],
    pi_G_expected = np_df[Parameter == "pi_G_expected", Value],
    pi_C_expected = np_df[Parameter == "pi_C_expected", Value]
  )
  cat(sprintf("✓ Loaded neutral parameters from CSV\n"))
  cat(sprintf("  α_G = %.6f, β_G = %.6f\n", neutral_params$alpha_G, neutral_params$beta_G))
  cat(sprintf("  α_C = %.6f, β_C = %.6f\n", neutral_params$alpha_C, neutral_params$beta_C))
} else {
  # Calculate from intronic SFS
  cat("\n--- Calculating neutral parameters from intronic SFS ---\n")
  sfs_G_file <- "./data/sfs_introns_G.csv"
  sfs_C_file <- "./data/sfs_introns_C.csv"
  
  if (!file.exists(sfs_G_file) || !file.exists(sfs_C_file)) {
    stop("Intronic SFS files not found.")
  }
  neutral_params <- load_and_estimate_neutral_params(sfs_G_file, sfs_C_file)
}

# 3. Load expression data
cat("\n--- Loading expression data ---\n")
exp_file <- "./data/compiled_expression_IM767.txt"

if (file.exists(exp_file)) {
  exp_complete <- fread(exp_file)
  
  # Get gene name column
  gene_col <- names(exp_complete)[1]
  setnames(exp_complete, gene_col, "Gene_name")
  
  # Calculate expression metrics
  numeric_cols <- setdiff(names(exp_complete), "Gene_name")
  exp_complete[, Mean_Log10_Exp := rowMeans(log10(.SD + 1)), .SDcols = numeric_cols]
  exp_complete[, Max_Log10_Exp := apply(log10(.SD + 1), 1, max), .SDcols = numeric_cols]
  
  # Define expression groups
  top_5_cutoff <- quantile(exp_complete$Max_Log10_Exp, probs = 0.95)
  bottom_5_cutoff <- quantile(exp_complete$Max_Log10_Exp, probs = 0.05)
  
  exp_complete[, Expression_Group := fifelse(
    Max_Log10_Exp >= top_5_cutoff, "Top 5%",
    fifelse(Max_Log10_Exp <= bottom_5_cutoff, "Bottom 5%", "Middle 90%")
  )]
  
  integrated_data <- exp_complete
  cat(sprintf("✓ Loaded expression data: %d genes\n", nrow(integrated_data)))
} else {
  stop("Expression file not found: ", exp_file)
}

# 4. Load preferred codons from ROC model
cat("\n--- Loading preferred codons ---\n")
pref_codons_file <- "./results/preferred_codons.txt"

if (file.exists(pref_codons_file)) {
  preferred_codons_vec <- scan(pref_codons_file, what = "character", quiet = TRUE)
  
  # Create data frame with AA mapping
  genetic_code <- Biostrings::GENETIC_CODE
  preferred_codons_df <- data.table(
    Codon = preferred_codons_vec,
    AA = as.character(genetic_code[preferred_codons_vec])
  )
  preferred_codons_df <- preferred_codons_df[AA != "*"]
  
  cat(sprintf("✓ Loaded %d preferred codons\n", nrow(preferred_codons_df)))
} else {
  stop("Preferred codons file not found.")
}

# 5. Get per-AA neutral parameters using Q-matrix approach
cat("\n--- Computing per-amino-acid neutral parameters ---\n")

# Calculate GC content from introns for Q-matrix
# Using alpha/(alpha+beta) as proxy for equilibrium frequency
pi_G <- neutral_params$alpha_G / (neutral_params$alpha_G + neutral_params$beta_G)
pi_C <- neutral_params$alpha_C / (neutral_params$alpha_C + neutral_params$beta_C)

# Assume symmetric for A/T
pi_A <- (1 - pi_G - pi_C) / 2
pi_T <- (1 - pi_G - pi_C) / 2

cat(sprintf("  Equilibrium frequencies: π_A=%.3f, π_C=%.3f, π_G=%.3f, π_T=%.3f\n",
            pi_A, pi_C, pi_G, pi_T))

# Build Q-matrix
Q <- solve_Q_matrix_for_consistency(pi_A, pi_C, pi_G, pi_T, kappa = 2)

# Calculate theta (4N*mu) from intronic data
theta_intron <- (neutral_params$alpha_G + neutral_params$beta_G + 
                 neutral_params$alpha_C + neutral_params$beta_C) / 4
cat(sprintf("  Estimated θ (from introns): %.6f\n", theta_intron))

# Get per-AA mutation rates
aa_params_qmatrix <- calculate_aa_specific_alpha_beta(
  Q = Q, 
  preferred_codons = preferred_codons_vec, 
  genetic_code = genetic_code,
  theta_intron = theta_intron
)

# Also get simple terminal-nucleotide based params
aa_params_simple <- get_aa_specific_neutral_params(neutral_params, preferred_codons_df)

# 6. Load CSP parameters from ROC model
cat("\n--- Loading CSP parameters ---\n")
csp_file <- "./results/MCMC_results/results_dM_fixed/run_1/Parameter_est/Cluster_1_Selection.csv"

if (file.exists(csp_file)) {
  mutation_file <- gsub("Selection", "Mutation", csp_file)
  csp_df <- load_csp_parameters(mutation_file, csp_file)
  cat(sprintf("✓ Loaded CSP parameters: %d codons\n", nrow(csp_df)))
  has_csp_data <- TRUE
} else {
  cat("⚠ CSP parameters not found\n")
  has_csp_data <- FALSE
}

# ==============================================================================
# INITIALIZE RESULTS LIST
# ==============================================================================

results_all <- list()

# ==============================================================================
# ANALYSIS 1: Re-estimate Gamma per Amino Acid with CORRECT α,β
# ==============================================================================

cat("\n")
cat(rep("█", 70), "\n", sep = "")
cat("ANALYSIS 1: Re-estimate Gamma per Amino Acid with CORRECT Parameters\n")
cat(rep("█", 70), "\n\n")

cat("THE PROBLEM:\n")
cat("  Previous gamma estimates used GLOBAL α,β from introns for ALL amino acids.\n")
cat("  But different AAs have G-ending vs C-ending preferred codons,\n")
cat("  which have DIFFERENT mutation parameters!\n\n")

cat("THE FIX:\n")
cat(sprintf("  G-ending preferred codons → α_G = %.6f, β_G = %.6f\n",
            neutral_params$alpha_G, neutral_params$beta_G))
cat(sprintf("  C-ending preferred codons → α_C = %.6f, β_C = %.6f\n\n",
            neutral_params$alpha_C, neutral_params$beta_C))

# Load raw frequency data per amino acid
raw_freq_file <- "./data/all_chromosomes.raw_frequencies.txt"

if (file.exists(raw_freq_file)) {
  
  cat("Loading raw frequencies per amino acid...\n")
  raw_freq <- fread(raw_freq_file, skip = 3)  # Skip comment lines
  setnames(raw_freq, c("AA", "Family", "Preferred_Freq"))
  
  # Filter to 4-fold and 2-fold sites (synonymous)
  raw_freq <- raw_freq[Family %in% c("2-fold", "4-fold")]
  
  cat(sprintf("  Loaded %d polymorphic sites\n", nrow(raw_freq)))
  cat(sprintf("  Amino acids: %s\n", paste(sort(unique(raw_freq$AA)), collapse = ", ")))
  
  # Get sample size from intronic SFS (most common n)
  sfs_G <- fread("./data/sfs_introns_G.csv")
  n_sample <- sfs_G[which.max(count), n]
  cat(sprintf("  Sample size (n) from introns: %d chromosomes\n\n", n_sample))
  
  # Convert frequencies to (k, n) format
  # k = round(Preferred_Freq * n)
  raw_freq[, k := round(Preferred_Freq * n_sample)]
  raw_freq[, n := n_sample]
  
  # Add terminal nucleotide from preferred codons
  # Map AA single letter to terminal nucleotide of preferred codon
  aa_terminal <- preferred_codons_df[, .(AA, Terminal_Nuc = substr(Codon, 3, 3))]
  raw_freq <- merge(raw_freq, aa_terminal, by = "AA", all.x = TRUE)
  
  # Check which AAs have which terminal nucleotide
  cat("=== Amino Acid Terminal Nucleotide Mapping ===\n")
  aa_summary <- raw_freq[, .(N_Sites = .N), by = .(AA, Terminal_Nuc)]
  print(aa_summary[order(Terminal_Nuc, AA)])
  
  # Estimate gamma for each amino acid using CORRECT α,β
  cat("\n=== Re-estimating Gamma per Amino Acid ===\n\n")
  
  gamma_per_aa <- raw_freq[, {
    
    aa_name <- AA[1]
    term_nuc <- Terminal_Nuc[1]
    
    # Get correct alpha/beta based on terminal nucleotide
    if (is.na(term_nuc)) {
      list(
        Gamma = NA_real_,
        SE = NA_real_,
        N_Sites = .N,
        Alpha_Used = NA_real_,
        Beta_Used = NA_real_
      )
    } else if (term_nuc == "G") {
      alpha <- neutral_params$alpha_G
      beta <- neutral_params$beta_G
    } else if (term_nuc == "C") {
      alpha <- neutral_params$alpha_C
      beta <- neutral_params$beta_C
    } else {
      alpha <- (neutral_params$alpha_G + neutral_params$alpha_C) / 2
      beta <- (neutral_params$beta_G + neutral_params$beta_C) / 2
    }
    
    # Estimate gamma using MLE
    gamma_est <- tryCatch({
      estimate_gamma_for_AA(
        counts = k,
        sample_sizes = n,
        alpha = alpha,
        beta = beta,
        S_interval = c(-10, 50)
      )
    }, error = function(e) {
      cat(sprintf("  Error for %s: %s\n", aa_name, e$message))
      NA_real_
    })
    
    list(
      Gamma = gamma_est,
      SE = NA_real_,  # Would need bootstrap
      N_Sites = .N,
      Alpha_Used = alpha,
      Beta_Used = beta
    )
  }, by = .(AA, Terminal_Nuc)]
  
  # Print results
  cat("\n=== Gamma Estimates per Amino Acid (CORRECTED) ===\n")
  print(gamma_per_aa[order(-Gamma)])
  
  # Summary by terminal nucleotide
  summary_by_nuc <- gamma_per_aa[!is.na(Gamma), .(
    N_AA = .N,
    Mean_Gamma = mean(Gamma),
    Median_Gamma = median(Gamma),
    SD_Gamma = sd(Gamma),
    Min_Gamma = min(Gamma),
    Max_Gamma = max(Gamma)
  ), by = Terminal_Nuc]
  
  cat("\n=== Summary by Terminal Nucleotide ===\n")
  print(summary_by_nuc)
  
  # Statistical test: G vs C ending
  g_gammas <- gamma_per_aa[Terminal_Nuc == "G" & !is.na(Gamma), Gamma]
  c_gammas <- gamma_per_aa[Terminal_Nuc == "C" & !is.na(Gamma), Gamma]
  
  if (length(g_gammas) >= 3 && length(c_gammas) >= 3) {
    wilcox_test <- wilcox.test(g_gammas, c_gammas)
    t_test <- t.test(g_gammas, c_gammas)
    
    cat(sprintf("\nG vs C ending codons:\n"))
    cat(sprintf("  Wilcoxon p-value: %.4f\n", wilcox_test$p.value))
    cat(sprintf("  t-test p-value: %.4f\n", t_test$p.value))
  }
  
  # Visualization
  p_aa <- ggplot(gamma_per_aa[!is.na(Gamma)], 
                 aes(x = reorder(AA, Gamma), y = Gamma, fill = Terminal_Nuc)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 1, linetype = "dotted", color = "red") +
    scale_fill_manual(values = c("G" = "#E41A1C", "C" = "#377EB8"),
                      name = "Terminal\nNucleotide",
                      na.value = "gray50") +
    coord_flip() +
    labs(
      title = "Selection Coefficient (γ) by Amino Acid - CORRECTED",
      subtitle = "Using AA-specific α,β based on terminal nucleotide",
      x = "Amino Acid",
      y = expression(gamma ~ "(4N"[e]*"s)")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  ggsave(file.path(figures_dir, "Analysis1_Gamma_per_AA_CORRECTED.pdf"),
         p_aa, width = 8, height = 10)
  cat(sprintf("\n✓ Plot saved: %s\n", file.path(figures_dir, "Analysis1_Gamma_per_AA_CORRECTED.pdf")))
  
  results_all$analysis1 <- list(
    gamma_per_aa = gamma_per_aa,
    summary_by_nuc = summary_by_nuc,
    raw_freq = raw_freq,
    plot = p_aa
  )
  
} else {
  cat("⚠ Raw frequency file not found: ", raw_freq_file, "\n")
}

# ==============================================================================
# INTERPRETATION AND RECOMMENDATIONS
# ==============================================================================

cat("\n")
cat(rep("═", 70), "\n", sep = "")
cat("INTERPRETATION\n")
cat(rep("═", 70), "\n\n")

cat("KEY METHODOLOGICAL ISSUE IDENTIFIED:\n")
cat("────────────────────────────────────\n")
cat("The original γ estimation used GLOBAL α,β from introns for ALL amino acids.\n")
cat("However, different amino acids have G-ending vs C-ending preferred codons,\n")
cat("and G and C have DIFFERENT mutation parameters!\n\n")

cat("THE FIX:\n")
cat("────────\n")
cat(sprintf("- G-ending preferred codons: Use α_G = %.4f, β_G = %.4f\n",
    neutral_params$alpha_G, neutral_params$beta_G))
cat(sprintf("- C-ending preferred codons: Use α_C = %.4f, β_C = %.4f\n\n",
    neutral_params$alpha_C, neutral_params$beta_C))

cat("POSSIBLE EXPLANATIONS FOR APPARENT INCONSISTENCY:\n")
cat("──────────────────────────────────────────────────\n")
cat("1. METHODOLOGICAL: Using global α,β instead of AA-specific values\n")
cat("   → Fixed by per-AA γ estimation in Analysis 1\n\n")

cat("2. BIOLOGICAL: Different timescales\n")
cat("   → ROC measures accumulated fixed differences (long timescale)\n")
cat("   → SFS measures current segregating variation (short timescale)\n\n")

cat("3. gBGC CONFOUNDING: GC-biased gene conversion\n")
cat("   → Creates apparent selection signal independent of expression\n")
cat("   → Would explain uniform γ across expression levels\n\n")

cat("4. SAMPLE SIZE: Low polymorphism may limit detection power\n")
cat("   → Especially for genes with few segregating sites per AA\n\n")

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

cat("\nSaving results...\n")

# Save all results
saveRDS(results_all, file.path(output_dir, "consistency_analysis_results.rds"))
cat(sprintf("✓ Results saved: %s\n", 
            file.path(output_dir, "consistency_analysis_results.rds")))

# Save per-AA gamma results if available
if (!is.null(results_all$analysis1)) {
  write.csv(
    results_all$analysis1$gamma_per_aa,
    file.path(output_dir, "gamma_per_aa_CORRECTED.csv"),
    row.names = FALSE
  )
  cat(sprintf("✓ Per-AA gamma saved: %s\n",
              file.path(output_dir, "gamma_per_aa_CORRECTED.csv")))
}

cat("\n")
cat(rep("═", 70), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(sprintf("Output directory: %s\n", output_dir))
cat(rep("═", 70), "\n")
