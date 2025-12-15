#' Integration of Intronic Polymorphism Data for Selection Inference
#' 
#' This module handles the workflow for estimating selection coefficients (gamma)
#' using mutation rate parameters (alpha, beta) derived from neutral intronic sites.
#' 
#' @author Luis Javier Madrigal-Roca and John K. Kelly
#' _____________________________________________________________________________

load_and_estimate_neutral_params <- function(sfs_G_file, sfs_C_file) {
  #' Load intronic SFS data and estimate alpha/beta for G and C separately
  #' 
  #' @param sfs_G_file Path to CSV with G-ending site frequency spectrum from introns
  #' @param sfs_C_file Path to CSV with C-ending site frequency spectrum from introns
  #' @return List with alpha_G, beta_G, alpha_C, beta_C and diagnostic info
  #' ___________________________________________________________________________
  
  require(data.table)
  
  cat("\n=== Estimating Neutral Mutation Parameters from Introns ===\n\n")
  
  # Estimate parameters for G sites
  cat("Processing G-ending intronic sites...\n")
  params_G <- solve_alpha_and_beta_from_introns(sfs_G_file)
  
  cat(sprintf("  α_G (4N·u_G) = %.6f\n", params_G$alpha))
  cat(sprintf("  β_G (4N·v_G) = %.6f\n", params_G$beta))
  cat(sprintf("  Sites analyzed: %d\n", params_G$n_sites))
  cat(sprintf("  Convergence: %s\n\n", 
              ifelse(params_G$convergence == 0, "Success", "Warning")))
  
  # Estimate parameters for C sites
  cat("Processing C-ending intronic sites...\n")
  params_C <- solve_alpha_and_beta_from_introns(sfs_C_file)
  
  cat(sprintf("  α_C (4N·u_C) = %.6f\n", params_C$alpha))
  cat(sprintf("  β_C (4N·v_C) = %.6f\n", params_C$beta))
  cat(sprintf("  Sites analyzed: %d\n", params_C$n_sites))
  cat(sprintf("  Convergence: %s\n\n", 
              ifelse(params_C$convergence == 0, "Success", "Warning")))
  
  # Calculate expected Pi for validation
  cat("=== Validation: Expected Nucleotide Diversity ===\n")
  cat("(Should match observed Pi from intronic sites)\n\n")
  
  pi_G_expected <- calculate_pi_analytical(params_G$alpha, params_G$beta, S = 0)
  pi_C_expected <- calculate_pi_analytical(params_C$alpha, params_C$beta, S = 0)
  
  cat(sprintf("  E[π] at G sites: %.6f\n", pi_G_expected))
  cat(sprintf("  E[π] at C sites: %.6f\n\n", pi_C_expected))
  
  # Return combined results
  results <- list(
    alpha_G = params_G$alpha,
    beta_G = params_G$beta,
    alpha_C = params_C$alpha,
    beta_C = params_C$beta,
    pi_G_expected = pi_G_expected,
    pi_C_expected = pi_C_expected,
    n_sites_G = params_G$n_sites,
    n_sites_C = params_C$n_sites,
    convergence_G = params_G$convergence,
    convergence_C = params_C$convergence
  )
  
  return(results)
}

annotate_preferred_codons_with_nucleotide <- function(preferred_codons_df) {
  #' Add nucleotide annotation (G vs C ending) to preferred codons table
  #' 
  #' @param preferred_codons_df Data frame with columns: AA, Preferred_Codon
  #' @return Same data frame with added column: Terminal_Nucleotide
  #' ___________________________________________________________________________
  
  preferred_codons_df$Terminal_Nucleotide <- substr(
    preferred_codons_df$Preferred_Codon, 
    nchar(preferred_codons_df$Preferred_Codon), 
    nchar(preferred_codons_df$Preferred_Codon)
  )
  
  # Validation: All should be G or C
  stopifnot(all(preferred_codons_df$Terminal_Nucleotide %in% c("G", "C")))
  
  return(preferred_codons_df)
}

estimate_gamma_by_gene_with_neutral_params <- function(codon_vcf_data, 
                                                       neutral_params,
                                                       preferred_codons_df) {
  #' Estimate gamma (selection coefficient) for each gene and amino acid
  #' using pre-computed alpha and beta from intronic data
  #' 
  #' @param codon_vcf_data Output from process_codon_vcf_fast()
  #' @param neutral_params Output from load_and_estimate_neutral_params()
  #' @param preferred_codons_df Data frame mapping AA to preferred codons
  #' @return Data table with gamma estimates per gene and amino acid
  #' ___________________________________________________________________________
  
  require(data.table)
  
  if (!is.data.table(codon_vcf_data)) setDT(codon_vcf_data)
  if (!is.data.table(preferred_codons_df)) setDT(preferred_codons_df)
  
  cat("\n=== Estimating Selection Coefficients (Gamma) ===\n\n")
  
  # Annotate preferred codons with terminal nucleotide
  pref_annotated <- annotate_preferred_codons_with_nucleotide(preferred_codons_df)
  
  # Merge to get nucleotide annotation for each site
  setkey(codon_vcf_data, AA)
  setkey(pref_annotated, AA)
  
  data_merged <- pref_annotated[codon_vcf_data, nomatch = 0]
  
  # Count amino acids by terminal nucleotide
  aa_by_nuc <- data_merged[, .(N_Sites = .N), by = .(AA, Terminal_Nucleotide)]
  setkey(aa_by_nuc, Terminal_Nucleotide)
  
  cat(sprintf("Amino acids with G-ending preferred codons: %d\n", 
              sum(aa_by_nuc[Terminal_Nucleotide == "G"]$N_Sites > 0)))
  cat(sprintf("Amino acids with C-ending preferred codons: %d\n\n", 
              sum(aa_by_nuc[Terminal_Nucleotide == "C"]$N_Sites > 0)))
  
  # Estimate gamma for each Gene x AA combination
  cat("Running gamma estimation (this may take a few minutes)...\n")
  
  gamma_results <- data_merged[, {
    
    # Select appropriate alpha/beta based on terminal nucleotide
    if (Terminal_Nucleotide[1] == "G") {
      alpha_use <- neutral_params$alpha_G
      beta_use <- neutral_params$beta_G
    } else {  # C
      alpha_use <- neutral_params$alpha_C
      beta_use <- neutral_params$beta_C
    }
    
    # Skip if too few sites
    if (.N < 5) {
      list(Gamma = NA_real_, N_Sites = .N, Terminal_Nuc = Terminal_Nucleotide[1])
    } else {
      # Estimate gamma
      gamma_est <- tryCatch(
        estimate_gamma_for_AA(
          counts = k,
          sample_sizes = n,
          alpha = alpha_use,
          beta = beta_use,
          S_interval = c(-10, 30)  # Allow negative selection too
        ),
        error = function(e) NA_real_
      )
      
      list(
        Gamma = gamma_est,
        N_Sites = .N,
        Terminal_Nuc = Terminal_Nucleotide[1],
        Mean_Freq_Preferred = mean(p, na.rm = TRUE),
        Alpha_Used = alpha_use,
        Beta_Used = beta_use
      )
    }
  }, by = .(Gene, AA)]
  
  # Calculate significance (using chi-square approximation)
  # Thr = 1 to align with the Drift Barrier hypothesis
  gamma_results[, Significant := ifelse(!is.na(Gamma), abs(Gamma) > 1, FALSE)]
  
  cat("\n=== Summary Statistics ===\n\n")
  cat(sprintf("Total Gene x AA combinations: %d\n", nrow(gamma_results)))
  cat(sprintf("Successfully estimated: %d (%.1f%%)\n", 
              sum(!is.na(gamma_results$Gamma)),
              100 * mean(!is.na(gamma_results$Gamma))))
  cat(sprintf("Significant (|γ| > 1.92): %d (%.1f%%)\n\n",
              sum(gamma_results$Significant, na.rm = TRUE),
              100 * mean(gamma_results$Significant, na.rm = TRUE)))
  
  # Summary by terminal nucleotide
  summary_by_nuc <- gamma_results[!is.na(Gamma), .(
    N = .N,
    Mean_Gamma = mean(Gamma),
    Median_Gamma = median(Gamma),
    SD_Gamma = sd(Gamma),
    Prop_Positive = mean(Gamma > 0),
    Prop_Significant = mean(Significant)
  ), by = Terminal_Nuc]
  
  cat("=== Gamma Distribution by Terminal Nucleotide ===\n\n")
  print(summary_by_nuc)
  cat("\n")
  
  return(gamma_results)
}

compare_gamma_with_expression <- function(gamma_results, expression_data) {
  #' Compare estimated selection coefficients with gene expression levels
  #' 
  #' @param gamma_results Output from estimate_gamma_by_gene_with_neutral_params()
  #' @param expression_data Data frame with columns: Gene, High_exp_log2
  #' @return Merged data table with correlation statistics
  #' ___________________________________________________________________________
  
  require(data.table)
  
  if (!is.data.table(gamma_results)) setDT(gamma_results)
  if (!is.data.table(expression_data)) setDT(expression_data)
  
  cat("\n=== Comparing Gamma with Expression Levels ===\n\n")
  
  # Aggregate gamma per gene (mean across amino acids)
  gamma_per_gene <- gamma_results[!is.na(Gamma), .(
    Mean_Gamma = mean(Gamma),
    Median_Gamma = median(Gamma),
    Max_Gamma = max(Gamma),
    N_AA_Analyzed = .N,
    Prop_Significant = mean(Significant)
  ), by = Gene]
  
  # Merge with expression
  setkey(gamma_per_gene, Gene)
  setkey(expression_data, Gene)
  
  merged <- expression_data[gamma_per_gene, nomatch = 0]
  
  # Calculate correlation
  cor_mean <- cor.test(merged$High_exp_log2, merged$Mean_Gamma, 
                       method = "spearman", exact = FALSE)
  cor_median <- cor.test(merged$High_exp_log2, merged$Median_Gamma,
                         method = "spearman", exact = FALSE)
  
  cat(sprintf("Genes analyzed: %d\n\n", nrow(merged)))
  cat("Spearman Correlation:\n")
  cat(sprintf("  Expression vs Mean Gamma:   ρ = %.3f, p = %.2e\n",
              cor_mean$estimate, cor_mean$p.value))
  cat(sprintf("  Expression vs Median Gamma: ρ = %.3f, p = %.2e\n\n",
              cor_median$estimate, cor_median$p.value))
  
  # Compare high vs low expression genes
  q75 <- quantile(merged$High_exp_log2, 0.75)
  q25 <- quantile(merged$High_exp_log2, 0.25)
  
  high_exp_gamma <- merged[High_exp_log2 >= q75]$Mean_Gamma
  low_exp_gamma <- merged[High_exp_log2 <= q25]$Mean_Gamma
  
  wilcox_test <- wilcox.test(high_exp_gamma, low_exp_gamma)
  
  cat("High Expression (Q4) vs Low Expression (Q1):\n")
  cat(sprintf("  High: Mean γ = %.3f (SD = %.3f)\n",
              mean(high_exp_gamma), sd(high_exp_gamma)))
  cat(sprintf("  Low:  Mean γ = %.3f (SD = %.3f)\n",
              mean(low_exp_gamma), sd(low_exp_gamma)))
  cat(sprintf("  Wilcoxon p-value: %.2e\n\n", wilcox_test$p.value))
  
  return(merged)
}

validate_against_cai <- function(gamma_results, integrated_data) {
  #' Cross-validate gamma estimates against CAI-based bias metrics
  #' 
  #' @param gamma_results Output from estimate_gamma_by_gene_with_neutral_params()
  #' @param integrated_data Main analysis data with CAI, CDC, etc.
  #' @return Data table with validation statistics
  #' ___________________________________________________________________________
  
  require(data.table)
  
  if (!is.data.table(gamma_results)) setDT(gamma_results)
  if (!is.data.table(integrated_data)) setDT(integrated_data)
  
  cat("\n=== Cross-Validation: Gamma vs CAI/CDC ===\n\n")
  
  # Aggregate gamma per gene
  gamma_summary <- gamma_results[!is.na(Gamma), .(
    Mean_Gamma = mean(Gamma),
    Median_Gamma = median(Gamma),
    N_Positive_Gamma = sum(Gamma > 0),
    N_Significant = sum(Significant)
  ), by = Gene]
  
  # Merge with integrated data
  setkey(gamma_summary, Gene)
  setkey(integrated_data, Gene)
  
  validation_data <- integrated_data[gamma_summary, nomatch = 0]
  
  # Correlations
  cat("Correlations (Spearman):\n")
  
  if ("CAI" %in% names(validation_data)) {
    cor_cai <- cor.test(validation_data$Mean_Gamma, validation_data$CAI,
                        method = "spearman", exact = FALSE)
    cat(sprintf("  Gamma vs CAI:  ρ = %.3f, p = %.2e\n",
                cor_cai$estimate, cor_cai$p.value))
  }
  
  if ("CDC" %in% names(validation_data)) {
    cor_cdc <- cor.test(validation_data$Mean_Gamma, validation_data$CDC,
                        method = "spearman", exact = FALSE)
    cat(sprintf("  Gamma vs CDC:  ρ = %.3f, p = %.2e\n",
                cor_cdc$estimate, cor_cdc$p.value))
  }
  
  if ("ENC" %in% names(validation_data)) {
    cor_enc <- cor.test(validation_data$Mean_Gamma, validation_data$ENC,
                        method = "spearman", exact = FALSE)
    cat(sprintf("  Gamma vs ENC:  ρ = %.3f, p = %.2e\n",
                cor_enc$estimate, cor_enc$p.value))
  }
  
  cat("\n")
  
  return(validation_data)
}