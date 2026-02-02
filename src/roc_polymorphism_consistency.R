#' ROC-SEMPPR vs Polymorphism Consistency Analysis Functions
#' 
#' Functions to test whether ROC-SEMPPR model results are consistent with
#' Site Frequency Spectrum (SFS) based selection estimates.
#' 
#' KEY INSIGHT: The original gamma estimation pooled all amino acids together,
#' but different amino acids have different mutation rates (alpha, beta).
#' For proper comparison with ROC model, gamma must be estimated PER amino acid
#' using amino-acid-specific alpha/beta values.
#' 
#' @author Luis Javier Madrigal-Roca
#' _____________________________________________________________________________

#' Calculate Q-matrix (mutation rate matrix) from nucleotide composition
#' 
#' Uses HKY85 model where transitions are weighted by kappa.
#' 
#' @param pi_A Frequency of adenine
#' @param pi_C Frequency of cytosine
#' @param pi_G Frequency of guanine  
#' @param pi_T Frequency of thymine
#' @param kappa Transition/transversion ratio (default: 2)
#' @return 4x4 normalized Q-matrix
solve_Q_matrix_for_consistency <- function(pi_A, pi_C, pi_G, pi_T, kappa = 2) {
  
  Q <- matrix(0, nrow = 4, ncol = 4, 
              dimnames = list(c("A","C","G","T"), c("A","C","G","T")))
  
  freqs <- c(A=pi_A, C=pi_C, G=pi_G, T=pi_T)
  bases <- c("A", "C", "G", "T")
  
  for (i in bases) {
    for (j in bases) {
      if (i == j) next
      
      is_transition <- (i=="A" & j=="G") | (i=="G" & j=="A") | 
        (i=="C" & j=="T") | (i=="T" & j=="C")
      
      rate <- freqs[j]
      if (is_transition) rate <- rate * kappa
      Q[i, j] <- rate
    }
  }
  
  diag(Q) <- -rowSums(Q)
  scaling_factor <- -sum(freqs * diag(Q))
  Q_normalized <- Q / scaling_factor
  
  return(Q_normalized)
}

#' Calculate per-amino-acid mutation rates (u, v) from Q-matrix
#' 
#' For each amino acid family, calculates:
#' - u: mutation rate from unpreferred to preferred codons
#' - v: mutation rate from preferred to unpreferred codons
#' 
#' @param Q 4x4 mutation rate matrix
#' @param preferred_codons Vector of preferred codons
#' @param genetic_code Named vector mapping codons to amino acids
#' @return Data.table with AA, u, v columns
get_aa_mutation_rates_from_Q <- function(Q, preferred_codons, genetic_code) {
  
  require(data.table)
  
  # Group codons by amino acid
  aa_to_codons <- split(names(genetic_code), genetic_code)
  aa_to_codons <- aa_to_codons[!names(aa_to_codons) %in% c("*", "STOP")]
  
  results <- list()
  
  for (aa in names(aa_to_codons)) {
    codons <- aa_to_codons[[aa]]
    if (length(codons) < 2) next  # Skip Met, Trp
    
    pref <- intersect(codons, preferred_codons)
    unpref <- setdiff(codons, pref)
    
    if (length(pref) == 0 || length(unpref) == 0) next
    
    # Calculate u (Unpreferred -> Preferred)
    u_rates <- c()
    for (c_un in unpref) {
      rate_out <- 0
      for (c_p in pref) {
        diffs <- 0
        nuc_from <- ""
        nuc_to <- ""
        for (i in 1:3) {
          if (substr(c_un, i, i) != substr(c_p, i, i)) {
            diffs <- diffs + 1
            nuc_from <- substr(c_un, i, i)
            nuc_to <- substr(c_p, i, i)
          }
        }
        if (diffs == 1) {
          rate_out <- rate_out + Q[nuc_from, nuc_to]
        }
      }
      u_rates <- c(u_rates, rate_out)
    }
    u <- mean(u_rates)
    
    # Calculate v (Preferred -> Unpreferred)
    v_rates <- c()
    for (c_p in pref) {
      rate_out <- 0
      for (c_un in unpref) {
        diffs <- 0
        nuc_from <- ""
        nuc_to <- ""
        for (i in 1:3) {
          if (substr(c_p, i, i) != substr(c_un, i, i)) {
            diffs <- diffs + 1
            nuc_from <- substr(c_p, i, i)
            nuc_to <- substr(c_un, i, i)
          }
        }
        if (diffs == 1) {
          rate_out <- rate_out + Q[nuc_from, nuc_to]
        }
      }
      v_rates <- c(v_rates, rate_out)
    }
    v <- mean(v_rates)
    
    results[[aa]] <- data.table(
      AA = aa,
      Preferred_Codon = paste(pref, collapse = ","),
      u = u,
      v = v
    )
  }
  
  rbindlist(results)
}

#' Calculate per-amino-acid alpha and beta from Q-matrix and theta
#' 
#' This is the PROPER way to get per-AA neutral parameters:
#' 1. Use Q-matrix to get relative mutation rates (u, v) per AA
#' 2. Scale by global theta (from introns) to get alpha = 4N*u, beta = 4N*v
#' 
#' @param Q 4x4 mutation rate matrix
#' @param preferred_codons Vector of preferred codons
#' @param genetic_code Named vector mapping codons to amino acids
#' @param theta_intron Global theta from intronic data (average of G and C)
#' @return Data.table with per-AA alpha and beta values
calculate_aa_specific_alpha_beta <- function(Q, preferred_codons, genetic_code, 
                                              theta_intron) {
  
  # Get per-AA mutation rates from Q-matrix
  aa_rates <- get_aa_mutation_rates_from_Q(Q, preferred_codons, genetic_code)
  
  # The Q-matrix rates are relative. To get alpha = 4N*u, we need to scale
  # by the total mutation rate (theta = 4N*mu_total)
  # 
  # For each AA: alpha = theta * (u / mu_total)
  #              beta  = theta * (v / mu_total)
  #
  # Since Q is normalized to 1 substitution per unit time,
  # we can use theta directly as the scaling factor
  
  aa_rates[, `:=`(
    Alpha = theta_intron * u,
    Beta = theta_intron * v
  )]
  
  # Add terminal nucleotide for reference
  aa_rates[, Terminal_Nuc := sapply(strsplit(Preferred_Codon, ","), function(x) {
    substr(x[1], 3, 3)
  })]
  
  cat("\n=== Per-AA Neutral Parameters (from Q-matrix) ===\n\n")
  cat(sprintf("Global theta (from introns): %.6f\n\n", theta_intron))
  
  print(aa_rates[, .(AA, Terminal_Nuc, u = round(u, 6), v = round(v, 6), 
                     Alpha = round(Alpha, 6), Beta = round(Beta, 6))])
  
  return(aa_rates)
}

#' Calculate alpha and beta for each amino acid from intronic data
#' 
#' SIMPLE VERSION: Uses terminal nucleotide to assign global G/C parameters.
#' For amino acids with G-ending preferred codons: use alpha_G, beta_G
#' For amino acids with C-ending preferred codons: use alpha_C, beta_C
#' 
#' @param neutral_params Output from load_and_estimate_neutral_params()
#' @param preferred_codons_df Data frame with AA and Preferred_Codon columns
#' @return Data.table with per-AA alpha and beta values
get_aa_specific_neutral_params <- function(neutral_params, preferred_codons_df) {
  
  require(data.table)
  
  if (!is.data.table(preferred_codons_df)) {
    preferred_codons_df <- as.data.table(preferred_codons_df)
  }
  
  # Handle column naming
  codon_col <- if ("Preferred_Codon" %in% names(preferred_codons_df)) {
    "Preferred_Codon"
  } else if ("Codon" %in% names(preferred_codons_df)) {
    "Codon"
  } else {
    stop("Need Preferred_Codon or Codon column")
  }
  
  # Get terminal nucleotide
  preferred_codons_df[, Terminal_Nuc := substr(get(codon_col), 3, 3)]
  
  # Assign alpha/beta based on terminal nucleotide
  aa_params <- preferred_codons_df[, .(
    AA = AA,
    Preferred_Codon = get(codon_col),
    Terminal_Nuc = Terminal_Nuc,
    Alpha = ifelse(Terminal_Nuc == "G", neutral_params$alpha_G, neutral_params$alpha_C),
    Beta = ifelse(Terminal_Nuc == "G", neutral_params$beta_G, neutral_params$beta_C)
  )]
  
  cat(sprintf("✓ Assigned neutral parameters for %d amino acids\n", nrow(aa_params)))
  cat(sprintf("  G-ending AAs: %d (α=%.4f, β=%.4f)\n",
              sum(aa_params$Terminal_Nuc == "G"),
              neutral_params$alpha_G, neutral_params$beta_G))
  cat(sprintf("  C-ending AAs: %d (α=%.4f, β=%.4f)\n",
              sum(aa_params$Terminal_Nuc == "C"),
              neutral_params$alpha_C, neutral_params$beta_C))
  
  return(aa_params)
}

#' Estimate gamma for each amino acid separately
#' 
#' This function estimates gamma (4Nes) for each amino acid using the
#' appropriate alpha/beta values for that amino acid's terminal nucleotide.
#' 
#' @param sfs_by_aa Data.table with columns: AA, k, n (SFS data per amino acid)
#' @param aa_params Output from get_aa_specific_neutral_params()
#' @param min_sites Minimum sites required for estimation (default: 20)
#' @return Data.table with gamma estimates per amino acid
estimate_gamma_per_amino_acid <- function(sfs_by_aa, aa_params, min_sites = 20) {
  
  require(data.table)
  require(gsl)
  
  if (!is.data.table(sfs_by_aa)) setDT(sfs_by_aa)
  if (!is.data.table(aa_params)) setDT(aa_params)
  
  cat("\n=== Estimating Gamma Per Amino Acid ===\n\n")
  
  # Merge SFS data with AA parameters
  setkey(sfs_by_aa, AA)
  setkey(aa_params, AA)
  
  results <- aa_params[, {
    
    aa <- AA
    alpha <- Alpha
    beta <- Beta
    
    # Get SFS for this AA
    aa_sfs <- sfs_by_aa[AA == aa]
    
    if (nrow(aa_sfs) < min_sites) {
      list(
        Gamma = NA_real_,
        SE_Gamma = NA_real_,
        N_Sites = nrow(aa_sfs),
        Status = "Too few sites"
      )
    } else {
      # Estimate gamma using MLE
      gamma_est <- tryCatch({
        estimate_gamma_for_AA(
          counts = aa_sfs$k,
          sample_sizes = aa_sfs$n,
          alpha = alpha,
          beta = beta,
          S_interval = c(-10, 50)
        )
      }, error = function(e) NA_real_)
      
      # Bootstrap SE
      if (!is.na(gamma_est)) {
        boot_gammas <- replicate(100, {
          idx <- sample(nrow(aa_sfs), replace = TRUE)
          tryCatch({
            estimate_gamma_for_AA(
              counts = aa_sfs$k[idx],
              sample_sizes = aa_sfs$n[idx],
              alpha = alpha,
              beta = beta,
              S_interval = c(-10, 50)
            )
          }, error = function(e) NA_real_)
        })
        se_gamma <- sd(boot_gammas, na.rm = TRUE)
      } else {
        se_gamma <- NA_real_
      }
      
      list(
        Gamma = gamma_est,
        SE_Gamma = se_gamma,
        N_Sites = nrow(aa_sfs),
        Status = ifelse(is.na(gamma_est), "Estimation failed", "Success")
      )
    }
  }, by = .(AA, Preferred_Codon, Terminal_Nuc, Alpha, Beta)]
  
  # Summary
  n_success <- sum(results$Status == "Success")
  cat(sprintf("Successfully estimated: %d/%d amino acids\n", 
              n_success, nrow(results)))
  
  if (n_success > 0) {
    cat("\nGamma summary (successful estimates):\n")
    cat(sprintf("  Mean: %.3f\n", mean(results$Gamma, na.rm = TRUE)))
    cat(sprintf("  Median: %.3f\n", median(results$Gamma, na.rm = TRUE)))
    cat(sprintf("  Range: [%.3f, %.3f]\n",
                min(results$Gamma, na.rm = TRUE),
                max(results$Gamma, na.rm = TRUE)))
  }
  
  return(results)
}

#' Extract expected codon frequencies from ROC model
#' 
#' Given CSP parameters (dM, dEta), calculate expected proportion of
#' preferred codon at different expression levels.
#' 
#' ROC model: P(preferred | phi) = exp(-dM_pref - dEta_pref * phi) / Z
#' 
#' @param csp_df CSP parameters with columns: AA, Codon, dM, dEta, is_optimal
#' @param phi_values Vector of expression values (linear scale)
#' @return Data.table with expected P(preferred) for each AA at each phi
get_roc_expected_frequencies <- function(csp_df, phi_values) {
  
  require(data.table)
  
  if (!is.data.table(csp_df)) setDT(csp_df)
  
  results <- lapply(unique(csp_df$AA), function(aa) {
    
    aa_df <- csp_df[AA == aa]
    
    # For each phi, calculate P(preferred)
    probs <- sapply(phi_values, function(phi) {
      
      # Log-unnormalized probabilities
      log_unnorm <- -aa_df$dM - aa_df$dEta * phi
      
      # Normalize
      max_log <- max(log_unnorm)
      log_Z <- max_log + log(sum(exp(log_unnorm - max_log)))
      
      probs <- exp(log_unnorm - log_Z)
      
      # Return probability of preferred codon
      sum(probs[aa_df$is_optimal])
    })
    
    data.table(
      AA = aa,
      Phi = phi_values,
      P_Preferred_Expected = probs
    )
  })
  
  rbindlist(results)
}

#' Calculate expected pi (diversity) from ROC model at given expression
#' 
#' Using the ROC-predicted allele frequency and the neutral params,
#' calculate what pi we expect at 4-fold sites for a given expression level.
#' 
#' @param p_pref Expected frequency of preferred allele from ROC model
#' @param alpha 4N*u (mutation unpreferred -> preferred)
#' @param beta 4N*v (mutation preferred -> unpreferred)
#' @param gamma Selection coefficient (4Nes)
#' @return Expected nucleotide diversity
calculate_expected_pi_from_roc <- function(p_pref, alpha, beta, gamma) {
  
  # At mutation-selection-drift equilibrium:
  # The expected heterozygosity depends on the balance
  
  # For a biallelic site under selection:
  # E[pi] = integral[2*p*(1-p) * f(p) dp] / integral[f(p) dp]
  # where f(p) is the Wright distribution
  
  # We already have this function
  calculate_pi_analytical(alpha, beta, gamma)
}

#' Test 1: Gene-level correlation between S_ROC and per-gene gamma
#' 
#' @param integrated_data Main data frame with S_ROC values
#' @param gamma_gene Data frame with per-gene gamma estimates
#' @param output_dir Directory to save plots
#' @return List with correlation results and plot
test_sroc_gamma_correlation <- function(integrated_data, gamma_gene, output_dir) {
  
  require(data.table)
  require(ggplot2)
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("ANALYSIS 1: Gene-Level Correlation Between S_ROC and Gamma\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  # Merge data
  if (!is.data.table(integrated_data)) setDT(integrated_data)
  if (!is.data.table(gamma_gene)) setDT(gamma_gene)
  
  merged <- merge(
    integrated_data[, .(Gene_name, S_ROC, Expression_Group, Max_Log10_Exp)],
    gamma_gene[, .(Gene_name = Gene_name, Gamma_Mean = Gamma_Weighted_Mean,
                   Selection_Intensity_Gamma = Selection_Intensity)],
    by = "Gene_name"
  )
  
  n_genes <- nrow(merged)
  cat(sprintf("Genes in both analyses: %d\n\n", n_genes))
  
  # Filter valid data
  merged <- merged[!is.na(S_ROC) & !is.na(Gamma_Mean) & S_ROC > 0]
  cat(sprintf("Genes with valid data: %d\n\n", nrow(merged)))
  
  # Spearman correlation
  cor_spearman <- cor.test(merged$S_ROC, merged$Gamma_Mean, 
                           method = "spearman", exact = FALSE)
  cor_kendall <- cor.test(merged$S_ROC, merged$Gamma_Mean, 
                          method = "kendall", exact = FALSE)
  
  cat("Overall Correlation:\n")
  cat(sprintf("  Spearman ρ = %.4f (p = %.2e)\n", 
              cor_spearman$estimate, cor_spearman$p.value))
  cat(sprintf("  Kendall τ = %.4f (p = %.2e)\n\n", 
              cor_kendall$estimate, cor_kendall$p.value))
  
  # Stratified by expression
  cat("Stratified Correlations:\n")
  strat_results <- merged[, {
    cor_test <- cor.test(S_ROC, Gamma_Mean, method = "spearman", exact = FALSE)
    list(
      N = .N,
      Rho = cor_test$estimate,
      P_value = cor_test$p.value,
      Mean_S_ROC = mean(S_ROC, na.rm = TRUE),
      Mean_Gamma = mean(Gamma_Mean, na.rm = TRUE)
    )
  }, by = Expression_Group]
  
  print(strat_results)
  
  # Create scatterplot
  p <- ggplot(merged, aes(x = log10(S_ROC + 0.01), y = Gamma_Mean)) +
    geom_point(aes(color = Expression_Group), alpha = 0.5, size = 1) +
    geom_smooth(method = "loess", color = "red", se = TRUE) +
    geom_vline(xintercept = log10(1), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = log10(5), linetype = "dotted", color = "gray50") +
    scale_color_manual(values = c("Top 5%" = "#E41A1C", 
                                  "Middle 90%" = "#999999",
                                  "Bottom 5%" = "#377EB8")) +
    labs(
      title = "Gene-Level Correlation: ROC Selection vs Polymorphism Gamma",
      subtitle = sprintf("Spearman ρ = %.3f (p = %.2e)", 
                        cor_spearman$estimate, cor_spearman$p.value),
      x = expression(log[10](S[ROC])),
      y = expression(gamma[polymorphism]),
      color = "Expression Group"
    ) +
    theme_bw() +
    annotate("text", x = log10(1), y = max(merged$Gamma_Mean) * 0.9,
             label = "S=1", hjust = -0.2) +
    annotate("text", x = log10(5), y = max(merged$Gamma_Mean) * 0.9,
             label = "S=5", hjust = -0.2)
  
  ggsave(file.path(output_dir, "Analysis1_SROC_vs_Gamma_correlation.pdf"),
         p, width = 10, height = 8)
  
  cat(sprintf("\n✓ Plot saved: %s\n", 
              file.path(output_dir, "Analysis1_SROC_vs_Gamma_correlation.pdf")))
  
  # Summary table
  summary_by_group <- merged[, .(
    N = .N,
    Mean_S_ROC = mean(S_ROC, na.rm = TRUE),
    SD_S_ROC = sd(S_ROC, na.rm = TRUE),
    Mean_Gamma = mean(Gamma_Mean, na.rm = TRUE),
    SD_Gamma = sd(Gamma_Mean, na.rm = TRUE)
  ), by = Expression_Group]
  
  cat("\n=== Summary by Expression Category ===\n")
  print(summary_by_group)
  
  # Interpretation
  cat("\n=== INTERPRETATION ===\n")
  if (cor_spearman$estimate > 0.3 && cor_spearman$p.value < 0.01) {
    cat("✓ CONSISTENT: Positive correlation between S_ROC and gamma.\n")
    cat("  Both methods identify similar genes under strong selection.\n")
  } else if (abs(cor_spearman$estimate) < 0.1) {
    cat("✗ INCONSISTENT: No correlation between S_ROC and gamma.\n")
    cat("  Methods may be measuring different processes.\n")
  } else {
    cat("⚠ PARTIALLY CONSISTENT: Weak correlation.\n")
    cat("  Some agreement but systematic differences exist.\n")
  }
  
  return(list(
    correlation = cor_spearman,
    correlation_kendall = cor_kendall,
    stratified = strat_results,
    summary = summary_by_group,
    merged_data = merged,
    plot = p
  ))
}

#' Test 2: Expression-stratified gamma comparison
#' 
#' @param vcf_data Codon VCF data with k, n, Gene, AA columns
#' @param integrated_data Main data with Expression_Group
#' @param aa_params Per-AA neutral parameters
#' @param output_dir Directory to save plots
#' @return List with stratified gamma estimates and tests
test_expression_stratified_gamma <- function(vcf_data, integrated_data,
                                             aa_params, output_dir) {
  
  require(data.table)
  require(ggplot2)
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("ANALYSIS 2: Expression-Stratified Gamma Comparison\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  if (!is.data.table(vcf_data)) setDT(vcf_data)
  if (!is.data.table(integrated_data)) setDT(integrated_data)
  
  # Map genes to expression groups
  gene_groups <- integrated_data[, .(Gene_name, Expression_Group)]
  
  # Merge
  vcf_data[, Gene_name := paste0("MgIM767.", Gene)]
  vcf_with_group <- merge(vcf_data, gene_groups, by = "Gene_name")
  
  cat("Sites per expression group:\n")
  print(vcf_with_group[, .N, by = Expression_Group])
  cat("\n")
  
  # Estimate gamma for each group x AA combination
  results_list <- list()
  
  for (grp in c("Top 5%", "Middle 90%", "Bottom 5%")) {
    
    cat(sprintf("\nEstimating gamma for %s genes...\n", grp))
    
    grp_data <- vcf_with_group[Expression_Group == grp]
    
    # Aggregate by AA
    aa_sfs <- grp_data[, .(k = sum(k), n = sum(n)), by = AA]
    
    # Estimate gamma per AA
    grp_gamma <- aa_params[, {
      aa_dt <- aa_sfs[AA == AA]
      
      if (nrow(aa_dt) == 0 || aa_dt$n < 100) {
        list(Gamma = NA_real_, N_Sites = 0)
      } else {
        gamma_est <- tryCatch({
          estimate_gamma_for_AA(
            counts = aa_dt$k,
            sample_sizes = aa_dt$n,
            alpha = Alpha,
            beta = Beta,
            S_interval = c(-10, 50)
          )
        }, error = function(e) NA_real_)
        
        list(Gamma = gamma_est, N_Sites = aa_dt$n)
      }
    }, by = .(AA, Terminal_Nuc, Alpha, Beta)]
    
    grp_gamma[, Expression_Group := grp]
    results_list[[grp]] <- grp_gamma
  }
  
  all_results <- rbindlist(results_list)
  
  # Summary
  summary_results <- all_results[!is.na(Gamma), .(
    N_AA = .N,
    Mean_Gamma = mean(Gamma),
    Median_Gamma = median(Gamma),
    SE_Gamma = sd(Gamma) / sqrt(.N),
    CI_Lower = mean(Gamma) - 1.96 * sd(Gamma) / sqrt(.N),
    CI_Upper = mean(Gamma) + 1.96 * sd(Gamma) / sqrt(.N)
  ), by = .(Expression_Group, Terminal_Nuc)]
  
  cat("\n=== Gamma Estimates by Expression Group and Terminal Nucleotide ===\n")
  print(summary_results)
  
  # Statistical tests
  cat("\n=== Statistical Tests ===\n")
  
  # Wilcoxon test: High vs Low
  high_gamma <- all_results[Expression_Group == "Top 5%" & !is.na(Gamma)]$Gamma
  low_gamma <- all_results[Expression_Group == "Bottom 5%" & !is.na(Gamma)]$Gamma
  
  if (length(high_gamma) > 3 && length(low_gamma) > 3) {
    wilcox_hl <- wilcox.test(high_gamma, low_gamma)
    cat(sprintf("High vs Low expression (Wilcoxon): W = %.0f, p = %.4f\n",
                wilcox_hl$statistic, wilcox_hl$p.value))
    
    # Effect size
    d <- (mean(high_gamma) - mean(low_gamma)) / 
      sqrt((var(high_gamma) + var(low_gamma)) / 2)
    cat(sprintf("Effect size (Cohen's d): %.3f\n", d))
  }
  
  # Create visualization
  p <- ggplot(summary_results, 
              aes(x = Expression_Group, y = Mean_Gamma, 
                  fill = Terminal_Nuc)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper),
                  position = position_dodge(width = 0.8),
                  width = 0.2) +
    scale_fill_manual(values = c("G" = "#E41A1C", "C" = "#377EB8"),
                      name = "Terminal\nNucleotide") +
    labs(
      title = "Selection Coefficient (γ) by Expression Level",
      subtitle = "Error bars = 95% CI",
      x = "Expression Group",
      y = expression(gamma ~ "(4N"[e]*"s)")
    ) +
    theme_bw() +
    theme(legend.position = "right")
  
  ggsave(file.path(output_dir, "Analysis2_Gamma_by_expression_group.pdf"),
         p, width = 8, height = 6)
  
  cat(sprintf("\n✓ Plot saved: %s\n",
              file.path(output_dir, "Analysis2_Gamma_by_expression_group.pdf")))
  
  # Interpretation
  cat("\n=== INTERPRETATION ===\n")
  
  gamma_high <- summary_results[Expression_Group == "Top 5%", Mean_Gamma]
  gamma_low <- summary_results[Expression_Group == "Bottom 5%", Mean_Gamma]
  
  if (length(gamma_high) > 0 && length(gamma_low) > 0) {
    if (mean(gamma_high) > mean(gamma_low) * 1.5) {
      cat("✓ CONSISTENT: γ is higher in high-expression genes.\n")
      cat("  Polymorphism data shows expression-dependent selection.\n")
    } else if (abs(mean(gamma_high) - mean(gamma_low)) < 0.3) {
      cat("✗ INCONSISTENT: γ is uniform across expression levels.\n")
      cat("  No evidence for stronger selection in highly expressed genes.\n")
    } else {
      cat("⚠ PARTIAL: Some difference in γ but not dramatic.\n")
    }
  }
  
  return(list(
    results = all_results,
    summary = summary_results,
    plot = p
  ))
}

#' Test 3: Diversity hump magnitude vs expression
#' 
#' Calculate Δπ = π_4fold - π_intron for each gene and test
#' if it correlates with expression.
#' 
#' @param pi_by_gene Data frame with pi values per gene by site type
#' @param integrated_data Main data with expression values
#' @param output_dir Directory to save plots
#' @return List with correlation results
test_diversity_hump_vs_expression <- function(pi_by_gene, integrated_data, output_dir) {
  
  require(data.table)
  require(ggplot2)
  require(mgcv)
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("ANALYSIS 3: Diversity Hump Magnitude vs Expression\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  if (!is.data.table(pi_by_gene)) setDT(pi_by_gene)
  if (!is.data.table(integrated_data)) setDT(integrated_data)
  
  # Merge pi data with expression
  merged <- merge(
    pi_by_gene,
    integrated_data[, .(Gene_name, Max_Log10_Exp, Expression_Group)],
    by = "Gene_name"
  )
  
  # Calculate diversity boost
  merged[, Delta_Pi := Pi_4fold_GC - Pi_Intron_GC]
  
  cat(sprintf("Genes with diversity data: %d\n", nrow(merged)))
  cat(sprintf("Mean Δπ: %.6f\n", mean(merged$Delta_Pi, na.rm = TRUE)))
  cat("\n")
  
  # Correlation tests
  cor_spearman <- cor.test(merged$Max_Log10_Exp, merged$Delta_Pi,
                           method = "spearman", exact = FALSE)
  cor_kendall <- cor.test(merged$Max_Log10_Exp, merged$Delta_Pi,
                          method = "kendall", exact = FALSE)
  
  cat("Correlation: Expression vs Δπ\n")
  cat(sprintf("  Spearman ρ = %.4f (p = %.2e)\n",
              cor_spearman$estimate, cor_spearman$p.value))
  cat(sprintf("  Kendall τ = %.4f (p = %.2e)\n\n",
              cor_kendall$estimate, cor_kendall$p.value))
  
  # GAM model
  gam_model <- gam(Delta_Pi ~ s(Max_Log10_Exp) + s(log10(Gene_Length)),
                   data = merged[!is.na(Delta_Pi) & Gene_Length > 0])
  
  cat("GAM Model Summary:\n")
  cat(sprintf("  Deviance explained: %.1f%%\n", 
              100 * summary(gam_model)$dev.expl))
  cat(sprintf("  Expression effect p-value: %.2e\n",
              summary(gam_model)$s.pv[1]))
  
  # Group comparison
  group_summary <- merged[!is.na(Delta_Pi), .(
    N = .N,
    Mean_Delta_Pi = mean(Delta_Pi),
    SD_Delta_Pi = sd(Delta_Pi),
    Median_Delta_Pi = median(Delta_Pi)
  ), by = Expression_Group]
  
  cat("\n=== Δπ by Expression Group ===\n")
  print(group_summary)
  
  # Kruskal-Wallis test
  kw_test <- kruskal.test(Delta_Pi ~ Expression_Group, 
                          data = merged[!is.na(Delta_Pi)])
  cat(sprintf("\nKruskal-Wallis test: χ² = %.2f, p = %.4f\n",
              kw_test$statistic, kw_test$p.value))
  
  # Create plots
  p1 <- ggplot(merged[!is.na(Delta_Pi)], 
               aes(x = Max_Log10_Exp, y = Delta_Pi)) +
    geom_point(aes(color = Expression_Group), alpha = 0.3, size = 0.5) +
    geom_smooth(method = "gam", formula = y ~ s(x), color = "red") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("Top 5%" = "#E41A1C",
                                  "Middle 90%" = "#999999",
                                  "Bottom 5%" = "#377EB8")) +
    labs(
      title = "Diversity Boost vs Expression",
      subtitle = sprintf("Spearman ρ = %.3f", cor_spearman$estimate),
      x = expression(log[10](Expression)),
      y = expression(Delta*pi ~ "(4-fold - intron)")
    ) +
    theme_bw()
  
  p2 <- ggplot(merged[!is.na(Delta_Pi)],
               aes(x = Expression_Group, y = Delta_Pi, fill = Expression_Group)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.2, outlier.alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("Top 5%" = "#E41A1C",
                                 "Middle 90%" = "#999999",
                                 "Bottom 5%" = "#377EB8")) +
    labs(
      title = "Diversity Boost by Expression Group",
      y = expression(Delta*pi)
    ) +
    theme_bw() +
    theme(legend.position = "none")
  
  combined <- cowplot::plot_grid(p1, p2, ncol = 2, rel_widths = c(1.2, 1))
  
  ggsave(file.path(output_dir, "Analysis3_Diversity_hump_vs_expression.pdf"),
         combined, width = 12, height = 5)
  
  cat(sprintf("\n✓ Plot saved: %s\n",
              file.path(output_dir, "Analysis3_Diversity_hump_vs_expression.pdf")))
  
  # Interpretation
  cat("\n=== INTERPRETATION ===\n")
  if (cor_spearman$estimate > 0.1 && cor_spearman$p.value < 0.05) {
    cat("✓ Δπ increases with expression.\n")
    cat("  Supports expression-dependent selection on polymorphisms.\n")
  } else if (abs(cor_spearman$estimate) < 0.05) {
    cat("✗ Δπ does not vary with expression.\n")
    cat("  Selection on polymorphisms appears uniform.\n")
  }
  
  return(list(
    correlation = cor_spearman,
    gam_model = gam_model,
    group_summary = group_summary,
    kw_test = kw_test,
    plots = combined
  ))
}

#' Test 4: S_ROC-stratified gamma estimation
#' 
#' Direct test: partition genes by S_ROC and estimate gamma for each partition.
#' 
#' @param vcf_data Codon VCF data
#' @param integrated_data Main data with S_ROC values
#' @param aa_params Per-AA neutral parameters
#' @param output_dir Directory to save plots
#' @return List with results
test_sroc_stratified_gamma <- function(vcf_data, integrated_data, 
                                       aa_params, output_dir) {
  
  require(data.table)
  require(ggplot2)
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("ANALYSIS 4: S_ROC-Stratified Gamma Estimation (MOST DIRECT TEST)\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  if (!is.data.table(vcf_data)) setDT(vcf_data)
  if (!is.data.table(integrated_data)) setDT(integrated_data)
  
  # Create S_ROC categories
  integrated_data[, SROC_Category := cut(
    S_ROC,
    breaks = c(-Inf, 0.1, 1, 5, Inf),
    labels = c("S<0.1 (Drift)", "0.1≤S<1 (Weak)", 
               "1≤S<5 (Moderate)", "S≥5 (Strong)"),
    right = FALSE
  )]
  
  cat("Gene distribution by S_ROC category:\n")
  print(table(integrated_data$SROC_Category))
  cat("\n")
  
  # Map genes to categories
  gene_sroc <- integrated_data[!is.na(SROC_Category), 
                                .(Gene_name, SROC_Category, S_ROC)]
  
  # Merge with VCF data
  vcf_data[, Gene_name := paste0("MgIM767.", Gene)]
  vcf_with_sroc <- merge(vcf_data, gene_sroc, by = "Gene_name")
  
  # Estimate gamma for each S_ROC category
  results_list <- list()
  
  for (cat_name in levels(vcf_with_sroc$SROC_Category)) {
    
    cat(sprintf("\nEstimating gamma for %s genes...\n", cat_name))
    
    cat_data <- vcf_with_sroc[SROC_Category == cat_name]
    
    if (nrow(cat_data) < 100) {
      cat("  Insufficient sites, skipping.\n")
      next
    }
    
    # Aggregate SFS by terminal nucleotide
    for (tn in c("G", "C")) {
      
      aa_subset <- aa_params[Terminal_Nuc == tn]
      cat_aa_data <- cat_data[AA %in% aa_subset$AA]
      
      if (nrow(cat_aa_data) < 50) next
      
      # Pool sites
      pooled_k <- sum(cat_aa_data$k)
      pooled_n <- sum(cat_aa_data$n)
      
      # Get average alpha/beta for this terminal nuc
      alpha_use <- mean(aa_subset$Alpha)
      beta_use <- mean(aa_subset$Beta)
      
      # Estimate gamma
      gamma_est <- tryCatch({
        estimate_gamma_for_AA(
          counts = cat_aa_data$k,
          sample_sizes = cat_aa_data$n,
          alpha = alpha_use,
          beta = beta_use,
          S_interval = c(-10, 50)
        )
      }, error = function(e) NA_real_)
      
      # Bootstrap CI
      boot_gammas <- replicate(500, {
        idx <- sample(nrow(cat_aa_data), replace = TRUE)
        tryCatch({
          estimate_gamma_for_AA(
            counts = cat_aa_data$k[idx],
            sample_sizes = cat_aa_data$n[idx],
            alpha = alpha_use,
            beta = beta_use,
            S_interval = c(-10, 50)
          )
        }, error = function(e) NA_real_)
      })
      
      ci <- quantile(boot_gammas, c(0.025, 0.975), na.rm = TRUE)
      
      results_list[[paste(cat_name, tn, sep = "_")]] <- data.table(
        SROC_Category = cat_name,
        Terminal_Nuc = tn,
        N_Sites = nrow(cat_aa_data),
        N_Genes = length(unique(cat_aa_data$Gene_name)),
        Mean_SROC = mean(cat_data$S_ROC),
        Gamma = gamma_est,
        CI_Lower = ci[1],
        CI_Upper = ci[2]
      )
    }
  }
  
  results <- rbindlist(results_list)
  
  cat("\n=== Gamma Estimates by S_ROC Category ===\n")
  print(results)
  
  # Trend test (Jonckheere-Terpstra)
  cat("\n=== Trend Test ===\n")
  
  # Simple linear regression of gamma on mean S_ROC
  if (nrow(results[!is.na(Gamma)]) >= 3) {
    lm_trend <- lm(Gamma ~ Mean_SROC, data = results)
    cat(sprintf("Linear trend: slope = %.4f, p = %.4f\n",
                coef(lm_trend)[2],
                summary(lm_trend)$coefficients[2, 4]))
    
    # Correlation
    cor_trend <- cor.test(results$Mean_SROC, results$Gamma, 
                          method = "spearman", exact = FALSE)
    cat(sprintf("Spearman correlation: ρ = %.3f, p = %.4f\n",
                cor_trend$estimate, cor_trend$p.value))
  }
  
  # Create visualization
  results[, SROC_Order := as.numeric(factor(SROC_Category,
                                            levels = c("S<0.1 (Drift)", "0.1≤S<1 (Weak)",
                                                       "1≤S<5 (Moderate)", "S≥5 (Strong)")))]
  
  p <- ggplot(results[!is.na(Gamma)], 
              aes(x = reorder(SROC_Category, SROC_Order), y = Gamma)) +
    geom_col(aes(fill = Terminal_Nuc), position = position_dodge(width = 0.8),
             width = 0.7) +
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper, group = Terminal_Nuc),
                  position = position_dodge(width = 0.8),
                  width = 0.2) +
    scale_fill_manual(values = c("G" = "#E41A1C", "C" = "#377EB8"),
                      name = "Terminal\nNucleotide") +
    labs(
      title = "Selection Coefficient (γ) by ROC Selection Intensity",
      subtitle = "If consistent: γ should increase with S_ROC category",
      x = expression(S[ROC] ~ "Category"),
      y = expression(gamma ~ "from Polymorphism (4N"[e]*"s)")
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 15, hjust = 1))
  
  ggsave(file.path(output_dir, "Analysis4_Gamma_by_SROC_category.pdf"),
         p, width = 10, height = 6)
  
  cat(sprintf("\n✓ Plot saved: %s\n",
              file.path(output_dir, "Analysis4_Gamma_by_SROC_category.pdf")))
  
  # Interpretation
  cat("\n=== INTERPRETATION ===\n")
  
  if (nrow(results[!is.na(Gamma)]) >= 3) {
    gamma_trend <- results[!is.na(Gamma)]
    gamma_low <- gamma_trend[SROC_Category == "S<0.1 (Drift)", mean(Gamma)]
    gamma_high <- gamma_trend[SROC_Category == "S≥5 (Strong)", mean(Gamma)]
    
    if (!is.na(gamma_low) && !is.na(gamma_high)) {
      if (gamma_high > gamma_low * 1.5) {
        cat("✓ CONSISTENT: γ increases with S_ROC category.\n")
        cat("  Genes with high ROC selection also show elevated γ from polymorphism.\n")
      } else if (abs(gamma_high - gamma_low) < 0.5) {
        cat("✗ INCONSISTENT: γ is similar across S_ROC categories.\n")
        cat("  ROC selection intensity does not predict polymorphism-based selection.\n")
      } else {
        cat("⚠ PARTIAL: Some trend but not as strong as expected.\n")
      }
    }
  }
  
  return(list(
    results = results,
    plot = p
  ))
}

#' Test 6: Low-expression genes vs intronic baseline
#' 
#' Validate that the mutation baseline (alpha, beta) from introns
#' matches the nucleotide composition of low-expression genes.
#' 
#' @param vcf_data Codon VCF data
#' @param integrated_data Main data
#' @param neutral_params Neutral parameters from introns
#' @param output_dir Directory to save plots
#' @return List with validation results
test_low_expression_baseline <- function(vcf_data, integrated_data,
                                         neutral_params, output_dir) {
  
  require(data.table)
  require(ggplot2)
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("ANALYSIS 6: Low-Expression Genes vs Intronic Baseline\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  if (!is.data.table(vcf_data)) setDT(vcf_data)
  if (!is.data.table(integrated_data)) setDT(integrated_data)
  
  # Get bottom 5% expression genes
  low_exp_genes <- integrated_data[Expression_Group == "Bottom 5%", Gene_name]
  
  cat(sprintf("Low-expression genes: %d\n\n", length(low_exp_genes)))
  
  # Get VCF data for low-expression genes
  vcf_data[, Gene_name := paste0("MgIM767.", Gene)]
  low_exp_vcf <- vcf_data[Gene_name %in% low_exp_genes]
  
  # Calculate nucleotide frequencies at 4-fold sites
  # p = k/n is the frequency of the preferred codon
  low_exp_vcf[, p := k / n]
  
  # Separate by terminal nucleotide
  # For G-ending preferred: p = freq of G
  # For C-ending preferred: p = freq of C
  
  mean_freq_G <- low_exp_vcf[Terminal_Nuc == "G", mean(p, na.rm = TRUE)]
  mean_freq_C <- low_exp_vcf[Terminal_Nuc == "C", mean(p, na.rm = TRUE)]
  
  # Expected frequency under neutrality (mutation-drift equilibrium)
  # E[p] = alpha / (alpha + beta)
  expected_freq_G <- neutral_params$alpha_G / 
    (neutral_params$alpha_G + neutral_params$beta_G)
  expected_freq_C <- neutral_params$alpha_C / 
    (neutral_params$alpha_C + neutral_params$beta_C)
  
  cat("=== Nucleotide Frequency Comparison ===\n\n")
  cat("G-ending sites:\n")
  cat(sprintf("  Observed in low-exp genes: %.4f\n", mean_freq_G))
  cat(sprintf("  Expected from introns:     %.4f\n", expected_freq_G))
  cat(sprintf("  Difference: %.4f\n\n", mean_freq_G - expected_freq_G))
  
  cat("C-ending sites:\n")
  cat(sprintf("  Observed in low-exp genes: %.4f\n", mean_freq_C))
  cat(sprintf("  Expected from introns:     %.4f\n", expected_freq_C))
  cat(sprintf("  Difference: %.4f\n\n", mean_freq_C - expected_freq_C))
  
  # Statistical test: Is the difference significant?
  # Use one-sample t-test comparing observed to expected
  
  g_sites <- low_exp_vcf[Terminal_Nuc == "G" & !is.na(p), p]
  c_sites <- low_exp_vcf[Terminal_Nuc == "C" & !is.na(p), p]
  
  if (length(g_sites) > 30) {
    t_test_G <- t.test(g_sites, mu = expected_freq_G)
    cat(sprintf("T-test (G sites): t = %.2f, p = %.4f\n",
                t_test_G$statistic, t_test_G$p.value))
  }
  
  if (length(c_sites) > 30) {
    t_test_C <- t.test(c_sites, mu = expected_freq_C)
    cat(sprintf("T-test (C sites): t = %.2f, p = %.4f\n",
                t_test_C$statistic, t_test_C$p.value))
  }
  
  # Create visualization
  comparison_df <- data.table(
    Nucleotide = rep(c("G", "C"), 2),
    Source = rep(c("Low-Exp Genes", "Intronic (Expected)"), each = 2),
    Frequency = c(mean_freq_G, mean_freq_C, expected_freq_G, expected_freq_C)
  )
  
  p <- ggplot(comparison_df, aes(x = Nucleotide, y = Frequency, fill = Source)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("Low-Exp Genes" = "#69b3a2",
                                 "Intronic (Expected)" = "#404080")) +
    labs(
      title = "Validation: Low-Expression Genes vs Intronic Baseline",
      subtitle = "If baseline is correct, bars should match",
      y = "Mean Frequency of Preferred Codon",
      x = "Preferred Codon Terminal Nucleotide"
    ) +
    theme_bw() +
    ylim(0, 1)
  
  ggsave(file.path(output_dir, "Analysis6_Low_expression_baseline_check.pdf"),
         p, width = 8, height = 6)
  
  cat(sprintf("\n✓ Plot saved: %s\n",
              file.path(output_dir, "Analysis6_Low_expression_baseline_check.pdf")))
  
  # Interpretation
  cat("\n=== INTERPRETATION ===\n")
  
  diff_G <- abs(mean_freq_G - expected_freq_G)
  diff_C <- abs(mean_freq_C - expected_freq_C)
  
  if (diff_G < 0.02 && diff_C < 0.02) {
    cat("✓ BASELINE VALID: Low-expression genes match intronic expectations.\n")
    cat("  The mutation parameters from introns are appropriate.\n")
  } else if (diff_G > 0.05 || diff_C > 0.05) {
    cat("✗ BASELINE MISMATCH: Large difference from intronic expectations.\n")
    cat("  This may explain inconsistencies between ROC and polymorphism analysis.\n")
    cat("  Consider:\n")
    cat("    1. Using low-expression gene data to re-estimate alpha/beta\n")
    cat("    2. Checking for selection even in 'low-expression' genes\n")
  } else {
    cat("⚠ PARTIAL MATCH: Some deviation from expectations.\n")
    cat("  Results should be interpreted with caution.\n")
  }
  
  return(list(
    observed = c(G = mean_freq_G, C = mean_freq_C),
    expected = c(G = expected_freq_G, C = expected_freq_C),
    comparison = comparison_df,
    plot = p
  ))
}

#' ROC Model Goodness of Fit Test
#' 
#' Compare observed codon frequencies per amino acid to ROC model predictions.
#' Uses chi-square or likelihood ratio test.
#' 
#' @param codon_freq_obs Observed codon frequencies per gene per AA
#' @param csp_df CSP parameters from ROC model
#' @param phi_values Expression values (MeanPhi) per gene
#' @return Data.table with GOF statistics per gene
test_roc_goodness_of_fit <- function(codon_freq_obs, csp_df, phi_values) {
  
  require(data.table)
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("ROC Model Goodness of Fit Test\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  if (!is.data.table(codon_freq_obs)) setDT(codon_freq_obs)
  if (!is.data.table(csp_df)) setDT(csp_df)
  if (!is.data.table(phi_values)) setDT(phi_values)
  
  # For each gene, calculate expected frequencies from ROC model
  # and compare to observed
  
  results <- codon_freq_obs[, {
    
    gene <- Gene[1]
    phi <- phi_values[Gene_name == gene, MeanPhi]
    
    if (length(phi) == 0 || is.na(phi)) {
      list(
        Chi_Sq = NA_real_,
        DF = NA_integer_,
        P_Value = NA_real_,
        G_Stat = NA_real_
      )
    } else {
      
      # Calculate expected counts for each codon
      expected_counts <- sapply(unique(AA), function(aa_code) {
        
        aa_csp <- csp_df[AA == aa_code]
        if (nrow(aa_csp) == 0) return(NULL)
        
        # Get observed counts for this AA
        aa_obs <- .SD[AA == aa_code]
        total_count <- sum(aa_obs$Count)
        
        if (total_count == 0) return(NULL)
        
        # Calculate expected probabilities from ROC model
        log_unnorm <- -aa_csp$dM - aa_csp$dEta * phi
        max_log <- max(log_unnorm)
        log_Z <- max_log + log(sum(exp(log_unnorm - max_log)))
        expected_probs <- exp(log_unnorm - log_Z)
        
        # Expected counts
        expected <- expected_probs * total_count
        
        # Observed counts
        observed <- aa_obs$Count
        
        list(observed = observed, expected = expected)
      }, simplify = FALSE)
      
      # Remove NULL entries
      expected_counts <- expected_counts[!sapply(expected_counts, is.null)]
      
      if (length(expected_counts) < 3) {
        list(
          Chi_Sq = NA_real_,
          DF = NA_integer_,
          P_Value = NA_real_,
          G_Stat = NA_real_
        )
      } else {
        # Pool all codons for chi-square test
        all_obs <- unlist(lapply(expected_counts, `[[`, "observed"))
        all_exp <- unlist(lapply(expected_counts, `[[`, "expected"))
        
        # Remove zero expected
        keep <- all_exp > 0.5
        all_obs <- all_obs[keep]
        all_exp <- all_exp[keep]
        
        # Chi-square statistic
        chi_sq <- sum((all_obs - all_exp)^2 / all_exp)
        df <- length(all_obs) - 1
        p_val <- pchisq(chi_sq, df, lower.tail = FALSE)
        
        # G-statistic (likelihood ratio)
        g_stat <- 2 * sum(all_obs * log(all_obs / all_exp), na.rm = TRUE)
        
        list(
          Chi_Sq = chi_sq,
          DF = df,
          P_Value = p_val,
          G_Stat = g_stat
        )
      }
    }
  }, by = Gene]
  
  # Summary
  cat(sprintf("Genes tested: %d\n", nrow(results[!is.na(Chi_Sq)])))
  cat(sprintf("Genes with significant deviation (p < 0.05): %d (%.1f%%)\n",
              sum(results$P_Value < 0.05, na.rm = TRUE),
              100 * mean(results$P_Value < 0.05, na.rm = TRUE)))
  
  cat("\nChi-Square Summary:\n")
  print(summary(results$Chi_Sq))
  
  return(results)
}

#' Enhanced ROC Model GOF Testing Per Amino Acid
#' 
#' This function compares observed vs expected codon frequencies at the 
#' amino acid level for binned expression levels. It provides both visual
#' diagnostics and statistical tests.
#' 
#' @param codon_freq_df Data.frame with Gene, AA, Codon, Count columns
#' @param csp_df CSP parameters from ROC model (AA, Codon, dM, dEta, is_optimal)
#' @param expression_df Data.frame with Gene, Exp_log10 columns
#' @param n_bins Number of expression bins (default: 10)
#' @param output_dir Directory for plots
#' @return List with GOF statistics and plots
test_roc_gof_per_amino_acid <- function(codon_freq_df, csp_df, expression_df,
                                         n_bins = 10, output_dir = NULL) {
  
  require(data.table)
  require(ggplot2)
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("ROC Model Goodness-of-Fit: Per Amino Acid Analysis\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  if (!is.data.table(codon_freq_df)) setDT(codon_freq_df)
  if (!is.data.table(csp_df)) setDT(csp_df)
  if (!is.data.table(expression_df)) setDT(expression_df)
  
  # Merge expression with codon frequencies
  merged <- merge(codon_freq_df, expression_df, by = "Gene")
  
  # Create expression bins
  merged[, Exp_Bin := cut(Exp_log10, 
                          breaks = quantile(Exp_log10, probs = seq(0, 1, length.out = n_bins + 1),
                                            na.rm = TRUE),
                          include.lowest = TRUE,
                          labels = FALSE)]
  
  # Calculate bin midpoints for prediction
  bin_midpoints <- merged[, .(Bin_Midpoint = mean(Exp_log10)), by = Exp_Bin]
  setkey(bin_midpoints, Exp_Bin)
  
  # Results storage
  gof_results <- list()
  
  # Process each amino acid
  for (aa in unique(csp_df$AA)) {
    
    aa_csp <- csp_df[AA == aa]
    aa_obs <- merged[AA == aa]
    
    if (nrow(aa_obs) < 100) next
    
    # Calculate observed frequencies by bin
    obs_by_bin <- aa_obs[, {
      total_count <- sum(Count)
      freq <- Count / total_count
      list(Codon = Codon, Observed = freq, Total_Count = total_count)
    }, by = .(Exp_Bin)]
    
    # Calculate expected frequencies for each bin
    exp_by_bin <- lapply(unique(obs_by_bin$Exp_Bin), function(bin) {
      
      phi <- bin_midpoints[Exp_Bin == bin, Bin_Midpoint]
      if (length(phi) == 0 || is.na(phi)) return(NULL)
      
      # Convert from log10 to linear scale for ROC model
      phi_linear <- 10^phi
      
      # ROC model prediction
      log_unnorm <- -aa_csp$dM - aa_csp$dEta * phi_linear
      max_log <- max(log_unnorm)
      log_Z <- max_log + log(sum(exp(log_unnorm - max_log)))
      expected <- exp(log_unnorm - log_Z)
      
      data.table(
        Exp_Bin = bin,
        Codon = aa_csp$Codon,
        Expected = expected
      )
    })
    exp_by_bin <- rbindlist(exp_by_bin[!sapply(exp_by_bin, is.null)])
    
    # Merge observed and expected
    comparison <- merge(obs_by_bin, exp_by_bin, by = c("Exp_Bin", "Codon"))
    
    if (nrow(comparison) < 10) next
    
    # Calculate chi-square per bin
    chi_sq_by_bin <- comparison[, {
      # Use total count to weight
      total <- Total_Count[1]
      exp_count <- Expected * total
      obs_count <- Observed * total
      
      # Chi-square contribution
      chi_sq <- sum((obs_count - exp_count)^2 / pmax(exp_count, 0.5))
      df <- .N - 1
      p_val <- pchisq(chi_sq, df, lower.tail = FALSE)
      
      list(Chi_Sq = chi_sq, DF = df, P_Value = p_val, N_Codons = .N)
    }, by = Exp_Bin]
    
    # Overall chi-square for this AA
    overall_chi_sq <- sum(chi_sq_by_bin$Chi_Sq)
    overall_df <- sum(chi_sq_by_bin$DF)
    overall_p <- pchisq(overall_chi_sq, overall_df, lower.tail = FALSE)
    
    # RMSE between observed and expected
    rmse <- sqrt(mean((comparison$Observed - comparison$Expected)^2))
    
    gof_results[[aa]] <- list(
      AA = aa,
      Overall_Chi_Sq = overall_chi_sq,
      Overall_DF = overall_df,
      Overall_P_Value = overall_p,
      RMSE = rmse,
      Chi_Sq_By_Bin = chi_sq_by_bin,
      Comparison = comparison
    )
    
    cat(sprintf("  %s: χ² = %.1f (df=%d, p=%.2e), RMSE = %.4f\n",
                aa, overall_chi_sq, overall_df, overall_p, rmse))
  }
  
  # Summary statistics
  summary_df <- data.table(
    AA = names(gof_results),
    Chi_Sq = sapply(gof_results, `[[`, "Overall_Chi_Sq"),
    DF = sapply(gof_results, `[[`, "Overall_DF"),
    P_Value = sapply(gof_results, `[[`, "Overall_P_Value"),
    RMSE = sapply(gof_results, `[[`, "RMSE")
  )
  
  cat("\n=== Summary ===\n")
  cat(sprintf("Amino acids tested: %d\n", nrow(summary_df)))
  cat(sprintf("Significant deviations (p < 0.05): %d (%.1f%%)\n",
              sum(summary_df$P_Value < 0.05),
              100 * mean(summary_df$P_Value < 0.05)))
  cat(sprintf("Mean RMSE: %.4f\n", mean(summary_df$RMSE)))
  cat(sprintf("Median RMSE: %.4f\n", median(summary_df$RMSE)))
  
  # Create diagnostic plot if output_dir provided
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    # RMSE by AA plot
    p1 <- ggplot(summary_df, aes(x = reorder(AA, RMSE), y = RMSE)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_hline(yintercept = median(summary_df$RMSE), 
                 linetype = "dashed", color = "red") +
      coord_flip() +
      labs(
        title = "ROC Model Goodness of Fit by Amino Acid",
        subtitle = "RMSE between observed and expected codon frequencies",
        x = "Amino Acid",
        y = "RMSE"
      ) +
      theme_minimal()
    
    ggsave(file.path(output_dir, "roc_gof_rmse_by_aa.pdf"), p1, 
           width = 8, height = 10)
    
    # P-value distribution
    p2 <- ggplot(summary_df, aes(x = P_Value)) +
      geom_histogram(bins = 20, fill = "steelblue", color = "white") +
      geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
      scale_x_continuous(limits = c(0, 1)) +
      labs(
        title = "GOF P-Value Distribution Across Amino Acids",
        subtitle = "Red line = α = 0.05",
        x = "P-Value (Chi-Square Test)",
        y = "Count"
      ) +
      theme_minimal()
    
    ggsave(file.path(output_dir, "roc_gof_pvalue_distribution.pdf"), p2,
           width = 8, height = 6)
    
    cat(sprintf("\nPlots saved to: %s\n", output_dir))
  }
  
  return(list(
    summary = summary_df,
    details = gof_results
  ))
}

#' Compare ROC Expected Frequencies to Observed Polymorphism Patterns
#' 
#' This is the KEY test: if ROC model predicts P(preferred) = 0.8 for
#' high-expression genes, do we see the same from polymorphism data?
#' 
#' @param csp_df CSP parameters
#' @param polymorphism_df Polymorphism data with allele frequencies
#' @param expression_df Expression data
#' @param neutral_params Neutral parameters (alpha, beta)
#' @return Comparison results
compare_roc_to_polymorphism <- function(csp_df, polymorphism_df, expression_df,
                                         neutral_params) {
  
  require(data.table)
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("ROC vs Polymorphism: Expected vs Observed Preferred Frequencies\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  if (!is.data.table(csp_df)) setDT(csp_df)
  if (!is.data.table(polymorphism_df)) setDT(polymorphism_df)
  if (!is.data.table(expression_df)) setDT(expression_df)
  
  # Define expression bins
  expression_df[, Exp_Bin := cut(Exp_log10,
                                  breaks = quantile(Exp_log10, 
                                                    probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1),
                                                    na.rm = TRUE),
                                  labels = c("0-10%", "10-25%", "25-50%", "50-75%", "75-90%", "90-100%"),
                                  include.lowest = TRUE)]
  
  results_list <- list()
  
  # For each expression bin
  for (bin in unique(expression_df$Exp_Bin)) {
    
    genes_in_bin <- expression_df[Exp_Bin == bin, Gene]
    mean_exp <- expression_df[Exp_Bin == bin, mean(Exp_log10, na.rm = TRUE)]
    phi_linear <- 10^mean_exp
    
    # ROC model expected frequency of preferred codon for each AA
    roc_expected <- csp_df[, {
      log_unnorm <- -dM - dEta * phi_linear
      max_log <- max(log_unnorm)
      log_Z <- max_log + log(sum(exp(log_unnorm - max_log)))
      probs <- exp(log_unnorm - log_Z)
      
      preferred_prob <- sum(probs[is_optimal])
      
      list(P_Preferred_ROC = preferred_prob)
    }, by = AA]
    
    # Observed from polymorphism data
    # Frequency of preferred allele at segregating sites
    poly_in_bin <- polymorphism_df[Gene %in% genes_in_bin]
    
    if (nrow(poly_in_bin) > 0) {
      obs_freq <- poly_in_bin[, .(
        Mean_P_Preferred_Obs = mean(P_Preferred, na.rm = TRUE),
        SE = sd(P_Preferred, na.rm = TRUE) / sqrt(.N),
        N_Sites = .N
      ), by = AA]
      
      # Merge
      comparison <- merge(roc_expected, obs_freq, by = "AA", all.x = TRUE)
      comparison[, Exp_Bin := bin]
      comparison[, Mean_Exp := mean_exp]
      
      results_list[[bin]] <- comparison
    }
  }
  
  if (length(results_list) == 0) {
    cat("No polymorphism data matched genes in expression bins.\n")
    return(NULL)
  }
  
  results <- rbindlist(results_list)
  
  # Summary
  cat("Comparison by Expression Bin:\n")
  summary_by_bin <- results[, .(
    Mean_ROC = mean(P_Preferred_ROC, na.rm = TRUE),
    Mean_Obs = mean(Mean_P_Preferred_Obs, na.rm = TRUE),
    N_AA = sum(!is.na(Mean_P_Preferred_Obs)),
    Total_Sites = sum(N_Sites, na.rm = TRUE)
  ), by = .(Exp_Bin, Mean_Exp)]
  
  print(summary_by_bin[order(Mean_Exp)])
  
  # Correlation
  valid <- results[!is.na(P_Preferred_ROC) & !is.na(Mean_P_Preferred_Obs)]
  if (nrow(valid) > 5) {
    cor_test <- cor.test(valid$P_Preferred_ROC, valid$Mean_P_Preferred_Obs,
                         method = "pearson")
    cat(sprintf("\nCorrelation (ROC expected vs Observed): r = %.3f (p = %.2e)\n",
                cor_test$estimate, cor_test$p.value))
  }
  
  return(results)
}

#' Generate Master Summary Table
#' 
#' @param results_list List of results from all analyses
#' @return Data.frame with summary
generate_summary_table <- function(results_list) {
  
  summary_table <- data.frame(
    Analysis = character(),
    Result = character(),
    Consistent = character(),
    Interpretation = character(),
    stringsAsFactors = FALSE
  )
  
  # Analysis 1
  if (!is.null(results_list$analysis1)) {
    rho <- results_list$analysis1$correlation$estimate
    consistent <- ifelse(rho > 0.3, "Yes", ifelse(rho > 0.1, "Partial", "No"))
    summary_table <- rbind(summary_table, data.frame(
      Analysis = "1. Gene-level S_ROC vs γ correlation",
      Result = sprintf("ρ = %.3f", rho),
      Consistent = consistent,
      Interpretation = ifelse(consistent == "Yes",
                              "Both methods identify same genes",
                              "Methods measure different things")
    ))
  }
  
  # Analysis 2
  if (!is.null(results_list$analysis2)) {
    gamma_diff <- with(results_list$analysis2$summary,
                       Mean_Gamma[Expression_Group == "Top 5%"] -
                         Mean_Gamma[Expression_Group == "Bottom 5%"])
    consistent <- ifelse(mean(gamma_diff, na.rm = TRUE) > 0.5, "Yes",
                         ifelse(mean(gamma_diff, na.rm = TRUE) > 0.2, "Partial", "No"))
    summary_table <- rbind(summary_table, data.frame(
      Analysis = "2. Expression-stratified γ",
      Result = sprintf("Δγ = %.2f", mean(gamma_diff, na.rm = TRUE)),
      Consistent = consistent,
      Interpretation = ifelse(consistent == "Yes",
                              "γ higher in high-expression genes",
                              "γ uniform across expression")
    ))
  }
  
  # Analysis 4
  if (!is.null(results_list$analysis4)) {
    summary_table <- rbind(summary_table, data.frame(
      Analysis = "4. S_ROC-stratified γ estimation",
      Result = sprintf("Median γ: High=%.2f, Mid=%.2f, Low=%.2f",
                       results_list$analysis4$summary$Gamma_Median[1],
                       results_list$analysis4$summary$Gamma_Median[2],
                       results_list$analysis4$summary$Gamma_Median[3]),
      Consistent = ifelse(
        results_list$analysis4$summary$Gamma_Median[1] >
          results_list$analysis4$summary$Gamma_Median[3] + 0.3,
        "Yes", "No"),
      Interpretation = "Most direct test of ROC-SFS consistency"
    ))
  }
  
  # Analysis 5 - Per-AA gamma
  if (!is.null(results_list$analysis5)) {
    summary_table <- rbind(summary_table, data.frame(
      Analysis = "5. Per-AA gamma estimation (KEY FIX)",
      Result = sprintf("Mean γ = %.2f (n=%d AAs)",
                       mean(results_list$analysis5$Gamma, na.rm = TRUE),
                       sum(!is.na(results_list$analysis5$Gamma))),
      Consistent = "KEY TEST",
      Interpretation = "Proper AA-specific α,β corrects pooling error"
    ))
  }
  
  # ROC GOF
  if (!is.null(results_list$roc_gof)) {
    summary_table <- rbind(summary_table, data.frame(
      Analysis = "ROC Model Goodness-of-Fit",
      Result = sprintf("Median RMSE = %.4f, %d%% significant",
                       median(results_list$roc_gof$summary$RMSE),
                       round(100 * mean(results_list$roc_gof$summary$P_Value < 0.05))),
      Consistent = ifelse(
        mean(results_list$roc_gof$summary$P_Value < 0.05) < 0.3,
        "Yes", "Partial"),
      Interpretation = "Tests if ROC model fits observed codon frequencies"
    ))
  }
  
  return(summary_table)
}
