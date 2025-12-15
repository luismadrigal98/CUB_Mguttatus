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
  #' @param preferred_codons_df Data frame with columns: AA, Codon (or Preferred_Codon)
  #' @return Same data frame with added column: Terminal_Nucleotide
  #' ___________________________________________________________________________
  
  require(data.table)
  if (!is.data.table(preferred_codons_df)) {
    preferred_codons_df <- as.data.table(preferred_codons_df)
  } else {
    preferred_codons_df <- copy(preferred_codons_df)
  }
  
  # Handle both column naming conventions
  codon_col <- if ("Preferred_Codon" %in% names(preferred_codons_df)) {
    "Preferred_Codon"
  } else if ("Codon" %in% names(preferred_codons_df)) {
    "Codon"
  } else {
    stop("preferred_codons_df must have either 'Codon' or 'Preferred_Codon' column")
  }
  
  # Handle AA column (might be 'AA' or 'aa')
  # Check if we have single-letter codes or full names
  aa_col <- if ("AA" %in% names(preferred_codons_df)) {
    "AA"
  } else if ("aa" %in% names(preferred_codons_df)) {
    "aa"
  } else {
    stop("preferred_codons_df must have either 'AA' or 'aa' column")
  }
  
  # Check if we need to convert AA names to single letters
  sample_aa <- preferred_codons_df[[aa_col]][1]
  needs_conversion <- nchar(sample_aa) > 1  # If multi-character, it's a full name
  
  if (needs_conversion) {
    cat("Converting amino acid names to single-letter codes...\n")
    
    # Create mapping from full names to single letters
    aa_name_to_letter <- c(
      "Ala" = "A", "Arg_2" = "R", "Arg_4" = "R", "Arg_6" = "R",
      "Asn" = "N", "Asp" = "D", "Cys" = "C", "Gln" = "Q", "Glu" = "E",
      "Gly" = "G", "His" = "H", "Ile" = "I", 
      "Leu_2" = "L", "Leu_4" = "L", "Leu_6" = "L",
      "Lys" = "K", "Met" = "M", "Phe" = "F", "Pro" = "P",
      "Ser_2" = "S", "Ser_4" = "S", "Ser_6" = "S",
      "Thr" = "T", "Trp" = "W", "Tyr" = "Y", "Val" = "V"
    )
    
    # Add the column if using 'aa'
    if (aa_col == "aa") {
      preferred_codons_df[, AA := aa_name_to_letter[get(aa_col)]]
    } else {
      preferred_codons_df[, AA := aa_name_to_letter[AA]]
    }
    
    # Remove unmapped values
    if (any(is.na(preferred_codons_df$AA))) {
      unmapped <- unique(preferred_codons_df[is.na(AA)][[aa_col]])
      cat(sprintf("⚠️  Warning: Could not map amino acids: %s\n", 
                  paste(unmapped, collapse = ", ")))
      preferred_codons_df <- preferred_codons_df[!is.na(AA)]
    }
  }
  
  # Extract terminal nucleotide
  preferred_codons_df[, Terminal_Nucleotide := substr(
    get(codon_col), 
    nchar(get(codon_col)), 
    nchar(get(codon_col))
  )]
  
  # Validation: All should be G or C
  if (!all(preferred_codons_df$Terminal_Nucleotide %in% c("G", "C"))) {
    # Show which codons are problematic
    bad_codons <- preferred_codons_df[!Terminal_Nucleotide %in% c("G", "C")]
    cat("\n⚠️  ERROR: Some preferred codons do not end in G or C:\n")
    print(bad_codons)
    stop("All preferred codons must end in G or C for this analysis")
  }
  
  # Standardize column names
  if (codon_col == "Codon") {
    setnames(preferred_codons_df, "Codon", "Preferred_Codon")
  }
  
  # Debug output
  cat("Preferred codons table after standardization:\n")
  print(head(preferred_codons_df[, .(AA, Preferred_Codon, Terminal_Nucleotide)]))
  cat(sprintf("\nUnique amino acids: %d\n", length(unique(preferred_codons_df$AA))))
  cat("Amino acids present (single-letter codes):\n")
  print(sort(unique(preferred_codons_df$AA)))
  cat("\n")
  
  return(preferred_codons_df)
}

estimate_gamma_by_gene_with_neutral_params <- function(codon_vcf_data, 
                                                       neutral_params,
                                                       preferred_codons_df) {
  #' Estimate gamma (selection coefficient) for each gene and amino acid
  #' using pre-computed alpha and beta from intronic data
  #' 
  #' @param codon_vcf_data Output from prepare_vcf_for_gamma_estimation()
  #'                       Must have columns: Gene, AA, Preferred_Codon, k, n, p
  #' @param neutral_params Output from load_and_estimate_neutral_params()
  #' @param preferred_codons_df Data frame mapping AA to preferred codons
  #'                            Columns: AA, Preferred_Codon (or Codon)
  #' @return Data table with gamma estimates per gene and amino acid
  #' ___________________________________________________________________________
  
  require(data.table)
  
  if (!is.data.table(codon_vcf_data)) setDT(codon_vcf_data)
  if (!is.data.table(preferred_codons_df)) setDT(preferred_codons_df)
  
  cat("\n=== Estimating Selection Coefficients (Gamma) ===\n\n")
  
  # Annotate preferred codons with terminal nucleotide
  pref_annotated <- annotate_preferred_codons_with_nucleotide(preferred_codons_df)
  
  # Keep only the columns we need for merging
  pref_minimal <- pref_annotated[, .(AA, Preferred_Codon, Terminal_Nucleotide)]
  
  # Add terminal nucleotide to VCF data by matching on BOTH AA and Preferred_Codon
  # This handles cases where AA might have multiple preferred codons (shouldn't happen but safe)
  setkey(codon_vcf_data, AA, Preferred_Codon)
  setkey(pref_minimal, AA, Preferred_Codon)
  
  data_merged <- merge(codon_vcf_data, pref_minimal, 
                      by = c("AA", "Preferred_Codon"), 
                      all.x = TRUE)
  
  # Check for missing matches
  n_missing <- sum(is.na(data_merged$Terminal_Nucleotide))
  if (n_missing > 0) {
    cat(sprintf("⚠️  WARNING: %d sites missing terminal nucleotide annotation\n", n_missing))
    cat("Amino acids without matches:\n")
    print(unique(data_merged[is.na(Terminal_Nucleotide), .(AA, Preferred_Codon)]))
    cat("\nRemoving these sites from analysis...\n\n")
    data_merged <- data_merged[!is.na(Terminal_Nucleotide)]
  }
  
  # Count amino acids by terminal nucleotide
  aa_by_nuc <- data_merged[, .(N_Sites = .N), by = .(AA, Terminal_Nucleotide)]
  
  cat(sprintf("Amino acids with G-ending preferred codons: %d\n", 
              nrow(aa_by_nuc[Terminal_Nucleotide == "G"])))
  cat(sprintf("Amino acids with C-ending preferred codons: %d\n\n", 
              nrow(aa_by_nuc[Terminal_Nucleotide == "C"])))
  
  cat("Terminal nucleotide distribution:\n")
  print(table(data_merged$Terminal_Nucleotide))
  cat("\n")
  
  # Estimate gamma for each Gene x AA combination
  cat("Running gamma estimation (this may take a few minutes)...\n")
  
  # Debug: Check first few rows
  cat("\nDebug - Sample of data going into gamma estimation:\n")
  print(head(data_merged[, .(Gene, AA, k, n, p, Terminal_Nucleotide)], 10))
  cat(sprintf("\nSample sizes (n) summary:\n"))
  print(summary(data_merged$n))
  cat(sprintf("\nPreferred codon counts (k) summary:\n"))
  print(summary(data_merged$k))
  
  # Check grouping
  cat(sprintf("\nNumber of rows in data_merged: %d\n", nrow(data_merged)))
  cat(sprintf("Number of unique Gene x AA combinations: %d\n", 
              nrow(unique(data_merged[, .(Gene, AA)]))))
  cat("\n")
  
  gamma_results <- data_merged[, {
    
    # Select appropriate alpha/beta based on terminal nucleotide
    term_nuc <- Terminal_Nucleotide[1]
    
    if (term_nuc == "G") {
      alpha_use <- neutral_params$alpha_G
      beta_use <- neutral_params$beta_G
    } else {  # C
      alpha_use <- neutral_params$alpha_C
      beta_use <- neutral_params$beta_C
    }
    
    # Always return same structure
    n_sites <- .N  # Number of polymorphic sites for this Gene x AA
    mean_p <- mean(p, na.rm = TRUE)
    
    # For single aggregated observation, we CAN estimate gamma
    # The k and n values represent the aggregated counts across all sites
    # We just need enough alleles sampled (n should be large)
    
    # Skip only if sample size is too small
    min_sample_size <- 100  # Need at least 100 alleles total
    
    if (n[1] < min_sample_size) {
      gamma_est <- NA_real_
      skip_reason <- "Sample size too small"
    } else {
      # Estimate gamma
      # NOTE: k and n are VECTORS if we have multiple observations
      # But in your case, each Gene x AA has ONE row with aggregated counts
      gamma_est <- tryCatch({
        estimate_gamma_for_AA(
          counts = k,
          sample_sizes = n,
          alpha = alpha_use,
          beta = beta_use,
          S_interval = c(-10, 30)
        )
      },
      error = function(e) {
        # Capture error for debugging
        if (!exists("first_error", envir = .GlobalEnv)) {
          assign("first_error", list(
            Gene = Gene[1],
            AA = AA[1],
            error = e$message,
            k = k,
            n = n,
            alpha = alpha_use,
            beta = beta_use
          ), envir = .GlobalEnv)
        }
        return(NA_real_)
      })
      skip_reason <- if (is.na(gamma_est)) "Optimization failed" else "Success"
    }
    
    # Return consistent structure
    list(
      Gamma = gamma_est,
      N_Sites = n_sites,
      Total_Alleles = sum(n),
      Terminal_Nuc = term_nuc,
      Mean_Freq_Preferred = mean_p,
      Alpha_Used = alpha_use,
      Beta_Used = beta_use
    )
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