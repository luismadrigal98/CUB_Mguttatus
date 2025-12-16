#' Integration of Intronic Polymorphism Data for Selection Inference
#' 
#' This module handles the workflow for estimating selection coefficients (gamma)
#' using mutation rate parameters (alpha, beta) derived from neutral intronic sites.
#' 
#' =============================================================================
#' COMPARISON WITH ANACODA SELECTION COEFFICIENTS
#' =============================================================================
#' 
#' Gamma (this approach):
#'   - Definition: gamma = 4N*s, selection coefficient favoring PREFERRED codon
#'   - Sign: POSITIVE = selection FOR preferred codon
#'           NEGATIVE = selection AGAINST preferred codon
#'   - Unit: Per amino acid position
#'   - Model: Wright-Fisher diffusion with biallelic mutation-selection balance
#'   - Interpretation: gamma > 1.92 means selection dominates drift (|4Nes| > 1)
#' 
#' AnaCoDa deltaEta (selection parameter):
#'   - Definition: Selection coefficient RELATIVE to reference codon
#'   - Sign: NEGATIVE for non-preferred codons (penalties)
#'           ZERO for reference codon
#'   - Unit: Per codon (not per AA position)
#'   - Model: ROC multinomial with mutation (deltaM) + selection (deltaEta)
#'   - Interpretation: More negative = stronger penalty = less preferred
#' 
#' Gene-Level Aggregation (for comparability):
#'   - AnaCoDa: Sum(|deltaEta_i| * Count_i) / Total_Synonymous_Codons
#'     --> Mean absolute selection per codon position
#'   
#'   - Gamma: Sum(|gamma_AA| * N_AA_positions) / Total_Synonymous_Codons
#'     --> Mean absolute selection per AA position
#'     OR: mean(|gamma_AA|) weighted by sample size
#'     
#' Expected Correlation:
#'   - Positive correlation: Both identify genes under selection for codon bias
#'   - Higher values in both = stronger selection
#'   - But magnitudes not directly comparable (different scales/models)
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

standardize_aa_names <- function(dt, aa_col = "AA") {
  #' Robustly convert amino acid names to single-letter codes
  #' Handles 3-letter codes (Ala, Ser) and split variants (Ser_2, Arg_4, Leu_6)
  #' 
  #' @param dt data.table to modify
  #' @param aa_col Name of the amino acid column (default: "AA")
  #' @return Modified data.table with standardized single-letter AA codes
  #' ___________________________________________________________________________
  
  require(data.table)
  
  # Comprehensive mapping including all split variants
  aa_map <- c(
    "Ala" = "A", "Arg" = "R", "Arg_2" = "R", "Arg_4" = "R", "Arg_6" = "R",
    "Asn" = "N", "Asp" = "D", "Cys" = "C", "Gln" = "Q", "Glu" = "E",
    "Gly" = "G", "His" = "H", "Ile" = "I", 
    "Leu" = "L", "Leu_2" = "L", "Leu_4" = "L", "Leu_6" = "L",
    "Lys" = "K", "Met" = "M", "Phe" = "F", "Pro" = "P",
    "Ser" = "S", "Ser_2" = "S", "Ser_4" = "S", "Ser_6" = "S",
    "Thr" = "T", "Trp" = "W", "Tyr" = "Y", "Val" = "V",
    # Already standardized (pass through)
    "A" = "A", "R" = "R", "N" = "N", "D" = "D", "C" = "C",
    "Q" = "Q", "E" = "E", "G" = "G", "H" = "H", "I" = "I",
    "L" = "L", "K" = "K", "M" = "M", "F" = "F", "P" = "P",
    "S" = "S", "T" = "T", "W" = "W", "Y" = "Y", "V" = "V"
  )
  
  # Check if column exists
  if (!aa_col %in% names(dt)) {
    stop(sprintf("Column '%s' not found in data.table", aa_col))
  }
  
  # Get original values for reporting
  original_values <- unique(dt[[aa_col]])
  
  # Update column in place
  dt[, (aa_col) := ifelse(get(aa_col) %in% names(aa_map), 
                          aa_map[get(aa_col)], 
                          get(aa_col))]
  
  # Report any unmapped values
  new_values <- unique(dt[[aa_col]])
  unmapped <- setdiff(original_values, names(aa_map))
  
  if (length(unmapped) > 0) {
    warning(sprintf("Unmapped amino acid codes found: %s", 
                    paste(unmapped, collapse = ", ")))
  }
  
  # Report conversion summary
  if (!all(original_values %in% c("A","R","N","D","C","Q","E","G","H","I",
                                   "L","K","M","F","P","S","T","W","Y","V"))) {
    cat(sprintf("✓ Standardized AA names: %s\n", 
                paste(original_values[1:min(3, length(original_values))], collapse = ", ")))
  }
  
  return(dt)
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
  
  # FIX 2: Standardize AA names before processing
  cat("Standardizing amino acid names...\n")
  codon_vcf_data <- standardize_aa_names(codon_vcf_data, "AA")
  preferred_codons_df <- standardize_aa_names(preferred_codons_df, "AA")
  
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
    
    # Get site-level data
    n_sites <- .N  # Number of polymorphic sites for this Gene x AA
    mean_p <- mean(p, na.rm = TRUE)
    
    # CRITICAL FIX: We now have MULTIPLE sites per Gene×AA
    # k and n are VECTORS with one element per codon position
    # The likelihood function will integrate over all sites
    
    # Quality check: Need sufficient sites AND reasonable sample sizes
    # FIX: Lowered from 5 to 1 to unlock data for genes with few sites
    min_sites <- 1        # Changed from 5 to 1 (85% of genes have <5 sites per AA)
    min_sample_per_site <- 50  # Each site needs ≥50 alleles
    
    # Check if we have enough data
    if (n_sites < min_sites) {
      gamma_est <- NA_real_
      skip_reason <- "Too few polymorphic sites"
    } else if (any(n < min_sample_per_site)) {
      gamma_est <- NA_real_
      skip_reason <- "Some sites have insufficient sample size"
    } else {
      # Estimate gamma using VECTORS of k and n across all sites
      # This properly accounts for variance between sites (drift)
      gamma_est <- tryCatch({
        estimate_gamma_for_AA(
          counts = k,              # Vector of preferred counts
          sample_sizes = n,        # Vector of sample sizes
          alpha = alpha_use,
          beta = beta_use,
          S_interval = c(0, 50)    # Positive only: comparable to AnaCoDa |4Nes|
        )
      },
      error = function(e) {
        # Capture error for debugging
        if (!exists("first_error", envir = .GlobalEnv)) {
          assign("first_error", list(
            Gene = Gene[1],
            AA = AA[1],
            error = e$message,
            n_sites = n_sites,
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
      N_Sites = n_sites,                    # Now correctly shows multiple sites
      Total_Alleles = sum(n),               # Sum across all sites
      Mean_Alleles_Per_Site = mean(n),      # NEW: Average sample size per site
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
  
  # CRITICAL VALIDATION: Check data structure
  cat("=== DATA STRUCTURE VALIDATION ===\n")
  cat(sprintf("Mean sites per Gene×AA: %.1f\n", mean(gamma_results$N_Sites, na.rm = TRUE)))
  cat(sprintf("Median sites per Gene×AA: %.0f\n", median(gamma_results$N_Sites, na.rm = TRUE)))
  cat(sprintf("Gene×AA with <5 sites: %d (%.1f%% - should be low!)\n",
              sum(gamma_results$N_Sites < 5, na.rm = TRUE),
              100 * mean(gamma_results$N_Sites < 5, na.rm = TRUE)))
  
  cat(sprintf("\nMean alleles per site: %.1f \n", 
              mean(gamma_results$Mean_Alleles_Per_Site, na.rm = TRUE)))
  
  if (mean(gamma_results$N_Sites, na.rm = TRUE) < 2) {
    warning("\n⚠️  CRITICAL ERROR: N_Sites is too low! Data is still collapsed!\n")
  }
  if (mean(gamma_results$Mean_Alleles_Per_Site, na.rm = TRUE) > 1000) {
    warning("\n⚠️  CRITICAL ERROR: Sample sizes are too large! Sites are being pooled!\n")
  }
  cat("\n")
  
  # Summary by terminal nucleotide
  summary_by_nuc <- gamma_results[!is.na(Gamma), .(
    N = .N,
    Mean_Gamma = mean(Gamma),
    Median_Gamma = median(Gamma),
    SD_Gamma = sd(Gamma),
    Prop_Positive = mean(Gamma > 0),
    Prop_Significant = mean(Significant),
    Mean_Sites = mean(N_Sites),
    Mean_Sample_Per_Site = mean(Mean_Alleles_Per_Site)
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

aggregate_gamma_per_gene <- function(gamma_results, codon_usage_df, genetic_code) {
  #' Aggregate gamma values to gene-level selection intensity
  #' Comparable to AnaCoDa's selection coefficient aggregation
  #' 
  #' Logic for comparability:
  #' - AnaCoDa: Sum(|deltaEta_i| * Count_i) / Total_Synonymous_Codons
  #'   where deltaEta_i is NEGATIVE for non-preferred codons
  #'   
  #' - Gamma: Sum(|gamma_AA| * N_AA_positions) / Total_Synonymous_Codons
  #'   where gamma is POSITIVE for selection favoring preferred codon
  #'   
  #' Both measure "mean selection intensity per codon position"
  #' Higher values = stronger selection for codon bias
  #' 
  #' @param gamma_results Output from estimate_gamma_by_gene_with_neutral_params()
  #'                      Must have columns: Gene, AA, Gamma
  #' @param codon_usage_df Codon usage matrix with amino acid counts per gene
  #'                       Must have: Gene_name column + codon columns
  #' @param genetic_code Named vector or data.table mapping codons to AAs
  #' @return Data table with gene-level selection metrics
  #' ___________________________________________________________________________
  
  require(data.table)
  
  if (!is.data.table(gamma_results)) gamma_results <- as.data.table(gamma_results)
  if (!is.data.table(codon_usage_df)) codon_usage_df <- as.data.table(codon_usage_df)
  
  # Standardize gene names
  gamma_results <- gamma_results |>
    dplyr::mutate(Gene_name = paste0("MgIM767.", Gene))
  
  cat("\n=== Aggregating Gamma to Gene-Level Selection Intensity ===\n\n")
  
  # Convert genetic code to mapping
  if (is.vector(genetic_code)) {
    codon_to_aa <- genetic_code
  } else if (is.data.table(genetic_code) || is.data.frame(genetic_code)) {
    codon_to_aa <- setNames(genetic_code$AA, genetic_code$Codon)
  } else {
    stop("genetic_code must be a named vector or data.table with Codon and AA columns")
  }
  
  # Get all codon columns from usage matrix
  codon_cols <- setdiff(names(codon_usage_df), "Gene_name")
  
  cat(sprintf("Converting codon counts to AA counts for %d genes...\n", 
              nrow(codon_usage_df)))
  
  # Convert codon usage to AA usage
  aa_usage <- codon_usage_df[, {
    
    gene_name <- Gene_name
    
    # Sum codons by amino acid
    aa_counts <- list()
    
    for (codon in codon_cols) {
      aa <- codon_to_aa[codon]
      
      # Skip if AA not found or is STOP
      if (is.na(aa) || aa == "STOP") next
      
      count <- .SD[[codon]]
      
      if (aa %in% names(aa_counts)) {
        aa_counts[[aa]] <- aa_counts[[aa]] + count
      } else {
        aa_counts[[aa]] <- count
      }
    }
    
    # Convert to data.table format
    data.table(
      AA = names(aa_counts),
      AA_Count = unlist(aa_counts)
    )
    
  }, by = Gene_name, .SDcols = codon_cols]
  
  # Calculate total synonymous codons per gene (exclude Met and Trp)
  gene_lengths <- aa_usage[!AA %in% c("M", "W"), .(
    Total_Synonymous_Codons = sum(AA_Count)
  ), by = Gene_name]
  
  # Merge gamma results with AA counts
  setkey(gamma_results, Gene_name, AA)
  setkey(aa_usage, Gene_name, AA)
  
  gamma_with_counts <- merge(gamma_results[!is.na(Gamma)], aa_usage, 
                             by = c("Gene_name", "AA"), 
                             all.x = TRUE)
  
  # Handle missing AA counts (shouldn't happen but safe)
  if (any(is.na(gamma_with_counts$AA_Count))) {
    cat("⚠️  Warning: Some AAs in gamma_results not found in codon_usage\n")
    gamma_with_counts <- gamma_with_counts[!is.na(AA_Count)]
  }
  
  # Calculate weighted aggregates per gene
  gamma_weighted <- gamma_with_counts[, .(
    
    # CORRECT FORMULA: Weight by actual AA counts from gene sequence
    # Selection_Load = Sum(|gamma_AA| * AA_Count_AA)
    Total_Selection_Load = sum(abs(Gamma) * AA_Count, na.rm = TRUE),
    
    # Weighted mean gamma (weighted by AA occurrence)
    Gamma_Weighted_Mean = weighted.mean(Gamma, w = AA_Count, na.rm = TRUE),
    
    # Simple mean gamma (all AAs equally weighted - for comparison)
    Gamma_Mean = mean(Gamma, na.rm = TRUE),
    
    # Median gamma (robust to outliers)
    Gamma_Median = median(Gamma, na.rm = TRUE),
    
    # Count statistics
    N_AA_Analyzed = .N,
    N_AA_Positions = sum(AA_Count),  # Total AA positions with gamma estimates
    N_Positive_Gamma = sum(Gamma > 0, na.rm = TRUE),
    N_Negative_Gamma = sum(Gamma < 0, na.rm = TRUE),
    N_Significant = sum(Significant, na.rm = TRUE),
    
    # Proportion metrics
    Prop_Positive = mean(Gamma > 0, na.rm = TRUE),
    Prop_Significant = mean(Significant, na.rm = TRUE),
    
    # Range
    Max_Gamma = max(Gamma, na.rm = TRUE),
    Min_Gamma = min(Gamma, na.rm = TRUE)
    
  ), by = Gene_name]
  
  # Merge with gene lengths to calculate Selection_Intensity
  setkey(gamma_weighted, Gene_name)
  setkey(gene_lengths, Gene_name)
  
  gamma_weighted <- gene_lengths[gamma_weighted, nomatch = 0]
  
  # Calculate Selection Intensity (normalized by total synonymous codons)
  # This is DIRECTLY comparable to AnaCoDa's S_coeff
  gamma_weighted[, Selection_Intensity := Total_Selection_Load / Total_Synonymous_Codons]
  
  # Summary
  cat(sprintf("Genes with gamma estimates: %d\n", nrow(gamma_weighted)))
  cat(sprintf("Median AAs per gene: %.0f\n", 
              median(gamma_weighted$N_AA_Analyzed)))
  cat(sprintf("Mean AA positions per gene: %.0f\n",
              mean(gamma_weighted$N_AA_Positions)))
  cat("\nSelection Intensity Distribution:\n")
  print(summary(gamma_weighted$Selection_Intensity))
  cat("\nWeighted Mean Gamma Distribution:\n")
  print(summary(gamma_weighted$Gamma_Weighted_Mean))
  cat("\n")
  
  # Compare distributions
  cat("Positive vs Negative Selection:\n")
  cat(sprintf("  Genes with net positive gamma: %d (%.1f%%)\n",
              sum(gamma_weighted$Gamma_Mean > 0),
              100 * mean(gamma_weighted$Gamma_Mean > 0)))
  cat(sprintf("  Genes with net negative gamma: %d (%.1f%%)\n\n",
              sum(gamma_weighted$Gamma_Mean < 0),
              100 * mean(gamma_weighted$Gamma_Mean < 0)))
  
  return(gamma_weighted)
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

compare_gamma_with_anacoda <- function(gamma_gene_level, anacoda_intensity) {
  #' Direct comparison of gamma-based vs AnaCoDa-based selection metrics
  #' 
  #' **DEPRECATED:** This function uses incorrect aggregation.
  #' Use contrast_gamma_anacoda() instead for mathematically correct comparison.
  #' 
  #' This function is kept for backward compatibility but will issue a warning
  #' and suggest using the correct function.
  #' 
  #' @param gamma_gene_level Output from aggregate_gamma_per_gene()
  #'                         Must have: Gene_name, Selection_Intensity
  #' @param anacoda_intensity AnaCoDa selection intensity data
  #'                          Must have: Gene_name, S_coeff
  #' @return Merged comparison table
  #' ___________________________________________________________________________
  
  require(data.table)
  
  if (!is.data.table(gamma_gene_level)) setDT(gamma_gene_level)
  if (!is.data.table(anacoda_intensity)) setDT(anacoda_intensity)
  
  cat("\n⚠️  WARNING: compare_gamma_with_anacoda() uses INCORRECT aggregation!\n")
  cat("   This function compares gene-level aggregates, not codon-level selection.\n")
  cat("   For mathematically correct comparison, use:\n")
  cat("     contrast_gamma_anacoda(gamma_results, codon_usage, preferred_codons,\n")
  cat("                            anacoda_intensity, genetic_code)\n\n")
  cat("   Continuing with simple aggregate comparison for diagnostic purposes...\n\n")
  
  cat("\n=== Comparing Gamma vs AnaCoDa Selection Coefficients ===\n")
  cat("   (Note: This is a simplified comparison, not the rigorous formula)\n\n")
  
  # Merge datasets
  setkey(gamma_gene_level, Gene_name)
  setkey(anacoda_intensity, Gene_name)
  
  comparison <- anacoda_intensity[gamma_gene_level, nomatch = 0]
  
  cat(sprintf("Genes in both analyses: %d\n\n", nrow(comparison)))
  
  # Key comparison: Selection_Intensity (gamma) vs S_coeff (AnaCoDa)
  # Both should be POSITIVE and measure "mean selection per codon"
  
  if ("S_coeff" %in% names(comparison) && "Selection_Intensity" %in% names(comparison)) {
    
    cor_intensity <- cor.test(comparison$Selection_Intensity, 
                              comparison$S_coeff,
                              method = "spearman", exact = FALSE)
    
    cat("Correlation: Gamma Selection_Intensity vs AnaCoDa S_coeff\n")
    cat(sprintf("  Spearman ρ = %.3f\n", cor_intensity$estimate))
    cat(sprintf("  p-value = %.2e\n\n", cor_intensity$p.value))
    
    # Interpretation guide
    if (cor_intensity$estimate > 0.5 && cor_intensity$p.value < 0.01) {
      cat("✓ Strong positive correlation - methods are consistent!\n")
      cat("  Both identify the same genes under strong codon bias selection.\n\n")
    } else if (cor_intensity$estimate > 0.3 && cor_intensity$p.value < 0.05) {
      cat("⚠ Moderate positive correlation - generally consistent.\n")
      cat("  Methods agree on general trends but may differ in details.\n\n")
    } else {
      cat("✗ Weak or no correlation - methods may be measuring different aspects.\n")
      cat("  Check for systematic differences in gene sets or AA coverage.\n\n")
    }
    
    # Quartile comparison
    q_gamma <- quantile(comparison$Selection_Intensity, c(0.25, 0.75))
    q_anacoda <- quantile(comparison$S_coeff, c(0.25, 0.75))
    
    cat("Distribution Comparison:\n")
    cat(sprintf("  Gamma Intensity:  Q1=%.3f, Median=%.3f, Q3=%.3f\n",
                q_gamma[1], 
                median(comparison$Selection_Intensity),
                q_gamma[2]))
    cat(sprintf("  AnaCoDa S_coeff:  Q1=%.3f, Median=%.3f, Q3=%.3f\n\n",
                q_anacoda[1],
                median(comparison$S_coeff),
                q_anacoda[2]))
  }
  
  # Also compare weighted means
  if ("Gamma_Weighted_Mean" %in% names(comparison) && "S_coeff" %in% names(comparison)) {
    
    cor_weighted <- cor.test(comparison$Gamma_Weighted_Mean,
                             comparison$S_coeff,
                             method = "spearman", exact = FALSE)
    
    cat("Correlation: Gamma Weighted_Mean vs AnaCoDa S_coeff\n")
    cat(sprintf("  Spearman ρ = %.3f\n", cor_weighted$estimate))
    cat(sprintf("  p-value = %.2e\n\n", cor_weighted$p.value))
  }
  
  return(comparison)
}

contrast_gamma_anacoda <- function(gamma_results, codon_usage, preferred_codons, 
                                   anacoda_intensity, genetic_code) {
  #' Mathematical contrast of polymorphism-based gamma vs AnaCoDa selection
  #' 
  #' Implements the rigorous aggregation formula:
  #'   S_poly = (1/L) * Sum_AA(Count_Unpref_AA * gamma_AA)
  #' 
  #' where:
  #'   - L = total gene length in codons
  #'   - Count_Unpref_AA = number of UNPREFERRED codons for each AA in the gene
  #'   - gamma_AA = selection coefficient favoring the preferred codon
  #' 
  #' This measures the "selection load" from using unpreferred codons,
  #' directly comparable to AnaCoDa's selection intensity (phi * deltaEta)
  #' 
  #' @param gamma_results Output from estimate_gamma_by_gene_with_neutral_params()
  #'                      Must have: Gene, AA, Gamma
  #' @param codon_usage Codon usage matrix (genes × codons)
  #'                    Must have: Gene_name, [codon columns]
  #' @param preferred_codons Table mapping AA to preferred codon
  #'                         Must have: AA, Preferred_Codon (or Codon)
  #' @param anacoda_intensity AnaCoDa selection intensity
  #'                          Must have: Gene_name, S_coeff
  #' @param genetic_code Mapping of codons to amino acids
  #'                     Must have: Codon, AA (or named vector)
  #' @return Merged data.table with S_poly, S_coeff, correlation stats, and plot
  #' ___________________________________________________________________________
  
  require(data.table)
  require(ggplot2)
  
  if (!is.data.table(gamma_results)) gamma_results <- as.data.table(gamma_results)
  if (!is.data.table(codon_usage)) codon_usage <- as.data.table(codon_usage)
  if (!is.data.table(preferred_codons)) preferred_codons <- as.data.table(preferred_codons)
  if (!is.data.table(anacoda_intensity)) anacoda_intensity <- as.data.table(anacoda_intensity)
  
  cat("\n=== Mathematical Contrast: Polymorphism vs AnaCoDa ===\n\n")
  
  # Convert genetic code to data.table if needed
  if (is.vector(genetic_code)) {
    genetic_code <- data.table(
      Codon = names(genetic_code),
      AA = as.character(genetic_code)
    )
  }
  if (!is.data.table(genetic_code)) genetic_code <- as.data.table(genetic_code)
  
  # Get preferred codon for each AA
  pref_col <- if ("Preferred_Codon" %in% names(preferred_codons)) {
    "Preferred_Codon"
  } else {
    "Codon"
  }
  
  aa_to_pref <- setNames(preferred_codons[[pref_col]], preferred_codons$AA)
  
  # Create codon-to-AA mapping
  codon_to_aa <- setNames(genetic_code$AA, genetic_code$Codon)
  
  # Get all codon columns from usage matrix
  codon_cols <- setdiff(names(codon_usage), "Gene_name")
  
  cat(sprintf("Processing %d genes with %d codons...\n", 
              nrow(codon_usage), length(codon_cols)))
  
  # Standardize gene names in gamma_results
  # gamma_results has 'Gene', codon_usage has 'Gene_name'
  gamma_for_lookup <- copy(gamma_results)
  if ("Gene" %in% names(gamma_for_lookup) && !"Gene_name" %in% names(gamma_for_lookup)) {
    setnames(gamma_for_lookup, "Gene", "Gene_name")
  }
  
  cat(sprintf("Gamma results: %d rows, %d unique genes\n",
              nrow(gamma_for_lookup),
              length(unique(gamma_for_lookup$Gene_name))))
  
  # Calculate S_poly for each gene
  cat("Calculating polymorphism-based selection load (S_poly)...\n")
  
  S_poly_results <- codon_usage[, {
    
    gene_name <- Gene_name
    
    # Calculate total gene length
    total_length <- sum(unlist(.SD), na.rm = TRUE)
    
    if (total_length == 0) {
      list(S_poly = NA_real_, Gene_Length = 0)
    } else {
      
      # Calculate selection load from unpreferred codons
      selection_load <- 0
      
      for (codon in codon_cols) {
        # Get AA for this codon
        aa <- codon_to_aa[codon]
        
        # Skip if AA not found or is Met/Trp/STOP
        if (is.na(aa) || aa %in% c("M", "W", "STOP")) next
        
        # Check if this codon is the PREFERRED one for its AA
        pref_codon <- aa_to_pref[aa]
        
        # If this is an UNPREFERRED codon
        if (!is.na(pref_codon) && codon != pref_codon) {
          
          # Get gamma for this AA in this gene (now using Gene_name)
          gamma_val <- gamma_for_lookup[Gene_name == gene_name & AA == aa, Gamma]
          
          if (length(gamma_val) > 0 && !is.na(gamma_val[1])) {
            # Count of this unpreferred codon
            count_unpref <- .SD[[codon]]
            
            # Add to selection load: count * gamma
            selection_load <- selection_load + (count_unpref * gamma_val[1])
          }
        }
      }
      
      # Normalize by gene length
      s_poly <- selection_load / total_length
      
      list(S_poly = s_poly, Gene_Length = total_length)
    }
    
  }, by = Gene_name, .SDcols = codon_cols]
  
  # Merge with AnaCoDa results
  setkey(S_poly_results, Gene_name)
  setkey(anacoda_intensity, Gene_name)
  
  contrast_data <- anacoda_intensity[S_poly_results, nomatch = 0]
  
  # Remove genes with missing data
  contrast_data <- contrast_data[!is.na(S_poly) & !is.na(S_coeff)]
  
  cat(sprintf("\nGenes in comparison: %d\n", nrow(contrast_data)))
  
  # Diagnostic: Check variance
  cat(sprintf("S_poly: Mean=%.4f, SD=%.4f, Range=[%.4f, %.4f]\n",
              mean(contrast_data$S_poly), sd(contrast_data$S_poly),
              min(contrast_data$S_poly), max(contrast_data$S_poly)))
  cat(sprintf("S_coeff: Mean=%.4f, SD=%.4f, Range=[%.4f, %.4f]\n\n",
              mean(contrast_data$S_coeff), sd(contrast_data$S_coeff),
              min(contrast_data$S_coeff), max(contrast_data$S_coeff)))
  
  # Calculate Spearman correlation
  if (nrow(contrast_data) > 10) {
    
    # Check if there's variance in both variables
    if (sd(contrast_data$S_poly) == 0) {
      warning("⚠️  S_poly has zero variance! All values are identical.")
      cat("This means no selection load was calculated.\n")
      cat("Possible issues:\n")
      cat("  1. No matching genes between gamma_results and codon_usage\n")
      cat("  2. All gamma values are NA\n")
      cat("  3. Column name mismatch\n\n")
      return(contrast_data)
    }
    
    if (sd(contrast_data$S_coeff) == 0) {
      warning("⚠️  S_coeff has zero variance! All AnaCoDa values are identical.")
      return(contrast_data)
    }
    
    cor_result <- cor.test(contrast_data$S_poly, 
                           contrast_data$S_coeff,
                           method = "spearman", 
                           exact = FALSE)
    
    cat("=== Correlation Analysis ===\n")
    cat(sprintf("Spearman ρ = %.4f\n", cor_result$estimate))
    cat(sprintf("p-value = %.2e\n\n", cor_result$p.value))
    
    # Interpretation
    if (!is.na(cor_result$estimate) && cor_result$estimate > 0.5 && cor_result$p.value < 0.01) {
      cat("✓ STRONG CONCORDANCE: Both methods identify the same selection signal!\n")
      cat("  Polymorphism-based inference validates AnaCoDa mechanistic model.\n\n")
    } else if (!is.na(cor_result$estimate) && cor_result$estimate > 0.3 && cor_result$p.value < 0.05) {
      cat("⚠ MODERATE CONCORDANCE: General agreement with some discrepancies.\n")
      cat("  Methods capture similar but not identical aspects of selection.\n\n")
    } else {
      cat("✗ WEAK CONCORDANCE: Methods may be measuring different signals.\n")
      cat("  Investigate systematic differences between approaches.\n\n")
    }
    
    # Summary statistics
    cat("=== Distribution Summary ===\n")
    cat(sprintf("S_poly (Polymorphism):\n"))
    cat(sprintf("  Mean = %.4f, Median = %.4f\n", 
                mean(contrast_data$S_poly), median(contrast_data$S_poly)))
    cat(sprintf("  Q1 = %.4f, Q3 = %.4f\n\n", 
                quantile(contrast_data$S_poly, 0.25), 
                quantile(contrast_data$S_poly, 0.75)))
    
    cat(sprintf("S_coeff (AnaCoDa):\n"))
    cat(sprintf("  Mean = %.4f, Median = %.4f\n", 
                mean(contrast_data$S_coeff), median(contrast_data$S_coeff)))
    cat(sprintf("  Q1 = %.4f, Q3 = %.4f\n\n", 
                quantile(contrast_data$S_coeff, 0.25), 
                quantile(contrast_data$S_coeff, 0.75)))
    
    # Create scatter plot
    cat("Generating scatter plot...\n")
    
    p <- ggplot(contrast_data, aes(x = S_coeff, y = S_poly)) +
      geom_point(alpha = 0.4, size = 1.5, color = "steelblue") +
      geom_smooth(method = "lm", color = "red", se = TRUE, linewidth = 1.2) +
      labs(
        title = "Validation: Polymorphism-Based vs AnaCoDa Selection Estimates",
        subtitle = sprintf("Spearman ρ = %.3f, p = %.2e | n = %d genes",
                           cor_result$estimate, 
                           cor_result$p.value,
                           nrow(contrast_data)),
        x = expression(bar(S)[AnaCoDa]~"(Selection Intensity)"),
        y = expression(bar(S)[poly]~"(Polymorphism-Based Load)")
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11),
        axis.title = element_text(size = 12),
        panel.grid.minor = element_blank()
      )
    
    # Save plot
    ggsave("./results/gamma_anacoda_contrast.pdf", 
           plot = p, 
           width = 8, height = 7)
    
    cat("✓ Plot saved: ./results/gamma_anacoda_contrast.pdf\n\n")
    
    # Add correlation stats to output
    contrast_data[, Spearman_rho := cor_result$estimate]
    contrast_data[, Spearman_p := cor_result$p.value]
    
  } else {
    warning("Not enough data for correlation analysis (n < 10)")
  }
  
  return(contrast_data)
}

estimate_gamma_gradient <- function(codon_vcf_data, neutral_params, 
                                    preferred_codons_df, n_bins = 10) {
  #' Test for gBGC signature: Is selection stronger at 5' end of genes?
  #' 
  #' Gene conversion and recombination-driven gBGC are expected to be 
  #' stronger near the 5' end due to recombination rate gradients.
  #' This function bins sites by relative position within genes and
  #' estimates a global gamma for each positional bin.
  #' 
  #' @param codon_vcf_data Output from prepare_vcf_for_gamma_estimation()
  #' @param neutral_params Output from load_and_estimate_neutral_params()
  #' @param preferred_codons_df Data frame with AA, Preferred_Codon, Terminal_Nucleotide
  #' @param n_bins Number of positional bins (default: 10 = deciles)
  #' @return Data table with gamma estimates per positional bin
  #' ___________________________________________________________________________
  
  require(data.table)
  require(ggplot2)
  
  if (!is.data.table(codon_vcf_data)) setDT(codon_vcf_data)
  if (!is.data.table(preferred_codons_df)) setDT(preferred_codons_df)
  
  cat("\n=== Testing for gBGC Gradient Along Genes ===\n\n")
  
  # 1. Calculate relative position for each site (0.0 to 1.0)
  cat("Calculating relative positions...\n")
  
  # Calculate gene lengths first
  gene_lengths <- codon_vcf_data[, .(Gene_Length = max(Codon_Pos)), by = Gene]
  
  # Merge and calculate relative positions
  vcf_with_pos <- merge(codon_vcf_data, gene_lengths, by = "Gene")
  vcf_with_pos[, Relative_Position := Codon_Pos / Gene_Length]
  
  # 2. Assign each site to a positional bin
  vcf_with_pos[, Position_Bin := cut(Relative_Position, 
                                      breaks = seq(0, 1, length.out = n_bins + 1),
                                      labels = 1:n_bins,
                                      include.lowest = TRUE)]
  
  cat(sprintf("Sites binned into %d positional categories\n", n_bins))
  cat(sprintf("Sites per bin: %.0f (median)\n\n", 
              median(table(vcf_with_pos$Position_Bin))))
  
  # 3. Add terminal nucleotide annotation
  cat("Adding terminal nucleotide annotations...\n")
  
  # Standardize AA names first
  vcf_with_pos <- standardize_aa_names(vcf_with_pos, "AA")
  preferred_codons_df <- standardize_aa_names(preferred_codons_df, "AA")
  
  # Annotate preferred codons if not already done
  if (!"Terminal_Nucleotide" %in% names(preferred_codons_df)) {
    preferred_codons_df <- annotate_preferred_codons_with_nucleotide(preferred_codons_df)
  }
  
  # Get minimal columns for merging
  pref_minimal <- preferred_codons_df[, .(AA, Preferred_Codon, Terminal_Nucleotide)]
  
  # Merge by both AA and Preferred_Codon
  setkey(vcf_with_pos, AA, Preferred_Codon)
  setkey(pref_minimal, AA, Preferred_Codon)
  
  vcf_annotated <- merge(vcf_with_pos, pref_minimal, 
                         by = c("AA", "Preferred_Codon"),
                         all.x = TRUE)
  
  # Check for missing annotations
  n_missing <- sum(is.na(vcf_annotated$Terminal_Nucleotide))
  if (n_missing > 0) {
    cat(sprintf("⚠️  Warning: %d sites missing terminal nucleotide annotation\n", n_missing))
    vcf_annotated <- vcf_annotated[!is.na(Terminal_Nucleotide)]
  }
  
  cat(sprintf("Sites with annotations: %d\n\n", nrow(vcf_annotated)))
  
  # 4. Estimate gamma for each bin (pooling across all genes)
  cat("Estimating gamma for each positional bin...\\n")
  
  gradient_results <- vcf_annotated[, {
    
    # Select appropriate alpha/beta based on terminal nucleotide
    term_nuc <- Terminal_Nucleotide[1]
    
    if (term_nuc == "G") {
      alpha_use <- neutral_params$alpha_G
      beta_use <- neutral_params$beta_G
    } else {  # C
      alpha_use <- neutral_params$alpha_C
      beta_use <- neutral_params$beta_C
    }
    
    # Estimate gamma pooling all sites in this bin
    n_sites_bin <- .N
    
    if (n_sites_bin < 10) {
      # Not enough sites in this bin
      list(
        Bin_Center = mean(Relative_Position),
        Gamma = NA_real_,
        N_Sites = n_sites_bin,
        Mean_Freq_Pref = NA_real_
      )
    } else {
      # Estimate gamma
      gamma_est <- tryCatch({
        estimate_gamma_for_AA(
          counts = k,
          sample_sizes = n,
          alpha = alpha_use,
          beta = beta_use,
          S_interval = c(0, 50)
        )
      },
      error = function(e) NA_real_)
      
      list(
        Bin_Center = mean(Relative_Position),
        Gamma = gamma_est,
        N_Sites = n_sites_bin,
        Mean_Freq_Pref = mean(p, na.rm = TRUE)
      )
    }
    
  }, by = Position_Bin]
  
  # 5. Test for gradient (correlation between position and gamma)
  if (sum(!is.na(gradient_results$Gamma)) >= 3) {
    
    valid_data <- gradient_results[!is.na(Gamma)]
    
    cor_test <- cor.test(valid_data$Bin_Center, valid_data$Gamma,
                         method = "spearman")
    
    cat("\\n=== Gradient Analysis Results ===\\n\\n")
    cat(sprintf("Spearman correlation (Position vs Gamma): ρ = %.3f\\n",
                cor_test$estimate))
    cat(sprintf("P-value: %.4f\\n\\n", cor_test$p.value))
    
    if (cor_test$estimate < 0 && cor_test$p.value < 0.05) {
      cat("✓ SIGNIFICANT NEGATIVE GRADIENT: Gamma decreases toward 3' end\\n")
      cat("  → Consistent with gBGC (stronger selection at 5' end)\\n\\n")
    } else if (cor_test$estimate > 0 && cor_test$p.value < 0.05) {
      cat("✓ SIGNIFICANT POSITIVE GRADIENT: Gamma increases toward 3' end\\n")
      cat("  → NOT consistent with typical gBGC pattern\\n\\n")
    } else {
      cat("✗ NO SIGNIFICANT GRADIENT detected\\n")
      cat("  → Selection appears uniform along gene length\\n\\n")
    }
    
    # 6. Plot the gradient
    p <- ggplot(gradient_results[!is.na(Gamma)], 
                aes(x = Bin_Center, y = Gamma)) +
      geom_point(aes(size = N_Sites), alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "dashed") +
      labs(
        title = "Gamma Gradient Along Gene Length (gBGC Test)",
        subtitle = sprintf("Spearman ρ = %.3f, p = %.4f", 
                          cor_test$estimate, cor_test$p.value),
        x = "Relative Position (0 = 5' end, 1 = 3' end)",
        y = "Gamma (Selection Coefficient)",
        size = "Sites per Bin"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    print(p)
    
    # Save plot
    if (!dir.exists("./results")) dir.create("./results")
    ggsave("./results/gamma_gradient_gBGC_test.pdf", 
           plot = p, 
           width = 8, height = 6)
    
    cat("✓ Plot saved: ./results/gamma_gradient_gBGC_test.pdf\\n\\n")
    
    # Add correlation stats to output
    gradient_results[, Spearman_rho := cor_test$estimate]
    gradient_results[, Spearman_p := cor_test$p.value]
    
  } else {
    warning("Not enough bins with valid gamma estimates (need ≥3)")
  }
  
  return(gradient_results)
}
