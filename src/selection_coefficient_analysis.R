#' Selection Coefficient Analysis - Per Amino Acid Family (Enhanced)
#' 
#' @description Estimate population-scaled selection coefficients (S = 4Nes)
#' for codon usage bias using the mutation-selection-drift balance model
#' from Hershberg & Petrov (2008), calculated independently for each amino acid family
#' with robust validation and diagnostic checks
#' 
#' @reference Hershberg R, Petrov DA. 2008. Selection on codon bias. 
#' Annu Rev Genet. 42:287-299.
#' 
#' @author Luis J. Madrigal-Roca
#' @date November 16, 2025
#' _____________________________________________________________________________


#' Get codon families (synonymous codons per amino acid)
#'
#' @param genetic_code Named vector mapping codons to amino acids
#' @param exclude_single Exclude Met and Trp (single codon AA) (default: TRUE)
#' @return List where each element is a vector of codons for one amino acid
get_codon_families <- function(genetic_code, exclude_single = TRUE) {
  
  # Remove STOP codons
  sense_codons <- genetic_code[genetic_code != "STOP"]
  
  # Optionally exclude single-codon amino acids
  if (exclude_single) {
    single_codon_aa <- c("Met", "Trp")
    sense_codons <- sense_codons[!sense_codons %in% single_codon_aa]
  }
  
  # Group codons by amino acid
  families <- split(names(sense_codons), sense_codons)
  
  return(families)
}


#' Identify preferred codon for each amino acid family
#'
#' @param preferred_codons Character vector of preferred codon names
#' @param genetic_code Named vector mapping codons to amino acids
#' @return Named list: amino acid -> preferred codon(s) in that family
get_preferred_by_family <- function(preferred_codons, genetic_code) {
  
  families <- get_codon_families(genetic_code, exclude_single = TRUE)
  
  preferred_by_aa <- lapply(families, function(family_codons) {
    # Find which codons in this family are preferred
    preferred_in_family <- intersect(family_codons, preferred_codons)
    return(preferred_in_family)
  })
  
  # Remove amino acids with no preferred codons defined
  preferred_by_aa <- preferred_by_aa[sapply(preferred_by_aa, length) > 0]
  
  return(preferred_by_aa)
}


#' Calculate proportion of preferred codons in ONE amino acid family for ONE gene
#'
#' @param gene_codon_counts Named vector of codon counts for one gene
#' @param family_codons Character vector of all codons in this AA family
#' @param preferred_codons Character vector of preferred codons in this AA family
#' @param min_codons Minimum number of codons required (default: 5)
#' @return Proportion of preferred codons (0 to 1) for this family, or NA if insufficient data
calculate_family_preferred_proportion <- function(gene_codon_counts, 
                                                  family_codons, 
                                                  preferred_codons,
                                                  min_codons = 5) {
  
  # Get counts for this family only
  family_counts <- gene_codon_counts[names(gene_codon_counts) %in% family_codons]
  family_counts <- as.numeric(family_counts[!is.na(family_counts)])
  
  total_count <- sum(family_counts)
  
  # Require minimum codon count for reliable estimation
  if (length(family_counts) == 0 || total_count < min_codons) {
    return(NA)
  }
  
  # Count preferred vs total in this family
  preferred_counts <- gene_codon_counts[names(gene_codon_counts) %in% preferred_codons]
  preferred_counts <- as.numeric(preferred_counts[!is.na(preferred_counts)])
  preferred_count <- sum(preferred_counts, na.rm = TRUE)
  
  if (total_count == 0) {
    return(NA)
  }
  
  # Add pseudocount to avoid exact 0 or 1 (which causes ln(0) or ln(∞))
  # Use Laplace smoothing: (k + α) / (n + α*K) where α=0.5, K=2 (preferred/unpreferred)
  alpha <- 0.5
  K <- 2
  P_adjusted <- (preferred_count + alpha) / (total_count + alpha * K)
  
  return(P_adjusted)
}


#' Estimate mutation bias per amino acid family from low-expression genes
#'
#' @param codon_usage Data frame with codon counts (Gene_name + codon columns)
#' @param expression_data Data frame with Gene_name and expression level
#' @param preferred_by_family Named list from get_preferred_by_family()
#' @param genetic_code Named vector mapping codons to amino acids
#' @param low_expr_quantile Quantile to define "low expression" (default: 0.10)
#' @param min_genes_per_family Minimum genes required to estimate M (default: 20)
#' @param min_codons Minimum codons per gene per family (default: 5)
#' @return Named list: amino acid -> list(mutation_bias, P_neutral, n_genes, CI)
estimate_mutation_bias_by_family <- function(codon_usage, expression_data, 
                                             preferred_by_family, genetic_code,
                                             low_expr_quantile = 0.10,
                                             min_genes_per_family = 20,
                                             min_codons = 5) {
  
  cat("\n=== Estimating Mutation Bias PER Amino Acid Family ===\n")
  
  # Merge with expression data
  merged <- codon_usage |>
    left_join(expression_data, by = "Gene_name")
  
  # Define low-expression genes
  expr_threshold <- quantile(merged$Expression, probs = low_expr_quantile, na.rm = TRUE)
  low_expr_genes <- merged |> filter(Expression <= expr_threshold)
  
  cat(sprintf("Using bottom %.0f%% of genes (n=%d) with expression <= %.2f\n",
              low_expr_quantile * 100, nrow(low_expr_genes), expr_threshold))
  
  # Get codon families
  families <- get_codon_families(genetic_code, exclude_single = TRUE)
  
  # Calculate mutation bias for each amino acid family
  mutation_params_by_aa <- lapply(names(preferred_by_family), function(aa) {
    
    family_codons <- families[[aa]]
    preferred_codons <- preferred_by_family[[aa]]
    
    # Calculate P(preferred) for this family in each low-expression gene
    P_values <- apply(low_expr_genes, 1, function(gene_row) {
      # Remove Gene_name and Expression columns
      gene_counts <- gene_row[!names(gene_row) %in% c("Gene_name", "Expression")]
      calculate_family_preferred_proportion(gene_counts, family_codons, 
                                            preferred_codons, min_codons)
    })
    
    # Remove NA values
    P_values <- P_values[!is.na(P_values)]
    
    if (length(P_values) < min_genes_per_family) {
      cat(sprintf("  WARNING: %s has only %d genes (< %d minimum) - SKIPPING\n", 
                  aa, length(P_values), min_genes_per_family))
      return(NULL)
    }
    
    # Average P(preferred) across low-expression genes
    P_neutral <- mean(P_values, na.rm = TRUE)
    P_neutral_se <- sd(P_values, na.rm = TRUE) / sqrt(length(P_values))
    
    # Calculate mutation bias: M = (1 - P) / P
    mutation_bias <- (1 - P_neutral) / P_neutral
    
    # Bootstrap 95% CI for M
    n_boot <- 1000
    boot_P <- replicate(n_boot, {
      boot_sample <- sample(P_values, replace = TRUE)
      mean(boot_sample)
    })
    boot_M <- (1 - boot_P) / boot_P
    M_CI <- quantile(boot_M, c(0.025, 0.975))
    
    cat(sprintf("  %s (n=%d): P_neutral=%.4f±%.4f, M=%.4f [%.4f, %.4f]\n", 
                aa, length(P_values), P_neutral, P_neutral_se, 
                mutation_bias, M_CI[1], M_CI[2]))
    
    return(list(
      amino_acid = aa,
      mutation_bias = mutation_bias,
      M_CI_lower = M_CI[1],
      M_CI_upper = M_CI[2],
      P_neutral = P_neutral,
      P_neutral_se = P_neutral_se,
      n_genes = length(P_values),
      family_codons = family_codons,
      preferred_codons = preferred_codons
    ))
  })
  
  names(mutation_params_by_aa) <- names(preferred_by_family)
  
  # Remove NULL entries (families with insufficient data)
  mutation_params_by_aa <- mutation_params_by_aa[!sapply(mutation_params_by_aa, is.null)]
  
  cat(sprintf("\nSuccessfully estimated mutation bias for %d amino acid families\n", 
              length(mutation_params_by_aa)))
  
  return(mutation_params_by_aa)
}


#' Calculate S for one amino acid family in one gene
#'
#' @param P_preferred Proportion of preferred codons in this family
#' @param mutation_bias M for this amino acid family
#' @return S = 4Nes for this family
calculate_S_family <- function(P_preferred, mutation_bias) {
  
  if (is.na(P_preferred) || is.na(mutation_bias)) {
    return(NA)
  }
  
  # After pseudocount adjustment, P should be in (0, 1)
  if (P_preferred <= 0 || P_preferred >= 1) {
    return(NA)
  }
  
  # Hershberg & Petrov equation: S = ln(M * (P / (1-P)))
  odds_ratio <- P_preferred / (1 - P_preferred)
  S <- log(mutation_bias * odds_ratio)
  
  return(S)
}


#' Validate mutation bias estimates against genomic composition
#'
#' @param mutation_params List from estimate_mutation_bias_by_family()
#' @param codon_usage Data frame with codon counts
#' @param genetic_code Genetic code mapping
validate_mutation_bias <- function(mutation_params, codon_usage, genetic_code) {
  
  cat("\n=== Validating Mutation Bias Estimates ===\n")
  
  # Calculate genome-wide GC content at 3rd codon position
  families <- get_codon_families(genetic_code, exclude_single = TRUE)
  
  # Get all codon counts (excluding Gene_name)
  codon_cols <- names(codon_usage)[names(codon_usage) != "Gene_name"]
  total_codons <- colSums(codon_usage[, ..codon_cols], na.rm = TRUE)
  
  # Calculate GC3 (GC at 3rd position for 4-fold degenerate sites)
  four_fold_families <- c("Ala", "Arg", "Gly", "Leu", "Pro", "Ser", "Thr", "Val")
  
  gc3_counts <- sapply(four_fold_families, function(aa) {
    if (!aa %in% names(families)) return(c(GC = 0, AT = 0))
    
    family_codons <- families[[aa]]
    gc3_codons <- family_codons[substr(family_codons, 3, 3) %in% c("G", "C")]
    at3_codons <- family_codons[substr(family_codons, 3, 3) %in% c("A", "T")]
    
    gc_count <- sum(total_codons[names(total_codons) %in% gc3_codons], na.rm = TRUE)
    at_count <- sum(total_codons[names(total_codons) %in% at3_codons], na.rm = TRUE)
    
    return(c(GC = gc_count, AT = at_count))
  })
  
  total_gc3 <- sum(gc3_counts["GC", ])
  total_at3 <- sum(gc3_counts["AT", ])
  genome_gc3 <- total_gc3 / (total_gc3 + total_at3)
  
  cat(sprintf("\nGenomic GC3 content (4-fold degenerate): %.4f\n", genome_gc3))
  cat(sprintf("Genomic AT3 content: %.4f\n", 1 - genome_gc3))
  
  # Compare with mutation bias estimates
  cat("\nMutation bias consistency check:\n")
  cat("If M > 1: mutation favors unpreferred (expected if preferred = GC-ending)\n")
  cat("If M < 1: mutation favors preferred (expected if preferred = AT-ending)\n\n")
  
  for (aa in names(mutation_params)) {
    params <- mutation_params[[aa]]
    preferred <- params$preferred_codons
    
    # Check if preferred codons are GC-rich or AT-rich at 3rd position
    pref_gc3 <- sum(substr(preferred, 3, 3) %in% c("G", "C"))
    pref_at3 <- sum(substr(preferred, 3, 3) %in% c("A", "T"))
    
    expectation <- ifelse(pref_gc3 > pref_at3, "M > 1", "M < 1")
    observed <- ifelse(params$mutation_bias > 1, "M > 1", "M < 1")
    match <- ifelse(expectation == observed, "✓", "✗")
    
    cat(sprintf("  %s: M=%.4f, preferred=%s, expect %s, observed %s %s\n",
                aa, params$mutation_bias, 
                paste(preferred, collapse=","),
                expectation, observed, match))
  }
}


#' Test sensitivity to low-expression quantile
#'
#' @param codon_usage Data frame with codon counts
#' @param expression_data Data frame with expression
#' @param preferred_by_family Preferred codons by family
#' @param genetic_code Genetic code
#' @param quantiles Vector of quantiles to test
test_quantile_sensitivity <- function(codon_usage, expression_data, 
                                      preferred_by_family, genetic_code,
                                      quantiles = c(0.05, 0.10, 0.15, 0.20)) {
  
  cat("\n=== Testing Sensitivity to Low-Expression Quantile ===\n")
  
  results_list <- lapply(quantiles, function(q) {
    cat(sprintf("\n--- Testing quantile = %.2f ---\n", q))
    
    params <- estimate_mutation_bias_by_family(
      codon_usage, expression_data, preferred_by_family, genetic_code,
      low_expr_quantile = q, min_genes_per_family = 10
    )
    
    # Extract M values
    M_values <- sapply(params, function(p) p$mutation_bias)
    
    data.frame(
      quantile = q,
      n_families = length(M_values),
      mean_M = mean(M_values, na.rm = TRUE),
      median_M = median(M_values, na.rm = TRUE),
      sd_M = sd(M_values, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  
  sensitivity_df <- do.call(rbind, results_list)
  
  cat("\n=== Quantile Sensitivity Summary ===\n")
  print(sensitivity_df, digits = 4)
  
  cat("\nInterpretation: M should be relatively stable across quantiles.\n")
  cat("Large variations suggest the assumption of s≈0 may not hold.\n")
  
  return(sensitivity_df)
}


#' Calculate selection coefficients per family for all genes
#'
#' @param codon_usage Data frame with Gene_name and codon counts
#' @param expression_data Data frame with Gene_name and Expression
#' @param preferred_codons Character vector of preferred codon names (w=1.0)
#' @param genetic_code Named vector mapping codons to amino acids
#' @param low_expr_quantile Quantile for estimating mutation bias (default: 0.10)
#' @param min_codons Minimum codons per gene per family (default: 5)
#' @param run_diagnostics Run validation and sensitivity tests (default: TRUE)
#' @return List with: gene_level, family_level, mutation_params, diagnostics
calculate_selection_coefficients_by_family <- function(codon_usage, expression_data, 
                                                       preferred_codons, genetic_code,
                                                       low_expr_quantile = 0.10,
                                                       min_codons = 5,
                                                       run_diagnostics = TRUE) {
  
  cat("\n=== Calculating Selection Coefficients Per Amino Acid Family ===\n")
  
  # Step 1: Get preferred codons organized by family
  preferred_by_family <- get_preferred_by_family(preferred_codons, genetic_code)
  
  cat(sprintf("\nAnalyzing %d amino acid families:\n", length(preferred_by_family)))
  for (aa in names(preferred_by_family)) {
    cat(sprintf("  %s: %s (preferred: %s)\n", 
                aa, 
                paste(get_codon_families(genetic_code)[[aa]], collapse=", "),
                paste(preferred_by_family[[aa]], collapse=", ")))
  }
  
  # Step 2: Run diagnostics if requested
  diagnostics <- list()
  
  if (run_diagnostics) {
    cat("\n" , rep("=", 70), "\n", sep="")
    cat("RUNNING DIAGNOSTIC TESTS\n")
    cat(rep("=", 70), "\n", sep="")
    
    # Test quantile sensitivity
    diagnostics$quantile_sensitivity <- test_quantile_sensitivity(
      codon_usage, expression_data, preferred_by_family, genetic_code
    )
  }
  
  # Step 3: Estimate mutation bias per family (using specified quantile)
  mutation_params <- estimate_mutation_bias_by_family(
    codon_usage, expression_data, preferred_by_family, genetic_code,
    low_expr_quantile, min_codons = min_codons
  )
  
  # Validate mutation bias
  if (run_diagnostics) {
    validate_mutation_bias(mutation_params, codon_usage, genetic_code)
  }
  
  # Step 4: Calculate S for each gene for each amino acid family
  cat("\n=== Calculating S per gene per amino acid family ===\n")
  
  families <- get_codon_families(genetic_code, exclude_single = TRUE)
  
  # Initialize results matrix: rows = genes, columns = amino acids
  S_matrix <- matrix(NA, nrow = nrow(codon_usage), ncol = length(mutation_params))
  rownames(S_matrix) <- codon_usage$Gene_name
  colnames(S_matrix) <- names(mutation_params)
  
  # Calculate S for each gene and each family
  pb <- txtProgressBar(min = 0, max = length(mutation_params), style = 3)
  
  for (idx in seq_along(mutation_params)) {
    aa <- names(mutation_params)[idx]
    
    family_codons <- mutation_params[[aa]]$family_codons
    preferred_codons_aa <- mutation_params[[aa]]$preferred_codons
    M <- mutation_params[[aa]]$mutation_bias
    
    # For each gene, calculate P and S for this family
    for (i in 1:nrow(codon_usage)) {
      gene_counts <- codon_usage[i, !names(codon_usage) %in% "Gene_name"]
      
      P_preferred <- calculate_family_preferred_proportion(
        gene_counts, family_codons, preferred_codons_aa, min_codons
      )
      
      S_matrix[i, aa] <- calculate_S_family(P_preferred, M)
    }
    
    setTxtProgressBar(pb, idx)
  }
  close(pb)
  
  # Step 5: Calculate weighted average S across families for each gene
  cat("\n\nCalculating weighted average S (weighted by AA frequency) for each gene...\n")
  
  # Count synonymous codons per family per gene for weighting
  weight_matrix <- matrix(0, nrow = nrow(codon_usage), ncol = length(mutation_params))
  rownames(weight_matrix) <- codon_usage$Gene_name
  colnames(weight_matrix) <- names(mutation_params)
  
  for (aa in names(mutation_params)) {
    family_codons <- mutation_params[[aa]]$family_codons
    
    for (i in 1:nrow(codon_usage)) {
      gene_counts <- codon_usage[i, !names(codon_usage) %in% "Gene_name"]
      family_count <- sum(as.numeric(gene_counts[names(gene_counts) %in% family_codons]), na.rm = TRUE)
      weight_matrix[i, aa] <- family_count
    }
  }
  
  # Calculate weighted average S
  gene_avg_S <- sapply(1:nrow(S_matrix), function(i) {
    s_vals <- S_matrix[i, ]
    weights <- weight_matrix[i, ]
    
    # Only use families with valid S and non-zero weight
    valid_idx <- !is.na(s_vals) & weights > 0
    
    if (sum(valid_idx) == 0) return(NA)
    
    weighted.mean(s_vals[valid_idx], weights[valid_idx])
  })
  
  # Also calculate unweighted for comparison
  gene_avg_S_unweighted <- rowMeans(S_matrix, na.rm = TRUE)
  gene_sd_S <- apply(S_matrix, 1, sd, na.rm = TRUE)
  gene_n_families <- rowSums(!is.na(S_matrix))
  
  # Bootstrap CI for each gene (optional - computationally expensive)
  # For now, use SE approximation: SE ≈ SD / sqrt(n_families)
  gene_se_S <- gene_sd_S / sqrt(gene_n_families)
  
  # Compile gene-level results
  gene_results <- data.frame(
    Gene_name = codon_usage$Gene_name,
    S_weighted = gene_avg_S,
    S_unweighted = gene_avg_S_unweighted,
    S_sd = gene_sd_S,
    S_se = gene_se_S,
    S_median = apply(S_matrix, 1, median, na.rm = TRUE),
    n_families = gene_n_families,
    stringsAsFactors = FALSE
  )
  
  # Merge with expression
  gene_results <- gene_results |>
    left_join(expression_data, by = "Gene_name")
  
  # Step 6: Summary statistics
  cat("\n=== Gene-Level S (Weighted by AA Composition) ===\n")
  cat(sprintf("Mean S (weighted): %.4f\n", mean(gene_results$S_weighted, na.rm = TRUE)))
  cat(sprintf("Mean S (unweighted): %.4f\n", mean(gene_results$S_unweighted, na.rm = TRUE)))
  cat(sprintf("Median S (weighted): %.4f\n", median(gene_results$S_weighted, na.rm = TRUE)))
  cat(sprintf("SD S: %.4f\n", sd(gene_results$S_weighted, na.rm = TRUE)))
  cat(sprintf("Range: %.4f to %.4f\n", 
              min(gene_results$S_weighted, na.rm = TRUE), 
              max(gene_results$S_weighted, na.rm = TRUE)))
  
  # Difference between weighted and unweighted
  weight_diff <- gene_results$S_weighted - gene_results$S_unweighted
  cat(sprintf("\nWeighted vs Unweighted difference:\n"))
  cat(sprintf("  Mean difference: %.4f\n", mean(weight_diff, na.rm = TRUE)))
  cat(sprintf("  SD difference: %.4f\n", sd(weight_diff, na.rm = TRUE)))
  cat(sprintf("  Max |difference|: %.4f\n", max(abs(weight_diff), na.rm = TRUE)))
  
  if (mean(abs(weight_diff), na.rm = TRUE) > 0.1) {
    cat("  → Large differences suggest amino acid composition effects are important!\n")
  } else {
    cat("  → Small differences suggest composition effects are minimal.\n")
  }
  
  # Count positive vs negative selection
  n_positive <- sum(gene_results$S_weighted > 0, na.rm = TRUE)
  n_negative <- sum(gene_results$S_weighted < 0, na.rm = TRUE)
  n_total <- sum(!is.na(gene_results$S_weighted))
  
  cat(sprintf("\nGenes with positive selection (S > 0): %d / %d (%.1f%%)\n",
              n_positive, n_total, 100 * n_positive / n_total))
  cat(sprintf("Genes with negative selection (S < 0): %d / %d (%.1f%%)\n",
              n_negative, n_total, 100 * n_negative / n_total))
  
  # Correlation with expression
  if ("Expression" %in% names(gene_results)) {
    complete_data <- gene_results[!is.na(gene_results$S_weighted) & !is.na(gene_results$Expression), ]
    
    cor_test_spearman <- cor.test(complete_data$S_weighted, log2(complete_data$Expression + 1), 
                                  method = "spearman", exact = FALSE)
    cor_test_pearson <- cor.test(complete_data$S_weighted, log2(complete_data$Expression + 1), 
                                 method = "pearson")
    
    cat(sprintf("\nSpearman correlation: S vs log2(Expression+1) = %.4f\n", 
                cor_test_spearman$estimate))
    cat(sprintf("  p-value: %.2e (n = %d genes)\n", 
                cor_test_spearman$p.value, nrow(complete_data)))
    
    cat(sprintf("Pearson correlation: r = %.4f, p = %.2e\n",
                cor_test_pearson$estimate, cor_test_pearson$p.value))
    
    if (cor_test_spearman$p.value < 0.001) cat("  *** Highly significant\n")
    else if (cor_test_spearman$p.value < 0.01) cat("  ** Very significant\n")
    else if (cor_test_spearman$p.value < 0.05) cat("  * Significant\n")
  }
  
  # Family-level statistics
  cat("\n=== Family-Level S Statistics ===\n")
  for (aa in colnames(S_matrix)) {
    S_vals <- S_matrix[, aa]
    params <- mutation_params[[aa]]
    
    cat(sprintf("%s: M=%.3f [%.3f,%.3f], S: mean=%.3f, median=%.3f, SD=%.3f (n=%d)\n",
                aa, 
                params$mutation_bias,
                params$M_CI_lower,
                params$M_CI_upper,
                mean(S_vals, na.rm = TRUE),
                median(S_vals, na.rm = TRUE),
                sd(S_vals, na.rm = TRUE),
                sum(!is.na(S_vals))))
  }
  
  # Check for outliers
  cat("\n=== Outlier Detection ===\n")
  S_values_all <- as.vector(S_matrix)
  S_values_all <- S_values_all[!is.na(S_values_all)]
  
  Q1 <- quantile(S_values_all, 0.25)
  Q3 <- quantile(S_values_all, 0.75)
  IQR <- Q3 - Q1
  outlier_threshold <- Q3 + 3 * IQR
  
  n_outliers <- sum(S_values_all > outlier_threshold | S_values_all < (Q1 - 3*IQR))
  cat(sprintf("Extreme outliers (beyond Q1/Q3 ± 3*IQR): %d / %d (%.2f%%)\n",
              n_outliers, length(S_values_all), 100*n_outliers/length(S_values_all)))
  
  if (n_outliers > 0) {
    extreme_idx <- which(gene_results$S_weighted > outlier_threshold | 
                           gene_results$S_weighted < (Q1 - 3*IQR))
    cat("\nTop outlier genes:\n")
    outlier_genes <- gene_results[extreme_idx, c("Gene_name", "S_weighted", "Expression")]
    outlier_genes <- outlier_genes[order(-abs(outlier_genes$S_weighted)), ]
    print(head(outlier_genes, 10))
  }
  
  # Convert S_matrix to long format for easier analysis
  family_results <- as.data.frame(S_matrix) |>
    tibble::rownames_to_column("Gene_name") |>
    tidyr::pivot_longer(cols = -Gene_name, 
                        names_to = "Amino_acid", 
                        values_to = "S") |>
    left_join(expression_data, by = "Gene_name")
  
  # Store mutation parameters and matrices as attributes
  attr(gene_results, "mutation_params") <- mutation_params
  attr(gene_results, "S_matrix") <- S_matrix
  attr(gene_results, "weight_matrix") <- weight_matrix
  
  return(list(
    gene_level = gene_results,
    family_level = family_results,
    mutation_params = mutation_params,
    diagnostics = diagnostics
  ))
}


#' Plot S distribution per amino acid family
#'
#' @param results List from calculate_selection_coefficients_by_family()
#' @param output_file Path to save plot
plot_S_by_family <- function(results, 
                             output_file = "./results/S_by_amino_acid_family.pdf") {
  
  library(ggplot2)
  
  plot_data <- results$family_level |> filter(!is.na(S))
  
  # Add sample sizes to labels
  n_per_family <- plot_data |>
    group_by(Amino_acid) |>
    summarise(n = sum(!is.na(S)))
  
  plot_data <- plot_data |>
    left_join(n_per_family, by = "Amino_acid") |>
    mutate(AA_label = paste0(Amino_acid, "\n(n=", n, ")"))
  
  p <- ggplot(plot_data, aes(x = reorder(AA_label, S, FUN = median), y = S, fill = Amino_acid)) +
    geom_boxplot(outlier.alpha = 0.3, outlier.size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(
      title = "Selection Coefficients by Amino Acid Family",
      subtitle = "Ordered by median S value",
      x = "Amino Acid",
      y = "S = 4Nes"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11)
    )
  
  ggsave(output_file, p, width = 12, height = 6)
  cat(sprintf("\nPlot saved to: %s\n", output_file))
  
  return(p)
}


#' Plot gene-level S vs expression
#'
#' @param results List from calculate_selection_coefficients_by_family()
#' @param output_file Path to save plot
plot_gene_S_vs_expression <- function(results,
                                      output_file = "./results/gene_S_vs_expression.pdf") {
  
  library(ggplot2)
  
  plot_data <- results$gene_level |>
    filter(!is.na(S_weighted), !is.na(Expression), Expression > 0)
  
  cor_val <- cor(plot_data$S_weighted, log2(plot_data$Expression + 1), method = "spearman")
  cor_test <- cor.test(plot_data$S_weighted, log2(plot_data$Expression + 1), 
                       method = "spearman", exact = FALSE)
  
  p_label <- ifelse(cor_test$p.value < 0.001, "p < 0.001",
                    sprintf("p = %.3f", cor_test$p.value))
  
  p <- ggplot(plot_data, aes(x = log2(Expression + 1), y = S_weighted)) +
    geom_point(alpha = 0.3, size = 1.5, color = "steelblue") +
    geom_smooth(method = "loess", color = "red", se = TRUE, linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    labs(
      title = "Gene-Level Selection vs Expression",
      subtitle = sprintf("S weighted by amino acid composition | Spearman ρ = %.3f (%s)", 
                         cor_val, p_label),
      x = "log2(Expression + 1)",
      y = "Weighted Mean S = 4Nes"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11)
    )
  
  ggsave(output_file, p, width = 8, height = 6)
  cat(sprintf("\nPlot saved to: %s\n", output_file))
  
  return(p)
}


#' Plot weighted vs unweighted S comparison
#'
#' @param results List from calculate_selection_coefficients_by_family()
#' @param output_file Path to save plot
plot_weighted_vs_unweighted <- function(results,
                                        output_file = "./results/S_weighted_vs_unweighted.pdf") {
  
  library(ggplot2)
  
  plot_data <- results$gene_level |>
    filter(!is.na(S_weighted), !is.na(S_unweighted))
  
  cor_val <- cor(plot_data$S_weighted, plot_data$S_unweighted, method = "pearson")
  
  p <- ggplot(plot_data, aes(x = S_unweighted, y = S_weighted)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
    labs(
      title = "Weighted vs Unweighted S",
      subtitle = sprintf("Pearson r = %.4f | Red line: y = x", cor_val),
      x = "S (unweighted mean across families)",
      y = "S (weighted by AA composition)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14)
    )
  
  ggsave(output_file, p, width = 7, height = 7)
  cat(sprintf("\nPlot saved to: %s\n", output_file))
  
  return(p)
}


#' Plot mutation bias validation
#'
#' @param results List from calculate_selection_coefficients_by_family()
#' @param output_file Path to save plot
plot_mutation_bias_validation <- function(results,
                                          output_file = "./results/mutation_bias_validation.pdf") {
  
  library(ggplot2)
  
  # Extract mutation parameters
  M_data <- lapply(results$mutation_params, function(p) {
    data.frame(
      Amino_acid = p$amino_acid,
      M = p$mutation_bias,
      M_lower = p$M_CI_lower,
      M_upper = p$M_CI_upper,
      P_neutral = p$P_neutral,
      stringsAsFactors = FALSE
    )
  })
  
  M_df <- do.call(rbind, M_data)
  
  p <- ggplot(M_df, aes(x = reorder(Amino_acid, M), y = M)) +
    geom_point(size = 3, color = "steelblue") +
    geom_errorbar(aes(ymin = M_lower, ymax = M_upper), width = 0.3, color = "steelblue") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(
      title = "Mutation Bias by Amino Acid Family",
      subtitle = "M = μ_unpreferred / μ_preferred | Error bars: 95% bootstrap CI",
      x = "Amino Acid",
      y = "Mutation Bias (M)"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10)
    ) +
    annotate("text", x = 1, y = max(M_df$M_upper) * 0.95, 
             label = "M > 1: mutation favors unpreferred\nM < 1: mutation favors preferred",
             hjust = 0, size = 3.5, color = "gray30")
  
  ggsave(output_file, p, width = 10, height = 6)
  cat(sprintf("\nPlot saved to: %s\n", output_file))
  
  return(p)
}


#' Plot S distribution histogram
#'
#' @param results List from calculate_selection_coefficients_by_family()
#' @param output_file Path to save plot
plot_S_distribution <- function(results,
                                output_file = "./results/S_distribution.pdf") {
  
  library(ggplot2)
  
  plot_data <- results$gene_level |> filter(!is.na(S_weighted))
  
  mean_S <- mean(plot_data$S_weighted)
  median_S <- median(plot_data$S_weighted)
  
  p <- ggplot(plot_data, aes(x = S_weighted)) +
    geom_histogram(bins = 60, fill = "steelblue", color = "white", alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    geom_vline(xintercept = median_S, linetype = "dashed", color = "darkgreen", linewidth = 1) +
    geom_vline(xintercept = mean_S, linetype = "dotted", color = "darkblue", linewidth = 1) +
    labs(
      title = "Distribution of Gene-Level Selection Coefficients",
      subtitle = sprintf("Median S = %.3f | Mean S = %.3f", median_S, mean_S),
      x = "S = 4Nes (weighted by AA composition)",
      y = "Number of genes"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11)
    ) +
    annotate("text", x = median_S, y = Inf, label = "Median", 
             vjust = 2, hjust = -0.1, size = 3.5, color = "darkgreen") +
    annotate("text", x = mean_S, y = Inf, label = "Mean", 
             vjust = 3.5, hjust = -0.1, size = 3.5, color = "darkblue")
  
  ggsave(output_file, p, width = 9, height = 6)
  cat(sprintf("\nPlot saved to: %s\n", output_file))
  
  return(p)
}


#' Plot correlation between family-level S values
#'
#' @param results List from calculate_selection_coefficients_by_family()
#' @param output_file Path to save plot
plot_family_correlations <- function(results,
                                     output_file = "./results/family_S_correlations.pdf") {
  
  library(ggplot2)
  library(reshape2)
  
  S_matrix <- attr(results$gene_level, "S_matrix")
  
  # Calculate correlation matrix
  cor_matrix <- cor(S_matrix, use = "pairwise.complete.obs", method = "spearman")
  
  # Melt for ggplot
  cor_melted <- melt(cor_matrix)
  
  p <- ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limit = c(-1, 1),
                         name = "Spearman ρ") +
    labs(
      title = "Correlation Between Family-Level S Values",
      subtitle = "Do amino acid families experience similar selection?",
      x = "", y = ""
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10)
    )
  
  ggsave(output_file, p, width = 10, height = 9)
  cat(sprintf("\nPlot saved to: %s\n", output_file))
  
  # Print interpretation
  avg_cor <- mean(cor_matrix[upper.tri(cor_matrix)])
  cat(sprintf("\nAverage pairwise correlation: %.3f\n", avg_cor))
  
  if (avg_cor > 0.5) {
    cat("  → Strong positive correlations suggest shared selective pressures\n")
  } else if (avg_cor > 0.2) {
    cat("  → Moderate correlations suggest some shared selective pressures\n")
  } else {
    cat("  → Weak correlations suggest independent selection on families\n")
  }
  
  return(p)
}


#' Create comprehensive diagnostic report
#'
#' @param results List from calculate_selection_coefficients_by_family()
#' @param output_dir Directory to save all plots and reports
create_diagnostic_report <- function(results, output_dir = "./results") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("\n=== Generating Comprehensive Diagnostic Report ===\n")
  
  # Generate all plots
  plot_S_by_family(results, file.path(output_dir, "S_by_amino_acid_family.pdf"))
  plot_gene_S_vs_expression(results, file.path(output_dir, "gene_S_vs_expression.pdf"))
  plot_weighted_vs_unweighted(results, file.path(output_dir, "S_weighted_vs_unweighted.pdf"))
  plot_mutation_bias_validation(results, file.path(output_dir, "mutation_bias_validation.pdf"))
  plot_S_distribution(results, file.path(output_dir, "S_distribution.pdf"))
  plot_family_correlations(results, file.path(output_dir, "family_S_correlations.pdf"))
  
  # Save data tables
  write.csv(results$gene_level, file.path(output_dir, "gene_level_selection.csv"), 
            row.names = FALSE)
  write.csv(results$family_level, file.path(output_dir, "family_level_selection.csv"), 
            row.names = FALSE)
  
  # Save mutation parameters
  M_data <- lapply(results$mutation_params, function(p) {
    data.frame(
      Amino_acid = p$amino_acid,
      Mutation_bias = p$mutation_bias,
      M_CI_lower = p$M_CI_lower,
      M_CI_upper = p$M_CI_upper,
      P_neutral = p$P_neutral,
      P_neutral_SE = p$P_neutral_se,
      N_genes = p$n_genes,
      Preferred_codons = paste(p$preferred_codons, collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  M_df <- do.call(rbind, M_data)
  write.csv(M_df, file.path(output_dir, "mutation_bias_estimates.csv"), row.names = FALSE)
  
  # Save diagnostic results if available
  if (!is.null(results$diagnostics$quantile_sensitivity)) {
    write.csv(results$diagnostics$quantile_sensitivity, 
              file.path(output_dir, "quantile_sensitivity.csv"), 
              row.names = FALSE)
  }
  
  # Create summary text report
  sink(file.path(output_dir, "analysis_summary.txt"))
  
  cat("=" , rep("=", 70), "\n", sep="")
  cat("SELECTION COEFFICIENT ANALYSIS - SUMMARY REPORT\n")
  cat(rep("=", 70), "\n", sep="")
  cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  cat("DATASET SUMMARY\n")
  cat(rep("-", 70), "\n", sep="")
  cat("Total genes analyzed:", nrow(results$gene_level), "\n")
  cat("Amino acid families:", length(results$mutation_params), "\n")
  cat("Genes with valid S:", sum(!is.na(results$gene_level$S_weighted)), "\n\n")
  
  cat("GENE-LEVEL SELECTION (S = 4Nes)\n")
  cat(rep("-", 70), "\n", sep="")
  cat("Mean S (weighted):", sprintf("%.4f", mean(results$gene_level$S_weighted, na.rm = TRUE)), "\n")
  cat("Median S (weighted):", sprintf("%.4f", median(results$gene_level$S_weighted, na.rm = TRUE)), "\n")
  cat("SD S:", sprintf("%.4f", sd(results$gene_level$S_weighted, na.rm = TRUE)), "\n")
  cat("Range:", sprintf("[%.4f, %.4f]", 
                        min(results$gene_level$S_weighted, na.rm = TRUE),
                        max(results$gene_level$S_weighted, na.rm = TRUE)), "\n\n")
  
  cat("Positive selection (S > 0):", 
      sprintf("%d genes (%.1f%%)", 
              sum(results$gene_level$S_weighted > 0, na.rm = TRUE),
              100 * mean(results$gene_level$S_weighted > 0, na.rm = TRUE)), "\n")
  cat("Negative selection (S < 0):",
      sprintf("%d genes (%.1f%%)",
              sum(results$gene_level$S_weighted < 0, na.rm = TRUE),
              100 * mean(results$gene_level$S_weighted < 0, na.rm = TRUE)), "\n\n")
  
  if ("Expression" %in% names(results$gene_level)) {
    cor_test <- cor.test(results$gene_level$S_weighted, 
                         log2(results$gene_level$Expression + 1),
                         method = "spearman", exact = FALSE)
    cat("Correlation with expression:\n")
    cat("  Spearman ρ =", sprintf("%.4f", cor_test$estimate), "\n")
    cat("  p-value =", sprintf("%.2e", cor_test$p.value), "\n\n")
  }
  
  cat("FAMILY-LEVEL MUTATION BIAS\n")
  cat(rep("-", 70), "\n", sep="")
  for (aa in names(results$mutation_params)) {
    p <- results$mutation_params[[aa]]
    cat(sprintf("%-10s M = %.4f [%.4f, %.4f], P_neutral = %.4f (n=%d genes)\n",
                paste0(aa, ":"), p$mutation_bias, p$M_CI_lower, p$M_CI_upper,
                p$P_neutral, p$n_genes))
  }
  
  cat("\n", rep("=", 70), "\n", sep="")
  cat("All results saved to:", output_dir, "\n")
  cat(rep("=", 70), "\n", sep="")
  
  sink()
  
  cat("\n=== Diagnostic Report Complete ===\n")
  cat("Files saved to:", output_dir, "\n")
  cat("  - 6 diagnostic plots (PDF)\n")
  cat("  - 4 data tables (CSV)\n")
  cat("  - 1 summary report (TXT)\n")
}


#' Main wrapper function with all diagnostics
#'
#' @param codon_usage Your codon usage data frame
#' @param expression_data Your expression data frame
#' @param preferred_codons Vector of preferred codons
#' @param genetic_code Standard genetic code
#' @param output_dir Output directory (default: "./results")
#' @param low_expr_quantile Quantile for low expression (default: 0.10)
#' @param run_diagnostics Run all validation tests (default: TRUE)
run_selection_analysis <- function(codon_usage, expression_data, 
                                   preferred_codons, genetic_code,
                                   output_dir = "./results",
                                   low_expr_quantile = 0.10,
                                   run_diagnostics = TRUE) {
  
  cat("\n")
  cat(rep("=", 80), "\n", sep="")
  cat("SELECTION COEFFICIENT ANALYSIS - PER AMINO ACID FAMILY\n")
  cat(rep("=", 80), "\n", sep="")
  cat("\nStarting analysis with enhanced validation and diagnostics...\n")
  
  # Run the main analysis
  results <- calculate_selection_coefficients_by_family(
    codon_usage = codon_usage,
    expression_data = expression_data,
    preferred_codons = preferred_codons,
    genetic_code = genetic_code,
    low_expr_quantile = low_expr_quantile,
    min_codons = 5,
    run_diagnostics = run_diagnostics
  )
  
  # Generate comprehensive report
  create_diagnostic_report(results, output_dir)
  
  cat("\n")
  cat(rep("=", 80), "\n", sep="")
  cat("ANALYSIS COMPLETE\n")
  cat(rep("=", 80), "\n", sep="")
  
  return(results)
}