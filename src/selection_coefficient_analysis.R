#' Selection Coefficient Analysis
#' 
#' @description Estimate population-scaled selection coefficients (S = 4Nes)
#' for codon usage bias using the mutation-selection-drift balance model
#' from Hershberg & Petrov (2008) Molecular Biology and Evolution
#' 
#' @reference Hershberg R, Petrov DA. 2008. Selection on codon bias. 
#' Annu Rev Genet. 42:287-299.
#' 
#' @author Luis J. Madrigal-Roca
#' @date November 1, 2025
#' _____________________________________________________________________________

#' Calculate proportion of preferred codons in a gene
#'
#' @param gene_codon_counts Named vector of codon counts for one gene
#' @param preferred_codons Character vector of preferred codon names
#' @param genetic_code Named vector mapping codons to amino acids
#' @return Proportion of preferred codons (0 to 1)
calculate_preferred_proportion <- function(gene_codon_counts, preferred_codons, genetic_code) {
  
  # Remove Gene_name if present
  if ("Gene_name" %in% names(gene_codon_counts)) {
    gene_codon_counts <- gene_codon_counts[names(gene_codon_counts) != "Gene_name"]
  }
  
  # Convert to numeric
  codon_names <- names(gene_codon_counts)
  gene_codon_counts <- as.numeric(gene_codon_counts)
  names(gene_codon_counts) <- codon_names
  
  # Get sense codons (exclude STOP, Met, Trp)
  sense_codons <- names(genetic_code)[genetic_code != "STOP"]
  single_codon_aa <- c("ATG", "TGG")
  sense_codons <- sense_codons[!sense_codons %in% single_codon_aa]
  
  # Filter to sense codons only
  gene_counts <- gene_codon_counts[names(gene_codon_counts) %in% sense_codons]
  gene_counts <- gene_counts[!is.na(gene_counts)]
  
  if (sum(gene_counts) == 0) {
    return(NA)
  }
  
  # Count preferred vs total
  preferred_count <- sum(gene_counts[names(gene_counts) %in% preferred_codons], na.rm = TRUE)
  total_count <- sum(gene_counts, na.rm = TRUE)
  
  return(preferred_count / total_count)
}


#' Estimate mutation bias from low-expression genes
#'
#' @param codon_usage Data frame with codon counts (Gene_name + codon columns)
#' @param expression_data Data frame with Gene_name and expression level
#' @param preferred_codons Character vector of preferred codon names
#' @param genetic_code Named vector mapping codons to amino acids
#' @param low_expr_quantile Quantile to define "low expression" (default: 0.10 = bottom 10%)
#' @return List with mutation_bias (M) and P_neutral
estimate_mutation_bias <- function(codon_usage, expression_data, preferred_codons, 
                                   genetic_code, low_expr_quantile = 0.10) {
  
  cat("\n=== Estimating Mutation Bias from Low-Expression Genes ===\n")
  
  # Merge with expression data
  merged <- codon_usage |>
    left_join(expression_data, by = "Gene_name")
  
  # Define low-expression genes (bottom X%)
  expr_threshold <- quantile(merged$Expression, probs = low_expr_quantile, na.rm = TRUE)
  low_expr_genes <- merged |> filter(Expression <= expr_threshold)
  
  cat(sprintf("Using bottom %.0f%% of genes (n=%d) with expression <= %.2f\n",
              low_expr_quantile * 100, nrow(low_expr_genes), expr_threshold))
  
  # Calculate P(preferred) for each low-expression gene
  P_values <- apply(low_expr_genes, 1, function(gene_row) {
    calculate_preferred_proportion(gene_row, preferred_codons, genetic_code)
  })
  
  # Average across all low-expression genes (assumes s≈0, so P reflects mutation only)
  P_neutral <- mean(P_values, na.rm = TRUE)
  
  # Calculate mutation bias: M = μ_p/μ_u = (1 - P_neutral) / P_neutral
  mutation_bias <- (1 - P_neutral) / P_neutral
  
  cat(sprintf("\nP(preferred) in low-expression genes: %.4f\n", P_neutral))
  cat(sprintf("Estimated mutation bias (M = μ_p/μ_u): %.4f\n", mutation_bias))
  
  if (mutation_bias > 1) {
    cat("  → Mutation favors unpreferred codons (AT bias in mutation)\n")
  } else {
    cat("  → Mutation favors preferred codons (GC bias in mutation)\n")
  }
  
  return(list(
    mutation_bias = mutation_bias,
    P_neutral = P_neutral,
    low_expr_threshold = expr_threshold,
    n_genes = nrow(low_expr_genes)
  ))
}


#' Calculate population-scaled selection coefficient (S = 4Nes) for a gene
#'
#' @param P_preferred Proportion of preferred codons in the gene
#' @param mutation_bias M = μ_p/μ_u from estimate_mutation_bias()
#' @return S = 4Nes (population-scaled selection coefficient)
calculate_S <- function(P_preferred, mutation_bias) {
  
  if (is.na(P_preferred) || P_preferred <= 0 || P_preferred >= 1) {
    return(NA)
  }
  
  # From Hershberg & Petrov equation:
  # P(preferred) = 1 / (1 + M * exp(-S))
  # Solving for S:
  # S = ln(M * (P / (1-P)))
  
  odds_ratio <- P_preferred / (1 - P_preferred)
  S <- log(mutation_bias * odds_ratio)
  
  return(S)
}


#' Calculate selection coefficients for all genes
#'
#' @param codon_usage Data frame with Gene_name and codon counts
#' @param expression_data Data frame with Gene_name and Expression
#' @param preferred_codons Character vector of preferred codon names (w=1.0)
#' @param genetic_code Named vector mapping codons to amino acids
#' @param low_expr_quantile Quantile for estimating mutation bias (default: 0.10)
#' @return Data frame with Gene_name, P_preferred, S, and expression data
calculate_selection_coefficients <- function(codon_usage, expression_data, 
                                            preferred_codons, genetic_code,
                                            low_expr_quantile = 0.10) {
  
  cat("\n=== Calculating Selection Coefficients (S = 4Nes) ===\n")
  
  # Step 1: Estimate mutation bias from low-expression genes
  mutation_params <- estimate_mutation_bias(codon_usage, expression_data, 
                                            preferred_codons, genetic_code,
                                            low_expr_quantile)
  
  M <- mutation_params$mutation_bias
  
  # Step 2: Calculate P(preferred) for all genes
  cat("\nCalculating preferred codon proportion for all genes...\n")
  
  P_preferred_values <- apply(codon_usage, 1, function(gene_row) {
    calculate_preferred_proportion(gene_row, preferred_codons, genetic_code)
  })
  
  # Step 3: Calculate S for all genes
  cat("Calculating S = 4Nes for all genes...\n")
  
  S_values <- sapply(P_preferred_values, function(P) {
    calculate_S(P, M)
  })
  
  # Combine results
  results <- data.frame(
    Gene_name = codon_usage$Gene_name,
    P_preferred = P_preferred_values,
    S = S_values,
    stringsAsFactors = FALSE
  )
  
  # Merge with expression data
  results <- results |>
    left_join(expression_data, by = "Gene_name")
  
  # Summary statistics
  cat("\n=== S (4Nes) Summary Statistics ===\n")
  cat(sprintf("Mean S: %.4f\n", mean(results$S, na.rm = TRUE)))
  cat(sprintf("Median S: %.4f\n", median(results$S, na.rm = TRUE)))
  cat(sprintf("SD S: %.4f\n", sd(results$S, na.rm = TRUE)))
  cat(sprintf("Range: %.4f to %.4f\n", 
              min(results$S, na.rm = TRUE), 
              max(results$S, na.rm = TRUE)))
  
  # Count genes with positive vs negative selection
  n_positive <- sum(results$S > 0, na.rm = TRUE)
  n_negative <- sum(results$S < 0, na.rm = TRUE)
  n_total <- sum(!is.na(results$S))
  
  cat(sprintf("\nGenes with positive selection (S > 0): %d / %d (%.1f%%)\n",
              n_positive, n_total, 100 * n_positive / n_total))
  cat(sprintf("Genes with negative selection (S < 0): %d / %d (%.1f%%)\n",
              n_negative, n_total, 100 * n_negative / n_total))
  
  # Correlation with expression (with significance test)
  if ("Expression" %in% names(results)) {
    # Filter complete cases
    complete_data <- results[!is.na(results$S) & !is.na(results$Expression), ]
    
    # Spearman correlation test
    cor_test <- cor.test(complete_data$S, log2(complete_data$Expression + 1), 
                        method = "spearman", exact = FALSE)
    
    cat(sprintf("\nSpearman correlation: S vs log2(Expression+1) = %.4f\n", cor_test$estimate))
    cat(sprintf("  p-value: %.2e\n", cor_test$p.value))
    cat(sprintf("  n = %d genes\n", nrow(complete_data)))
    
    if (cor_test$p.value < 0.001) {
      cat("  *** Highly significant (p < 0.001)\n")
    } else if (cor_test$p.value < 0.01) {
      cat("  ** Very significant (p < 0.01)\n")
    } else if (cor_test$p.value < 0.05) {
      cat("  * Significant (p < 0.05)\n")
    } else {
      cat("  Not significant (p >= 0.05)\n")
    }
    
    # Also test Pearson for comparison
    cor_test_pearson <- cor.test(complete_data$S, log2(complete_data$Expression + 1), 
                                 method = "pearson")
    cat(sprintf("\nPearson correlation: r = %.4f, p = %.2e\n", 
                cor_test_pearson$estimate, cor_test_pearson$p.value))
  }
  
  # Store mutation bias parameters
  attr(results, "mutation_bias") <- M
  attr(results, "P_neutral") <- mutation_params$P_neutral
  
  return(results)
}


#' Calculate individual selection coefficient (s) given Ne
#'
#' @param S Population-scaled selection coefficient (S = 4Nes)
#' @param Ne Effective population size
#' @return Individual selection coefficient s
calculate_s_from_S <- function(S, Ne) {
  return(S / (4 * Ne))
}


#' Analyze selection coefficients across different Ne values
#'
#' @param selection_results Data frame from calculate_selection_coefficients()
#' @param Ne_values Vector of Ne values to test (default: literature range)
#' @return Data frame with s estimates for different Ne values
analyze_Ne_sensitivity <- function(selection_results, 
                                   Ne_values = c(1e5, 2e5, 3e5, 5e5, 1e6)) {
  
  cat("\n=== Sensitivity Analysis: s across different Ne values ===\n")
  cat("Testing Ne values (from literature):\n")
  
  results_list <- lapply(Ne_values, function(Ne) {
    s_values <- calculate_s_from_S(selection_results$S, Ne)
    
    data.frame(
      Ne = Ne,
      mean_s = mean(s_values, na.rm = TRUE),
      median_s = median(s_values, na.rm = TRUE),
      sd_s = sd(s_values, na.rm = TRUE),
      min_s = min(s_values, na.rm = TRUE),
      max_s = max(s_values, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  
  sensitivity_df <- do.call(rbind, results_list)
  
  # Print table
  cat("\n")
  print(sensitivity_df, digits = 6)
  
  # Check if s ≈ 1/Ne (the weak selection regime)
  cat("\n=== Testing Weak Selection Regime (s ≈ 1/Ne) ===\n")
  sensitivity_df$s_times_Ne <- sensitivity_df$median_s * sensitivity_df$Ne
  
  cat("If selection is weak, s·Ne should be ~ 1\n")
  cat("Observed s·Ne values:\n")
  print(sensitivity_df[, c("Ne", "median_s", "s_times_Ne")], digits = 4)
  
  return(sensitivity_df)
}


#' Create visualization of S vs expression
#'
#' @param selection_results Data frame from calculate_selection_coefficients()
#' @param output_file Path to save plot (default: "./results/S_vs_expression.pdf")
plot_S_vs_expression <- function(selection_results, 
                                 output_file = "./results/S_vs_expression.pdf") {
  
  library(ggplot2)
  
  # Filter out NA values
  plot_data <- selection_results |>
    filter(!is.na(S), !is.na(Expression), Expression > 0)
  
  # Calculate correlation
  cor_val <- cor(plot_data$S, log2(plot_data$Expression + 1), 
                 method = "spearman")
  
  p <- ggplot(plot_data, aes(x = log2(Expression + 1), y = S)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "loess", color = "red", se = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = "Selection Strength vs Gene Expression",
      subtitle = sprintf("Spearman ρ = %.3f", cor_val),
      x = "log2(Expression + 1)",
      y = "S = 4Nes (Population-scaled selection coefficient)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  ggsave(output_file, p, width = 8, height = 6)
  cat(sprintf("\nPlot saved to: %s\n", output_file))
  
  return(p)
}


#' Create histogram of S values
#'
#' @param selection_results Data frame from calculate_selection_coefficients()
#' @param output_file Path to save plot
plot_S_distribution <- function(selection_results,
                               output_file = "./results/S_distribution.pdf") {
  
  library(ggplot2)
  
  plot_data <- selection_results |> filter(!is.na(S))
  
  p <- ggplot(plot_data, aes(x = S)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
    geom_vline(xintercept = median(plot_data$S, na.rm = TRUE), 
               linetype = "dashed", color = "darkgreen", size = 1) +
    labs(
      title = "Distribution of Selection Coefficients",
      subtitle = sprintf("Median S = %.3f", median(plot_data$S, na.rm = TRUE)),
      x = "S = 4Nes (Population-scaled selection coefficient)",
      y = "Number of genes"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12)
    )
  
  ggsave(output_file, p, width = 8, height = 6)
  cat(sprintf("\nPlot saved to: %s\n", output_file))
  
  return(p)
}
