#!/usr/bin/env Rscript
# calculate_cai.R
# Calculate Codon Adaptation Index (CAI) for genes
# CAI measures the degree of bias towards codons preferred in highly expressed genes

#' Calculate relative adaptiveness (w) for codons from reference set
#'
#' @param reference_codon_counts Data frame with codon counts from highly expressed genes
#'                               Columns: Gene_name, then one column per codon
#' @param genetic_code Data frame mapping codons to amino acids
#' @return Data frame with codon, amino_acid, frequency, and relative_adaptiveness (w)
calculate_relative_adaptiveness <- function(reference_codon_counts, genetic_code) {
  
  # Sum codon counts across all reference genes
  codon_sums <- colSums(reference_codon_counts[, -1], na.rm = TRUE)
  
  # Create data frame with codon frequencies
  codon_freq <- data.frame(
    codon = names(codon_sums),
    count = as.numeric(codon_sums),
    stringsAsFactors = FALSE
  )
  
  # Add amino acid information
  codon_freq <- codon_freq |>
    left_join(genetic_code, by = "codon")
  
  # For each amino acid, find the maximum frequency (most preferred codon)
  # and calculate relative adaptiveness w = freq_codon / freq_max
  codon_freq <- codon_freq |>
    group_by(amino_acid) |>
    mutate(
      max_count = max(count, na.rm = TRUE),
      relative_adaptiveness = count / max_count
    ) |>
    ungroup()
  
  # Handle stop codons (set w = 1 since they don't contribute to CAI)
  codon_freq$relative_adaptiveness[codon_freq$amino_acid == "STOP"] <- 1.0
  
  return(codon_freq)
}


#' Calculate CAI for a single gene
#'
#' @param gene_codon_counts Named vector of codon counts for one gene
#' @param w_values Named vector of relative adaptiveness values (names = codons)
#' @return CAI value (0 to 1)
calculate_gene_cai <- function(gene_codon_counts, w_values) {
  
  # Remove gene name if present
  if ("Gene_name" %in% names(gene_codon_counts)) {
    gene_codon_counts <- gene_codon_counts[names(gene_codon_counts) != "Gene_name"]
  }
  
  # Get codons that are actually used in this gene (count > 0)
  used_codons <- names(gene_codon_counts)[gene_codon_counts > 0]
  
  if (length(used_codons) == 0) {
    return(NA)  # No codons found
  }
  
  # Get w values for used codons
  w_used <- w_values[used_codons]
  
  # Get counts for used codons
  counts_used <- gene_codon_counts[used_codons]
  
  # CAI = geometric mean of w values, weighted by codon counts
  # CAI = exp( (1/L) * sum(n_i * ln(w_i)) )
  # where L = total number of codons, n_i = count of codon i, w_i = adaptiveness of codon i
  
  total_codons <- sum(counts_used)
  
  if (total_codons == 0) {
    return(NA)
  }
  
  # Calculate weighted log sum
  # Using log(w) and handling potential zeros
  log_w <- log(w_used)
  log_w[is.infinite(log_w)] <- 0  # Handle w=0 cases (shouldn't happen but be safe)
  
  weighted_log_sum <- sum(counts_used * log_w)
  
  # CAI = exp(weighted_log_sum / total_codons)
  cai <- exp(weighted_log_sum / total_codons)
  
  return(cai)
}


#' Calculate CAI for all genes
#'
#' @param codon_counts Data frame with codon counts for all genes (first column = Gene_name)
#' @param reference_genes Character vector of gene names for reference set (highly expressed)
#' @param genetic_code Data frame mapping codons to amino acids
#' @return Data frame with Gene_name and CAI
calculate_cai <- function(codon_counts, reference_genes, genetic_code) {
  
  cat("\n=== Calculating Codon Adaptation Index (CAI) ===\n")
  
  # 1. Extract reference set (highly expressed genes)
  reference_set <- codon_counts |>
    filter(Gene_name %in% reference_genes)
  
  cat(sprintf("Reference set: %d highly expressed genes\n", nrow(reference_set)))
  
  # 2. Calculate relative adaptiveness (w) for each codon
  cat("Calculating relative adaptiveness (w) for each codon...\n")
  w_table <- calculate_relative_adaptiveness(reference_set, genetic_code)
  
  # Create named vector for easy lookup
  w_values <- setNames(w_table$relative_adaptiveness, w_table$codon)
  
  # Print some statistics about optimal codons
  cat("\nOptimal codons (w = 1.0) in reference set:\n")
  optimal <- w_table |>
    filter(relative_adaptiveness == 1.0, amino_acid != "STOP") |>
    select(amino_acid, codon, count) |>
    arrange(amino_acid)
  print(optimal)
  
  # 3. Calculate CAI for all genes
  cat(sprintf("\nCalculating CAI for %d genes...\n", nrow(codon_counts)))
  
  cai_values <- apply(codon_counts, 1, function(gene_row) {
    calculate_gene_cai(gene_row, w_values)
  })
  
  # Create result data frame
  result <- data.frame(
    Gene_name = codon_counts$Gene_name,
    CAI = cai_values,
    stringsAsFactors = FALSE
  )
  
  # Summary statistics
  cat("\n=== CAI Summary Statistics ===\n")
  cat(sprintf("Mean CAI: %.4f\n", mean(result$CAI, na.rm = TRUE)))
  cat(sprintf("Median CAI: %.4f\n", median(result$CAI, na.rm = TRUE)))
  cat(sprintf("SD CAI: %.4f\n", sd(result$CAI, na.rm = TRUE)))
  cat(sprintf("Range: %.4f - %.4f\n", min(result$CAI, na.rm = TRUE), max(result$CAI, na.rm = TRUE)))
  
  # Check correlation between reference genes and CAI
  reference_cai <- result |> filter(Gene_name %in% reference_genes)
  cat(sprintf("\nMean CAI in reference set: %.4f\n", mean(reference_cai$CAI, na.rm = TRUE)))
  
  return(list(
    cai_values = result,
    w_table = w_table,
    reference_genes = reference_genes
  ))
}
