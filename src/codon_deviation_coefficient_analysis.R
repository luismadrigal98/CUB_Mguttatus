#' @description Collection of functions to implement the calculation and statistical testing
#' of CUB based on the codon_deviation_coefficient, described in detail in:
#' 
#' @reference Zhang et al. BMC Bioinformatics 2012, 13:43, 
#' http://www.biomedcentral.com/1471-2105/13/43
#' 
#' @note
#' This metric takes into account the GC content and the purine contents as
#' background nucleotide composition from observed positional GC and purine contents.
#' Analyzes 59 codons (excluding stop codons, Met, and Trp which have no synonymous choice).
#' 
#' @author (script): Luis J. Madrigal-Roca
#' @date 10/31/2025
#' _____________________________________________________________________________

# Codon Deviation Coefficient (CDC) Implementation in R
# Based on Zhang et al. (2012) BMC Bioinformatics

# Load required library
if (!require("Biostrings")) {
  if (!require("BiocManager")) install.packages("BiocManager")
  BiocManager::install("Biostrings")
}
library(Biostrings)

# Define genetic code (59 sense codons excluding single-codon amino acids) - compatible with main analysis
get_sense_codons <- function(genetic_code) {
  # Extract codons that are not stop codons
  sense_codons <- names(genetic_code)[genetic_code != "STOP"]
  
  # Remove single-codon amino acids (no synonymous codon choice)
  # Methionine (ATG) and Tryptophan (TGG) have only one codon each
  single_codon_aa <- c("ATG", "TGG")  # Met and Trp
  sense_codons <- sense_codons[!sense_codons %in% single_codon_aa]
  
  return(sense_codons)
}

#' Calculate nucleotide composition at each codon position from codon counts
#'
#' @param codon_counts Named vector of codon counts for one gene
#' @param genetic_code Named vector mapping codons to amino acids
#' @return List with GC and purine contents at each codon position
get_positional_composition_from_counts <- function(codon_counts, genetic_code) {
  
  # Get sense codons (exclude stop codons)
  sense_codons <- get_sense_codons(genetic_code)
  
  # Filter to sense codons only
  codon_counts <- codon_counts[names(codon_counts) %in% sense_codons]
  codon_counts <- codon_counts[!is.na(codon_counts)]
  
  if (length(codon_counts) == 0 || sum(codon_counts) == 0) {
    warning("No valid codon counts found")
    return(list(S1 = 0.5, S2 = 0.5, S3 = 0.5, R1 = 0.5, R2 = 0.5, R3 = 0.5))
  }
  
  # Initialize nucleotide counts for each position
  pos1_counts <- c(Ai = 0, Ti = 0, Gi = 0, Ci = 0)
  pos2_counts <- c(Ai = 0, Ti = 0, Gi = 0, Ci = 0)
  pos3_counts <- c(Ai = 0, Ti = 0, Gi = 0, Ci = 0)
  
  # Count nucleotides at each position weighted by codon frequency
  for (codon in names(codon_counts)) {
    if (nchar(codon) == 3) {
      count <- codon_counts[codon]
      
      # Skip if count is NA or 0
      if (is.na(count) || count == 0) next
      
      bases <- strsplit(codon, "")[[1]]
      
      # Map base letters to variable names to avoid R conflicts
      base1 <- paste0(bases[1], "i")
      base2 <- paste0(bases[2], "i")
      base3 <- paste0(bases[3], "i")
      
      # Add count to each position, with NA protection
      pos1_counts[base1] <- pos1_counts[base1] + count
      pos2_counts[base2] <- pos2_counts[base2] + count
      pos3_counts[base3] <- pos3_counts[base3] + count
    }
  }
  
  # Calculate GC content (S) for each position
  calc_gc <- function(counts) {
    total <- sum(counts)
    if (total == 0) return(0.5)
    gc <- (counts["Gi"] + counts["Ci"]) / total
    # Check for NA
    if (is.na(gc)) {
      warning(sprintf("NA in GC calculation: Gi=%s, Ci=%s, total=%s", 
                     counts["Gi"], counts["Ci"], total))
      return(0.5)
    }
    return(gc)
  }
  
  # Calculate purine content (R) for each position
  calc_purine <- function(counts) {
    total <- sum(counts)
    if (total == 0) return(0.5)
    purine <- (counts["Ai"] + counts["Gi"]) / total
    # Check for NA
    if (is.na(purine)) {
      warning(sprintf("NA in purine calculation: Ai=%s, Gi=%s, total=%s", 
                     counts["Ai"], counts["Gi"], total))
      return(0.5)
    }
    return(purine)
  }
  
  list(
    S1 = calc_gc(pos1_counts), S2 = calc_gc(pos2_counts), S3 = calc_gc(pos3_counts),
    R1 = calc_purine(pos1_counts), R2 = calc_purine(pos2_counts), R3 = calc_purine(pos3_counts)
  )
}

#' Calculate expected nucleotide frequencies at each position
#'
#' @param S GC content at position i
#' @param R Purine content at position i
#' @return Named vector with expected frequencies of A, T, G, C
calc_expected_nucleotides <- function(S, R) {
  c(
    Ai = (1 - S) * R,
    Ti = (1 - S) * (1 - R),
    Gi = S * R,
    Ci = S * (1 - R)
  )
}

#' Calculate expected codon usage based on positional composition
#'
#' @param comp List of positional compositions (S1, S2, S3, R1, R2, R3)
#' @param genetic_code Named vector mapping codons to amino acids
#' @return Named vector of expected codon frequencies
calc_expected_codon_usage <- function(comp, genetic_code) {
  # Get sense codons
  sense_codons <- get_sense_codons(genetic_code)
  
  # Get expected nucleotide frequencies for each position
  nuc1 <- calc_expected_nucleotides(comp$S1, comp$R1)
  nuc2 <- calc_expected_nucleotides(comp$S2, comp$R2)
  nuc3 <- calc_expected_nucleotides(comp$S3, comp$R3)
  
  # Calculate expected usage for all sense codons
  expected <- numeric(length(sense_codons))
  names(expected) <- sense_codons
  
  for (codon in sense_codons) {
    bases <- strsplit(codon, "")[[1]]
    base1 <- paste0(bases[1], "i")
    base2 <- paste0(bases[2], "i")
    base3 <- paste0(bases[3], "i")
    expected[codon] <- nuc1[base1] * nuc2[base2] * nuc3[base3]
  }
  
  # Normalize
  expected / sum(expected)
}

#' Calculate observed codon usage from codon counts
#'
#' @param codon_counts Named vector of codon counts for one gene
#' @param genetic_code Named vector mapping codons to amino acids
#' @return Named vector of observed codon frequencies
calc_observed_codon_usage <- function(codon_counts, genetic_code) {
  # Get sense codons
  sense_codons <- get_sense_codons(genetic_code)
  
  # Filter to sense codons and remove NAs
  observed <- codon_counts[names(codon_counts) %in% sense_codons]
  observed <- observed[!is.na(observed)]
  
  # Check if we have any counts
  if (length(observed) == 0 || sum(observed) == 0) {
    warning("No valid codon counts found")
    # Return uniform distribution
    uniform <- rep(1/length(sense_codons), length(sense_codons))
    names(uniform) <- sense_codons
    return(uniform)
  }
  
  # Ensure all sense codons are represented (fill missing with 0)
  full_observed <- numeric(length(sense_codons))
  names(full_observed) <- sense_codons
  full_observed[names(observed)] <- observed
  
  # Normalize to frequencies
  total <- sum(full_observed)
  if (total == 0) {
    # Return uniform distribution if no counts
    return(rep(1/length(sense_codons), length(sense_codons)))
  }
  
  full_observed / total
}

#' Calculate CDC using cosine distance metric
#'
#' @param expected Named vector of expected codon frequencies
#' @param observed Named vector of observed codon frequencies
#' @return CDC value (0 = no bias, 1 = maximum bias)
calc_cdc <- function(expected, observed) {
  # Ensure vectors are aligned
  codons <- intersect(names(expected), names(observed))
  expected <- expected[codons]
  observed <- observed[codons]
  
  # Calculate cosine similarity
  numerator <- sum(expected * observed)
  denominator <- sqrt(sum(expected^2)) * sqrt(sum(observed^2))
  cosine_sim <- numerator / denominator
  
  # CDC is 1 - cosine similarity
  1 - cosine_sim
}

#' Generate random codon counts based on positional composition
#'
#' @param comp List of positional compositions
#' @param total_codons Total number of codons to generate
#' @param genetic_code Named vector mapping codons to amino acids
#' @return Named vector of random codon counts
generate_random_codon_counts <- function(comp, total_codons, genetic_code) {
  nuc1 <- calc_expected_nucleotides(comp$S1, comp$R1)
  nuc2 <- calc_expected_nucleotides(comp$S2, comp$R2)
  nuc3 <- calc_expected_nucleotides(comp$S3, comp$R3)
  
  # Check for NA values in nucleotide frequencies
  if (any(is.na(c(nuc1, nuc2, nuc3)))) {
    warning("NA values in nucleotide frequencies, using uniform distribution")
    sense_codons <- get_sense_codons(genetic_code)
    uniform_counts <- rep(total_codons %/% length(sense_codons), length(sense_codons))
    names(uniform_counts) <- sense_codons
    return(uniform_counts)
  }
  
  # Get sense codons
  sense_codons <- get_sense_codons(genetic_code)
  
  # Calculate expected frequencies for each codon
  expected_freqs <- numeric(length(sense_codons))
  names(expected_freqs) <- sense_codons
  
  for (codon in sense_codons) {
    bases <- strsplit(codon, "")[[1]]
    base1 <- paste0(bases[1], "i")
    base2 <- paste0(bases[2], "i")
    base3 <- paste0(bases[3], "i")
    expected_freqs[codon] <- nuc1[base1] * nuc2[base2] * nuc3[base3]
  }
  
  # Check for NA or zero sum in expected frequencies
  if (any(is.na(expected_freqs)) || sum(expected_freqs) == 0) {
    warning("Invalid expected frequencies, using uniform distribution")
    uniform_counts <- rep(total_codons %/% length(sense_codons), length(sense_codons))
    names(uniform_counts) <- sense_codons
    return(uniform_counts)
  }
  
  # Normalize frequencies
  expected_freqs <- expected_freqs / sum(expected_freqs)
  
  # Generate random counts using multinomial distribution
  random_counts <- rmultinom(1, total_codons, expected_freqs)[, 1]
  names(random_counts) <- sense_codons
  
  return(random_counts)
}

#' Calculate CDC with bootstrap significance test for single gene
#'
#' @param codon_counts Named vector of codon counts for one gene
#' @param genetic_code Named vector mapping codons to amino acids
#' @param n_bootstrap Number of bootstrap replicates (default: 1000)
#' @return List containing CDC value, p-value, and bootstrap distribution
calculate_cdc_single <- function(codon_counts, genetic_code, n_bootstrap = 1000) {
  
  # Handle data.frame/data.table row input
  if (is.data.frame(codon_counts)) {
    # Convert single row data.frame to named vector
    # Get column names excluding Gene_name
    codon_cols <- names(codon_counts)[names(codon_counts) != "Gene_name"]
    
    # Extract values as numeric vector - universal approach for data.table and data.frame
    # Convert to standard data.frame first to avoid data.table syntax issues
    codon_counts_df <- as.data.frame(codon_counts)
    codon_values <- as.numeric(unlist(codon_counts_df[1, codon_cols]))
    names(codon_values) <- codon_cols
    codon_counts <- codon_values
    
  } else {
    # Handle named vector input
    if ("Gene_name" %in% names(codon_counts)) {
      codon_counts <- codon_counts[names(codon_counts) != "Gene_name"]
    }
    
    # Convert to numeric while preserving names
    count_names <- names(codon_counts)
    codon_counts <- as.numeric(codon_counts)
    names(codon_counts) <- count_names
  }
  
  # Get total number of codons
  total_codons <- sum(codon_counts, na.rm = TRUE)
  
  if (total_codons == 0) {
    warning("No codon counts found")
    return(list(CDC = NA, p_value = NA, bootstrap_distribution = NULL))
  }
  
  # Calculate positional composition from codon counts
  comp <- get_positional_composition_from_counts(codon_counts, genetic_code)
  
  # Check for invalid composition values
  comp_values <- unlist(comp)
  if (any(is.na(comp_values)) || any(comp_values < 0) || any(comp_values > 1)) {
    warning("Invalid positional composition values")
    return(list(CDC = NA, p_value = NA, bootstrap_distribution = NULL))
  }
  
  # Calculate expected and observed codon usage
  expected <- calc_expected_codon_usage(comp, genetic_code)
  observed <- calc_observed_codon_usage(codon_counts, genetic_code)
  
  # Check for NA values in expected or observed
  if (any(is.na(expected)) || any(is.na(observed))) {
    warning("NA values in expected or observed codon usage")
    return(list(CDC = NA, p_value = NA, bootstrap_distribution = NULL))
  }
  
  # Calculate CDC
  cdc_value <- calc_cdc(expected, observed)
  
  # Bootstrap resampling
  bootstrap_cdc <- numeric(n_bootstrap)
  
  for (i in 1:n_bootstrap) {
    # Generate random codon counts with same composition and total
    random_counts <- generate_random_codon_counts(comp, total_codons, genetic_code)
    
    # Calculate CDC for random counts
    random_observed <- calc_observed_codon_usage(random_counts, genetic_code)
    bootstrap_cdc[i] <- calc_cdc(expected, random_observed)
  }
  
  # Calculate two-sided p-value
  p_lower <- sum(bootstrap_cdc <= cdc_value) / n_bootstrap
  p_upper <- sum(bootstrap_cdc >= cdc_value) / n_bootstrap
  p_value <- 2 * min(p_lower, p_upper)
  p_value <- min(p_value, 1)
  
  list(
    CDC = cdc_value,
    p_value = p_value,
    bootstrap_distribution = bootstrap_cdc,
    expected_usage = expected,
    observed_usage = observed,
    composition = comp,
    total_codons = total_codons
  )
}

#' Calculate CDC for all genes in codon usage data frame
#'
#' @param codon_usage_df Data frame with Gene_name column and codon count columns
#' @param genetic_code Named vector mapping codons to amino acids
#' @param n_bootstrap Number of bootstrap replicates per gene (default: 1000)
#' @return Data frame with Gene_name, CDC, and p_value columns
calculate_cdc_all <- function(codon_usage_df, genetic_code, n_bootstrap = 1000) {
  
  cat(sprintf("\n=== Calculating CDC for %d genes ===\n", nrow(codon_usage_df)))
  cat(sprintf("Bootstrap replicates per gene: %d\n", n_bootstrap))
  cat(sprintf("Data structure: %s\n", class(codon_usage_df)[1]))
  cat(sprintf("Columns: %s\n", paste(head(names(codon_usage_df), 10), collapse = ", ")))
  
  # Initialize results
  results <- data.frame(
    Gene_name = codon_usage_df$Gene_name,
    CDC = numeric(nrow(codon_usage_df)),
    p_value = numeric(nrow(codon_usage_df)),
    stringsAsFactors = FALSE
  )
  
  # Calculate CDC for each gene
  for (i in seq_len(nrow(codon_usage_df))) {
    if (i %% 1000 == 0) {
      cat(sprintf("  Processed %d / %d genes (%.1f%%)\n", 
                  i, nrow(codon_usage_df), 100 * i / nrow(codon_usage_df)))
    }
    
    # Extract codon counts for this gene
    # Convert single row to named vector to avoid data.table issues
    gene_counts <- codon_usage_df[i, , drop = FALSE]
    
    # Calculate CDC
    cdc_result <- calculate_cdc_single(gene_counts, genetic_code, n_bootstrap)
    
    # Store results
    results$CDC[i] <- cdc_result$CDC
    results$p_value[i] <- cdc_result$p_value
  }
  
  cat("CDC calculation complete!\n")
  
  # Summary statistics
  cat("\n=== CDC Summary Statistics ===\n")
  cat(sprintf("Mean CDC: %.4f\n", mean(results$CDC, na.rm = TRUE)))
  cat(sprintf("Median CDC: %.4f\n", median(results$CDC, na.rm = TRUE)))
  cat(sprintf("SD CDC: %.4f\n", sd(results$CDC, na.rm = TRUE)))
  cat(sprintf("Range: %.4f - %.4f\n", min(results$CDC, na.rm = TRUE), max(results$CDC, na.rm = TRUE)))
  
  # Significance summary
  sig_count <- sum(results$p_value < 0.05, na.rm = TRUE)
  cat(sprintf("\nSignificant CDC values (p < 0.05): %d / %d (%.1f%%)\n", 
              sig_count, sum(!is.na(results$p_value)), 
              100 * sig_count / sum(!is.na(results$p_value))))
  
  return(results)
}

#' Plot bootstrap distribution with observed CDC
#'
#' @param cdc_result Result object from calculate_cdc_single()
plot_cdc_bootstrap <- function(cdc_result) {
  if (is.null(cdc_result$bootstrap_distribution)) {
    warning("No bootstrap distribution available")
    return(NULL)
  }
  
  hist(cdc_result$bootstrap_distribution, 
       breaks = 50,
       main = "Bootstrap Distribution of CDC",
       xlab = "CDC Value",
       ylab = "Frequency",
       col = "lightblue",
       border = "white")
  abline(v = cdc_result$CDC, col = "red", lwd = 2, lty = 2)
  legend("topright", 
         legend = c(
           paste("Observed CDC =", round(cdc_result$CDC, 4)),
           paste("P-value =", round(cdc_result$p_value, 4))
         ),
         col = c("red", NA),
         lty = c(2, NA),
         lwd = c(2, NA),
         bty = "n")
}

#' Integrate CDC analysis with main analysis pipeline
#'
#' @param codon_usage Data frame from main analysis (with Gene_name column)
#' @param genetic_code Named vector from main analysis
#' @param expression_data Optional: data frame with Gene_name and expression info
#' @param n_bootstrap Number of bootstrap replicates (default: 100 for speed)
#' @return Data frame with CDC results, optionally merged with expression data
integrate_cdc_analysis <- function(codon_usage, genetic_code, expression_data = NULL, n_bootstrap = 100) {
  
  cat("\n=== Integrating CDC Analysis with Main Pipeline ===\n")
  
  # Remove .1 suffix from gene names if present (to match main analysis)
  if (any(grepl("\\.1$", codon_usage$Gene_name))) {
    cat("Removing .1 suffix from gene names to match main analysis format\n")
    codon_usage$Gene_name <- sub("\\.1$", "", codon_usage$Gene_name)
  }
  
  # Calculate CDC for all genes
  cdc_results <- calculate_cdc_all(codon_usage, genetic_code, n_bootstrap)
  
  # Merge with expression data if provided
  if (!is.null(expression_data)) {
    cat("Merging CDC results with expression data\n")
    final_results <- merge(expression_data, cdc_results, by = "Gene_name", all.x = TRUE)
    
    # Check merge success
    merged_count <- sum(!is.na(final_results$CDC))
    cat(sprintf("Successfully merged %d / %d genes with expression data\n", 
                merged_count, nrow(expression_data)))
    
    return(final_results)
  } else {
    return(cdc_results)
  }
}

#' Quick CDC analysis for testing (smaller bootstrap)
#'
#' @param codon_usage Data frame with codon counts
#' @param genetic_code Named vector mapping codons to amino acids
#' @param n_genes Number of genes to analyze (for testing)
#' @return Data frame with CDC results for subset of genes
quick_cdc_test <- function(codon_usage, genetic_code, n_genes = 100) {
  
  cat(sprintf("\n=== Quick CDC Test (%d genes, 50 bootstrap replicates each) ===\n", n_genes))
  
  # Take random subset
  if (nrow(codon_usage) > n_genes) {
    test_genes <- sample(seq_len(nrow(codon_usage)), n_genes)
    codon_subset <- codon_usage[test_genes, ]
  } else {
    codon_subset <- codon_usage
  }
  
  # Calculate CDC with reduced bootstrap
  cdc_results <- calculate_cdc_all(codon_subset, genetic_code, n_bootstrap = 50)
  
  return(cdc_results)
}