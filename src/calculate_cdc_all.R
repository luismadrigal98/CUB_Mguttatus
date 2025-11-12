calculate_cdc_all <- function(codon_usage_df, genetic_code, n_bootstrap = 1000, n_cores = NULL) {
  
  cat(sprintf("\n=== Calculating CDC for %d genes ===\n", nrow(codon_usage_df)))
  cat(sprintf("Bootstrap replicates per gene: %d\n", n_bootstrap))
  cat(sprintf("Data structure: %s\n", class(codon_usage_df)[1]))
  cat(sprintf("Columns: %s\n", paste(head(names(codon_usage_df), 10), collapse = ", ")))
  
  # Setup parallel processing
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 1
  }
  n_cores <- max(1, min(n_cores, parallel::detectCores()))
  
  cat(sprintf("Using %d cores for parallel processing\n", n_cores))
  
  # Convert to standard data.frame to avoid data.table issues in parallel
  codon_usage_df <- as.data.frame(codon_usage_df)
  
  # Setup parallel backend
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl))
  
  # Export necessary objects and data to cluster
  parallel::clusterExport(cl, c("calculate_cdc_single", "get_positional_composition_from_counts",
                                "calc_expected_codon_usage", "calc_observed_codon_usage",
                                "calc_cdc", "generate_random_codon_counts",
                                "calc_expected_nucleotides", "get_sense_codons",
                                "genetic_code"),
                          envir = environment())
  
  # Calculate CDC for each gene in parallel with progress
  cat("Processing genes...\n")
  start_time <- Sys.time()
  
  results_list <- parallel::parLapplyLB(cl, seq_len(nrow(codon_usage_df)), function(i) {
    # Progress reporting handled by main thread, not workers
    
    # Extract codon counts for this gene
    gene_counts <- codon_usage_df[i, , drop = FALSE]
    
    # Calculate CDC
    cdc_result <- calculate_cdc_single(gene_counts, genetic_code, n_bootstrap)
    
    # Return just the essentials
    list(CDC = cdc_result$CDC, p_value = cdc_result$p_value)
  })
  
  # Combine results
  results <- data.frame(
    Gene_name = codon_usage_df$Gene_name,
    CDC = sapply(results_list, function(x) x$CDC),
    p_value = sapply(results_list, function(x) x$p_value),
    stringsAsFactors = FALSE
  )
  
  # Apply FDR correction for multiple testing
  results$p_adj <- p.adjust(results$p_value, method = "BH")
  
  elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(sprintf("CDC calculation complete in %.1f seconds (%.2f genes/sec)!\n", 
              elapsed_time, nrow(codon_usage_df) / elapsed_time))
  
  # Summary statistics
  cat("\n=== CDC Summary Statistics ===\n")
  cat(sprintf("Mean CDC: %.4f\n", mean(results$CDC, na.rm = TRUE)))
  cat(sprintf("Median CDC: %.4f\n", median(results$CDC, na.rm = TRUE)))
  cat(sprintf("SD CDC: %.4f\n", sd(results$CDC, na.rm = TRUE)))
  cat(sprintf("Range: %.4f - %.4f\n", min(results$CDC, na.rm = TRUE), max(results$CDC, na.rm = TRUE)))
  
  # Significance summary - uncorrected and FDR-corrected
  sig_count_raw <- sum(results$p_value < 0.05, na.rm = TRUE)
  sig_count_fdr <- sum(results$p_adj < 0.05, na.rm = TRUE)
  n_valid <- sum(!is.na(results$p_value))
  
  cat(sprintf("\nSignificant CDC values (uncorrected p < 0.05): %d / %d (%.1f%%)\n", 
              sig_count_raw, n_valid, 100 * sig_count_raw / n_valid))
  cat(sprintf("Significant CDC values (FDR-adjusted q < 0.05): %d / %d (%.1f%%)\n", 
              sig_count_fdr, n_valid, 100 * sig_count_fdr / n_valid))
  
  return(results)
}