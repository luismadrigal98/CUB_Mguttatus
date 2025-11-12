function(codon_usage, genetic_code, expression_data = NULL, n_bootstrap = 100, n_cores = NULL) {
  
  cat("\n=== Integrating CDC Analysis with Main Pipeline ===\n")
  
  # Remove .1 suffix from gene names if present (to match main analysis)
  if (any(grepl("\\.1$", codon_usage$Gene_name))) {
    cat("Removing .1 suffix from gene names to match main analysis format\n")
    codon_usage$Gene_name <- sub("\\.1$", "", codon_usage$Gene_name)
  }
  
  # Calculate CDC for all genes
  cdc_results <- calculate_cdc_all(codon_usage, genetic_code, n_bootstrap, n_cores)
  
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