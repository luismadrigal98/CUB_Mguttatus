count_preferred_codons <- function(gene_list, codon_usage_data, preferred_list) {
  
  # Get codon usage for these genes - convert to data.frame to avoid data.table issues
  gene_usage <- codon_usage_data |>
    dplyr::filter(Gene_name %in% gene_list) |>
    as.data.frame()
  
  # For each gene, calculate proportion of preferred codons used
  preferred_codons_vec <- preferred_list$Codon
  
  results <- data.frame(
    Gene_name = gene_usage$Gene_name,
    Total_Codons = 0,
    Preferred_Codons = 0,
    Preferred_Proportion = 0
  )
  
  codon_cols <- setdiff(names(gene_usage), "Gene_name")
  
  for (i in 1:nrow(gene_usage)) {
    gene <- gene_usage$Gene_name[i]
    
    # Total codons (all synonymous codons)
    total <- sum(as.numeric(gene_usage[i, codon_cols]), na.rm = TRUE)
    
    # Preferred codons
    preferred_cols <- intersect(preferred_codons_vec, codon_cols)
    preferred <- sum(as.numeric(gene_usage[i, preferred_cols]), na.rm = TRUE)
    
    results$Total_Codons[i] <- total
    results$Preferred_Codons[i] <- preferred
    results$Preferred_Proportion[i] <- if(total > 0) preferred / total else NA
  }
  
  return(results)
}