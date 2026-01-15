count_preferred_by_aa <- function(codon_data, preferred_codons, genetic_code) {
  
  # Convert to data.frame if it's a data.table (avoid data.table indexing issues)
  if ("data.table" %in% class(codon_data)) {
    codon_data <- as.data.frame(codon_data)
  }
  
  # Get codon columns (exclude Gene_name and Expression_Group)
  codon_cols <- setdiff(names(codon_data), c("Gene_name", "Expression_Group"))
  
  # Initialize results
  aa_summary <- data.frame()
  
  for (aa in unique(genetic_code)) {
    if (aa == "STOP") next
    
    # Get all codons for this amino acid
    codons_for_aa <- names(genetic_code)[genetic_code == aa]
    codons_for_aa <- codons_for_aa[codons_for_aa %in% codon_cols]
    
    # Skip Met and Trp (non-synonymous)
    if (length(codons_for_aa) <= 1) next
    
    # Get preferred codons for this amino acid
    preferred_for_aa <- intersect(codons_for_aa, preferred_codons)
    unpreferred_for_aa <- setdiff(codons_for_aa, preferred_codons)
    
    # Count total occurrences
    if (length(preferred_for_aa) > 0) {
      preferred_count <- sum(rowSums(codon_data[, preferred_for_aa, drop = FALSE], na.rm = TRUE), na.rm = TRUE)
    } else {
      preferred_count <- 0
    }
    
    if (length(unpreferred_for_aa) > 0) {
      unpreferred_count <- sum(rowSums(codon_data[, unpreferred_for_aa, drop = FALSE], na.rm = TRUE), na.rm = TRUE)
    } else {
      unpreferred_count <- 0
    }
    
    total_count <- preferred_count + unpreferred_count
    
    # Calculate proportion
    prop_preferred <- ifelse(total_count > 0, preferred_count / total_count, NA)
    
    aa_summary <- rbind(aa_summary,
                        data.frame(
                          Amino_Acid = aa,
                          N_synonymous = length(codons_for_aa),
                          Preferred_codons = paste(preferred_for_aa, collapse = ","),
                          Preferred_count = preferred_count,
                          Unpreferred_count = unpreferred_count,
                          Total_count = total_count,
                          Prop_preferred = prop_preferred,
                          stringsAsFactors = FALSE
                        ))
  }
  
  return(aa_summary)
}