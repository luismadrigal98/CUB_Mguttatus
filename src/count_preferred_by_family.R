count_preferred_by_family <- function(codon_data, preferred_codons, families) {
  
  # Convert to data.frame if it's a data.table
  if ("data.table" %in% class(codon_data)) {
    codon_data <- as.data.frame(codon_data)
  }
  
  # Get codon columns
  codon_cols <- setdiff(names(codon_data), c("Gene_name", "Expression_Group"))
  
  # Initialize results
  family_summary <- data.frame()
  
  for (family_name in names(families)) {
    codons_in_family <- families[[family_name]]
    codons_in_family <- codons_in_family[codons_in_family %in% codon_cols]
    
    if (length(codons_in_family) == 0) next
    
    # Split into preferred vs unpreferred
    preferred_in_family <- intersect(codons_in_family, preferred_codons)
    unpreferred_in_family <- setdiff(codons_in_family, preferred_codons)
    
    # Count occurrences
    if (length(preferred_in_family) > 0) {
      preferred_count <- sum(rowSums(codon_data[, preferred_in_family, drop = FALSE], na.rm = TRUE), na.rm = TRUE)
    } else {
      preferred_count <- 0
    }
    
    if (length(unpreferred_in_family) > 0) {
      unpreferred_count <- sum(rowSums(codon_data[, unpreferred_in_family, drop = FALSE], na.rm = TRUE), na.rm = TRUE)
    } else {
      unpreferred_count <- 0
    }
    
    total_count <- preferred_count + unpreferred_count
    prop_preferred <- ifelse(total_count > 0, preferred_count / total_count, NA)
    
    # Parse family name
    aa <- sub("_.*", "", family_name)
    degeneracy <- sub(".*_", "", family_name)
    
    family_summary <- rbind(family_summary,
                            data.frame(
                              Amino_Acid = aa,
                              Degeneracy = degeneracy,
                              Family = family_name,
                              N_codons = length(codons_in_family),
                              Preferred_codons = paste(preferred_in_family, collapse = ","),
                              Preferred_count = preferred_count,
                              Unpreferred_count = unpreferred_count,
                              Total_count = total_count,
                              Prop_preferred = prop_preferred,
                              stringsAsFactors = FALSE
                            ))
  }
  
  return(family_summary)
}