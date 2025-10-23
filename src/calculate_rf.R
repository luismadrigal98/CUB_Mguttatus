calculate_rf <- function(codon_counts, genetic_code)
{
  #' Calculate Relative Frequency of each synonymous codon
  #' 
  #' @description RF is the proportion of usage of a given codon from a set of
  #' synonynmous one. It ranges from 0 to 1.
  #' 
  #' @param codon_counts Data table with codon counts per gene (from codon_quant)
  #' @param genetic_code Named vector mapping codons to amino acids
  #' 
  #' @return Data table with RF values per gene
  #' ___________________________________________________________________________
  
  library(data.table)
  
  # Get codon columns (exclude Gene_name)
  codon_cols <- setdiff(names(codon_counts), "Gene_name")
  
  # Calculate RSCU for each gene
  rf_results <- copy(codon_counts)
  
  # Group codons by amino acid
  aa_groups <- split(names(genetic_code), genetic_code)
  
  # Remove STOP codons from analysis
  aa_groups <- aa_groups[names(aa_groups) != "STOP"]
  
  for(aa in names(aa_groups))
  {
    # Get the set of synonymous codons for this AA
    syn_codons <- aa_groups[[aa]]
    
    # 5. Perform data.table operations on ALL genes at once
    
    # A) Calculate T (total_aa_count) for every gene
    #    rowSums() is applied to the subset of data (.SD) defined by .SDcols
    rf_results[, total_aa_count := rowSums(.SD), .SDcols = syn_codons]
    
    # B) Calculate RF = X_i / T for each synonymous codon
    #    This is a fast loop over a few column NAMES, not rows.
    for(codon in syn_codons)
    {
      # We use fifelse() for a fast, safe division by zero.
      # If expected_freq is 0, RSCU is 0, otherwise calculate X_i / E_i
      rf_results[, (codon) := fifelse(total_aa_count == 0, 
                                        0, 
                                        .SD[[1]] / total_aa_count), 
                   .SDcols = codon]
    }
  }
  
  # Clean up the temporary helper columns
  rf_results[, c("total_aa_count") := NULL]
  
  return(rf_results)
}