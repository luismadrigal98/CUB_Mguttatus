calculate_rscu <- function(codon_counts, genetic_code)
{
  #' Calculate Relative Synonymous Codon Usage (RSCU)
  #' 
  #' @description RSCU is the ratio of the observed frequency of a codon to the 
  #' expected frequency if all synonymous codons for an amino acid were used 
  #' equally. RSCU > 1 indicates positive bias, < 1 negative bias.
  #' 
  #' @param codon_counts Data table with codon counts per gene (from codon_quant)
  #' @param genetic_code Named vector mapping codons to amino acids
  #' 
  #' @return Data table with RSCU values per gene
  #' ___________________________________________________________________________
  
  library(data.table)
  
  # Get codon columns (exclude Gene_name)
  codon_cols <- setdiff(names(codon_counts), "Gene_name")
  
  # Calculate RSCU for each gene
  rscu_results <- copy(codon_counts)
  
  # Group codons by amino acid
  aa_groups <- split(names(genetic_code), genetic_code)
  
  # Remove STOP codons from analysis
  aa_groups <- aa_groups[names(aa_groups) != "STOP"]
  
  for(aa in names(aa_groups))
  {
    # Get the set of synonymous codons for this AA
    syn_codons <- aa_groups[[aa]]
    n_synonymous <- length(syn_codons)
    
    # 5. Perform data.table operations on ALL genes at once
    
    # A) Calculate T (total_aa_count) for every gene
    #    rowSums() is applied to the subset of data (.SD) defined by .SDcols
    rscu_results[, total_aa_count := rowSums(.SD), .SDcols = syn_codons]
    
    # B) Calculate E (expected_freq) for every gene
    rscu_results[, expected_freq := total_aa_count / n_synonymous]
    
    # C) Calculate RSCU = X_i / E_i for each synonymous codon
    #    This is a fast loop over a few column NAMES, not rows.
    for(codon in syn_codons)
    {
      # We use fifelse() for a fast, safe division by zero.
      # If expected_freq is 0, RSCU is 0, otherwise calculate X_i / E_i
      rscu_results[, (codon) := fifelse(expected_freq == 0, 
                                        0, 
                                        .SD[[1]] / expected_freq), 
                   .SDcols = codon]
    }
  }
  
  # Clean up the temporary helper columns
  rscu_results[, c("total_aa_count", "expected_freq") := NULL]
  
  return(rscu_results)
}