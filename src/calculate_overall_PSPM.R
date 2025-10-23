calculate_overall_PSPM <- function(rf_data, genetic_code)
{
  #' Function to calculate the overall Position-Specific Probability Matrix (PSPM)
  #' across all genes in the dataset.
  #' 
  #' @param rf_data Data frame with relative frequencies (RF) of codons per gene.
  #' @param genetic_code Data frame defining the genetic code with codons and their
  #' corresponding amino acids.
  #' 
  #' @return List of overall PSPM matrices per amino acid.
  #' ___________________________________________________________________________
  
  # Calculate average RF across all genes for each codon
  rf_avg <- rf_data |>
    dplyr::select(-Gene_name) |>
    colMeans()
  
  # Iterate over each amino acid in the genetic code
  unique_aas <- unique(genetic_code)
  
  overall_PSPM <- lapply(unique_aas, function(aa)
  {
    # Get synonymous codons for the amino acid
    syn_codons <- names(genetic_code)[genetic_code == aa]
    
    # Get RF values for these codons
    rf_values <- rf_avg[syn_codons]
    
    # Re-normalize
    total_rf_for_aa <- sum(rf_values)
    if (total_rf_for_aa > 0) {
      rf_values <- rf_values / total_rf_for_aa
    }
    
    # Calculate PSPM for these synonymous codons
    PSPM(syn_codons, rf_values)
  })
  
  names(overall_PSPM) <- unique_aas
  
  return(overall_PSPM)
}