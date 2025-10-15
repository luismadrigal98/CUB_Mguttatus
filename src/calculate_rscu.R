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
  
  for(gene_idx in 1:nrow(codon_counts))
  {
    for(aa in names(aa_groups))
    {
      synonymous_codons <- aa_groups[[aa]]
      
      # Get counts for this amino acid's synonymous codons
      codon_values <- as.numeric(codon_counts[gene_idx, synonymous_codons, with = FALSE])
      total_aa_count <- sum(codon_values)
      
      if(total_aa_count > 0)
      {
        n_synonymous <- length(synonymous_codons)
        expected_freq <- total_aa_count / n_synonymous
        
        # Calculate RSCU for each codon
        for(i in seq_along(synonymous_codons))
        {
          codon <- synonymous_codons[i]
          observed <- codon_values[i]
          rscu_value <- ifelse(expected_freq > 0, 
                               (observed / expected_freq) * n_synonymous / n_synonymous,
                               0)
          rscu_results[gene_idx, (codon) := observed / expected_freq]
        }
      }
      else
      {
        # If no counts, set RSCU to 0
        for(codon in synonymous_codons)
        {
          rscu_results[gene_idx, (codon) := 0]
        }
      }
    }
  }
  
  return(rscu_results)
}
