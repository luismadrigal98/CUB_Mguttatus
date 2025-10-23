calculate_enc <- function(codon_counts, genetic_code)
{
  #' Calculate Effective Number of Codons (ENC)
  #' 
  #' @description ENC quantifies codon bias, ranging from 20 (extreme bias, 
  #' one codon per amino acid) to 61 (no bias, all codons used equally).
  #' Based on Wright (1990) method.
  #' 
  #' @param codon_counts Data table with codon counts per gene
  #' @param genetic_code Named vector mapping codons to amino acids
  #' 
  #' @return Data frame with Gene_name and ENC values
  #' ___________________________________________________________________________
  
  library(data.table)
  
  # Group codons by amino acid
  aa_groups <- split(names(genetic_code), genetic_code)
  aa_groups <- aa_groups[names(aa_groups) != "STOP" & 
                           names(aa_groups) != "Trp" &
                           names(aa_groups) != "Met"]
  
  # Classify amino acids by number of synonymous codons
  aa_by_n <- list(
    "2" = names(aa_groups)[sapply(aa_groups, length) == 2],
    "3" = names(aa_groups)[sapply(aa_groups, length) == 3],
    "4" = names(aa_groups)[sapply(aa_groups, length) == 4],
    "6" = names(aa_groups)[sapply(aa_groups, length) == 6]
  )
  
  enc_values <- numeric(nrow(codon_counts))
  
  for(gene_idx in 1:nrow(codon_counts))
  {
    F_values <- list()
    
    # Calculate F for each degeneracy class
    for(n in names(aa_by_n))
    {
      n_int <- as.integer(n)
      aa_list <- aa_by_n[[n]]
      
      if(length(aa_list) == 0) next
      
      F_sum <- 0
      n_aa <- 0
      
      for(aa in aa_list)
      {
        codons <- aa_groups[[aa]]
        counts <- as.numeric(codon_counts[gene_idx, codons, with = FALSE])
        total <- sum(counts)
        
        if(total > 0)
        {
          # Calculate homozygosity
          p_squared_sum <- sum((counts / total)^2)
          F_aa <- (total * p_squared_sum - 1) / (total - 1)
          
          if(!is.na(F_aa) && is.finite(F_aa))
          {
            F_sum <- F_sum + F_aa
            n_aa <- n_aa + 1
          }
        }
      }
      
      if(n_aa > 0)
      {
        F_values[[n]] <- F_sum / n_aa
      }
      else
      {
        F_values[[n]] <- 0
      }
    }
    
    # Calculate ENC
    # ENC = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6
    enc <- 2  # Met and Trp (non-degenerate)
    
    if(!is.null(F_values[["2"]]) && F_values[["2"]] > 0)
    {
      enc <- enc + 9 / F_values[["2"]]
    }
    else
    {
      enc <- enc + 9
    }
    
    if(!is.null(F_values[["3"]]) && F_values[["3"]] > 0)
    {
      enc <- enc + 1 / F_values[["3"]]
    }
    else
    {
      enc <- enc + 1
    }
    
    if(!is.null(F_values[["4"]]) && F_values[["4"]] > 0)
    {
      enc <- enc + 5 / F_values[["4"]]
    }
    else
    {
      enc <- enc + 5
    }
    
    if(!is.null(F_values[["6"]]) && F_values[["6"]] > 0)
    {
      enc <- enc + 3 / F_values[["6"]]
    }
    else
    {
      enc <- enc + 3
    }
    
    enc_values[gene_idx] <- enc
  }
  
  result <- data.frame(
    Gene_name = codon_counts$Gene_name,
    ENC = enc_values
  )
  
  return(result)
}
