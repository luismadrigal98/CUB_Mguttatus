CUB_g_test <- function(codon_counts, genetic_code)
{
  #' G-test for codon usage bias
  #' 
  #' @description Performs likelihood ratio test (G-test) to assess whether
  #' codon usage deviates significantly from equal usage within each amino acid.
  #' Tests the null hypothesis that all synonymous codons are used equally.
  #' 
  #' @param codon_counts Data table with codon counts per gene
  #' @param genetic_code Named vector mapping codons to amino acids
  #' 
  #' @return Data frame with G-test results per gene:
  #' - Gene_name: Gene identifier
  #' - G_statistic: G-test statistic
  #' - df: Degrees of freedom
  #' - p_value: P-value from chi-square distribution
  #' - significant: Logical, TRUE if p < 0.05
  #' ___________________________________________________________________________
  
  require(data.table)
  
  # Convert to data.table
  dt <- as.data.table(codon_counts)
  
  # Get codon columns
  codon_cols <- setdiff(names(dt), "Gene_name")
  
  # Create codon to AA mapping
  codon_to_aa <- data.table(
    Codon = names(genetic_code),
    AA = as.character(genetic_code)
  )
  
  # Remove STOP codons and non degenerated amino acids
  codon_to_aa <- codon_to_aa[AA != "STOP" & AA != "Trp" & AA != "Met"]
  
  # Count synonymous codon family sizes
  aa_degeneracy <- codon_to_aa[, .(n_codons = .N), by = AA]
  
  # Merge with codon to AA mapping
  codon_to_aa <- merge(codon_to_aa, aa_degeneracy, by = "AA")
  
  # Only keep amino acids with > 1 codon (synonymous)
  codon_to_aa <- codon_to_aa[n_codons > 1]
  
  # Initialize results
  results <- data.table(
    Gene_name = dt$Gene_name,
    G_statistic = numeric(nrow(dt)),
    df = integer(nrow(dt)),
    p_value = numeric(nrow(dt)),
    significant = logical(nrow(dt))
  )
  
  # Calculate G-test for each gene
  for(i in 1:nrow(dt))
  {
    gene_name <- dt$Gene_name[i]
    
    # Get codon counts for this gene
    gene_counts <- as.numeric(dt[i, ..codon_cols])
    names(gene_counts) <- codon_cols
    
    # G statistic and df for this gene
    G <- 0
    df_total <- 0
    
    # Group by amino acid
    for(aa in unique(codon_to_aa$AA))
    {
      # Get codons for this amino acid
      aa_codons <- codon_to_aa[AA == aa, Codon]
      
      # Get observed counts
      observed <- gene_counts[aa_codons]
      
      # Skip if no counts for this amino acid
      total <- sum(observed)
      if(total == 0) next
      
      # Expected counts (equal usage)
      n_codons <- length(observed)
      expected <- rep(total / n_codons, n_codons)
      
      # Calculate G statistic contribution
      # G = 2 * sum(O * ln(O/E))
      # Only include non-zero observed counts
      obs_nonzero <- observed > 0
      if(sum(obs_nonzero) > 0)
      {
        G_aa <- 2 * sum(observed[obs_nonzero] * log(observed[obs_nonzero] / expected[obs_nonzero]))
        G <- G + G_aa
      }
      
      # Degrees of freedom = number of codons - 1
      df_total <- df_total + (n_codons - 1)
    }
    
    # Store results
    results$G_statistic[i] <- G
    results$df[i] <- df_total
    
    # Calculate p-value from chi-square distribution
    if(df_total > 0 && G >= 0)
    {
      results$p_value[i] <- pchisq(G, df = df_total, lower.tail = FALSE)
      results$significant[i] <- results$p_value[i] < 0.05
    }
    else
    {
      results$p_value[i] <- NA
      results$significant[i] <- FALSE
    }
  }
  
  return(results)
}
