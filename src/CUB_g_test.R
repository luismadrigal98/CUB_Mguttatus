CUB_g_test <- function(codon_counts, genetic_code, mode = "by_gene",
                       correct.p_values = T)
{
  #' G-test for codon usage bias
  #' 
  #' @description Performs likelihood ratio test (G-test) to assess whether
  #' codon usage deviates significantly from equal usage within each amino acid.
  #' Tests the null hypothesis that all synonymous codons are used equally.
  #' 
  #' @param codon_counts Data table with codon counts per gene
  #' @param genetic_code Named vector mapping codons to amino acids
  #' @param mode "by_gene" to perform test per gene (default).
  #' "by_genome" to perform test on pooled counts across all genes.
  #' "by_aminoacid" to perform test per amino acid across all genes.
  #' "heterogeneity_per_aa" to assess heterogeneity of codon usage per amino acid across genes.
  #' @param correct.p_values Logical, whether to apply multiple testing correction
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
  
  # ============================================================================
  # Mode 1: By Gene - Test each gene individually
  # ============================================================================
  
  if(mode == "by_gene")
  {
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
    
    # Apply multiple testing correction if requested
    if(correct.p_values)
    {
      results$p_value_adj <- p.adjust(results$p_value, method = "BH")
      results$significant <- results$p_value_adj < 0.05
    }
    
    return(results)
  }
  
  # ============================================================================
  # Mode 2: By Genome - Pool all genes together
  # ============================================================================
  
  if(mode == "by_genome")
  {
    # Sum counts across all genes
    pooled_counts <- colSums(dt[, ..codon_cols])
    
    G_total <- 0
    df_total <- 0
    
    # Test each amino acid
    for(aa in unique(codon_to_aa$AA))
    {
      aa_codons <- codon_to_aa[AA == aa, Codon]
      observed <- pooled_counts[aa_codons]
      
      total <- sum(observed)
      if(total == 0) next
      
      n_codons <- length(observed)
      expected <- rep(total / n_codons, n_codons)
      
      obs_nonzero <- observed > 0
      if(sum(obs_nonzero) > 0)
      {
        G_aa <- 2 * sum(observed[obs_nonzero] * log(observed[obs_nonzero] / expected[obs_nonzero]))
        G_total <- G_total + G_aa
      }
      
      df_total <- df_total + (n_codons - 1)
    }
    
    p_value <- pchisq(G_total, df = df_total, lower.tail = FALSE)
    
    results <- data.table(
      Test = "Genome_wide",
      G_statistic = G_total,
      df = df_total,
      p_value = p_value,
      significant = p_value < 0.05
    )
    
    return(results)
  }
  
  # ============================================================================
  # Mode 3: By Amino Acid - Test each amino acid across all genes
  # ============================================================================
  
  if(mode == "by_aminoacid")
  {
    results_list <- list()
    
    for(aa in unique(codon_to_aa$AA))
    {
      aa_codons <- codon_to_aa[AA == aa, Codon]
      
      # Pool counts for this amino acid across all genes
      aa_counts <- colSums(dt[, ..aa_codons])
      
      total <- sum(aa_counts)
      if(total == 0) next
      
      n_codons <- length(aa_codons)
      expected <- rep(total / n_codons, n_codons)
      
      obs_nonzero <- aa_counts > 0
      G <- 0
      if(sum(obs_nonzero) > 0)
      {
        G <- 2 * sum(aa_counts[obs_nonzero] * log(aa_counts[obs_nonzero] / expected[obs_nonzero]))
      }
      
      df <- n_codons - 1
      p_value <- pchisq(G, df = df, lower.tail = FALSE)
      
      results_list[[aa]] <- data.table(
        Amino_Acid = aa,
        N_codons = n_codons,
        Total_count = total,
        G_statistic = G,
        df = df,
        p_value = p_value,
        significant = p_value < 0.05
      )
    }
    
    results <- rbindlist(results_list)
    
    # Apply multiple testing correction if requested
    if(correct.p_values)
    {
      results$p_value_adj <- p.adjust(results$p_value, method = "BH")
      results$significant <- results$p_value_adj < 0.05
    }
    
    return(results)
  }
  
  # ============================================================================
  # Mode 4: Heterogeneity per AA - Test if codon usage varies across genes
  # ============================================================================
  
  if(mode == "heterogeneity_per_aa")
  {
    results_list <- list()
    
    for(aa in unique(codon_to_aa$AA))
    {
      aa_codons <- codon_to_aa[AA == aa, Codon]
      n_codons <- length(aa_codons)
      
      # Get counts for this AA across all genes
      aa_data <- dt[, c("Gene_name", aa_codons), with = FALSE]
      
      # Calculate total G (pooled test)
      pooled_counts <- colSums(aa_data[, ..aa_codons])
      total_pooled <- sum(pooled_counts)
      
      if(total_pooled == 0) next
      
      expected_pooled <- rep(total_pooled / n_codons, n_codons)
      obs_nonzero_pooled <- pooled_counts > 0
      
      G_pooled <- 0
      if(sum(obs_nonzero_pooled) > 0)
      {
        G_pooled <- 2 * sum(pooled_counts[obs_nonzero_pooled] * 
                           log(pooled_counts[obs_nonzero_pooled] / expected_pooled[obs_nonzero_pooled]))
      }
      
      # Calculate G for each gene
      G_total <- 0
      n_genes_with_data <- 0
      
      for(j in 1:nrow(aa_data))
      {
        gene_counts <- as.numeric(aa_data[j, ..aa_codons])
        total_gene <- sum(gene_counts)
        
        if(total_gene == 0) next
        
        n_genes_with_data <- n_genes_with_data + 1
        expected_gene <- rep(total_gene / n_codons, n_codons)
        obs_nonzero_gene <- gene_counts > 0
        
        if(sum(obs_nonzero_gene) > 0)
        {
          G_gene <- 2 * sum(gene_counts[obs_nonzero_gene] * 
                           log(gene_counts[obs_nonzero_gene] / expected_gene[obs_nonzero_gene]))
          G_total <- G_total + G_gene
        }
      }
      
      # Heterogeneity G = G_total - G_pooled
      G_heterogeneity <- G_total - G_pooled
      
      # df = (n_genes - 1) * (n_codons - 1)
      df_heterogeneity <- (n_genes_with_data - 1) * (n_codons - 1)
      
      if(df_heterogeneity > 0 && G_heterogeneity >= 0)
      {
        p_value <- pchisq(G_heterogeneity, df = df_heterogeneity, lower.tail = FALSE)
        
        results_list[[aa]] <- data.table(
          Amino_Acid = aa,
          N_codons = n_codons,
          N_genes = n_genes_with_data,
          G_pooled = G_pooled,
          G_total = G_total,
          G_heterogeneity = G_heterogeneity,
          df = df_heterogeneity,
          p_value = p_value,
          significant = p_value < 0.05
        )
      }
    }
    
    results <- rbindlist(results_list)
    
    # Apply multiple testing correction if requested
    if(correct.p_values)
    {
      results$p_value_adj <- p.adjust(results$p_value, method = "BH")
      results$significant <- results$p_value_adj < 0.05
    }
    
    return(results)
  }
  
  stop("Invalid mode. Choose: 'by_gene', 'by_genome', 'by_aminoacid', or 'heterogeneity_per_aa'")
}