calculate_gc_content <- function(codon_counts)
{
  #' Calculate GC content metrics for codon usage bias analysis
  #' 
  #' @description Calculates GC, GC1, GC2, GC3, GC12, and GC3s content.
  #' GC3s excludes Met and Trp (non-degenerate at 3rd position).
  #' 
  #' @param codon_counts Data table with codon counts per gene
  #' 
  #' @return Data frame with Gene_name and GC content metrics
  #' ___________________________________________________________________________
  
  library(data.table)
  
  # Convert to data.table if not already
  dt <- as.data.table(codon_counts)
  
  codon_cols <- setdiff(names(dt), "Gene_name")
  
  # Pre-compute base composition for each codon (vectorized, done once)
  codon_bases <- data.table(
    Codon = codon_cols,
    Base1 = substr(codon_cols, 1, 1),
    Base2 = substr(codon_cols, 2, 2),
    Base3 = substr(codon_cols, 3, 3),
    Is_GC1 = as.integer(substr(codon_cols, 1, 1) %in% c("G", "C")),
    Is_GC2 = as.integer(substr(codon_cols, 2, 2) %in% c("G", "C")),
    Is_GC3 = as.integer(substr(codon_cols, 3, 3) %in% c("G", "C")),
    Is_Syn3 = as.integer(!(codon_cols %in% c("ATG", "TGG")))  # Exclude Met, Trp
  )
  setkey(codon_bases, Codon)
  
  # Melt the data to long format (Gene x Codon x Count)
  dt_long <- melt(dt, id.vars = "Gene_name", 
                  variable.name = "Codon", 
                  value.name = "Count",
                  variable.factor = FALSE)
  
  # Merge with codon base composition
  dt_long <- merge(dt_long, codon_bases, by = "Codon")
  
  # Calculate GC counts per gene using vectorized operations
  results <- dt_long[, .(
    # Overall GC
    GC_count = sum(Count * (Is_GC1 + Is_GC2 + Is_GC3)),
    Total_count = sum(Count * 3),
    
    # Position-specific
    GC1_count = sum(Count * Is_GC1),
    Total1 = sum(Count),
    
    GC2_count = sum(Count * Is_GC2),
    Total2 = sum(Count),
    
    GC3_count = sum(Count * Is_GC3),
    Total3 = sum(Count),
    
    # GC3s (synonymous sites only)
    GC3s_count = sum(Count * Is_GC3 * Is_Syn3),
    Total3s = sum(Count * Is_Syn3)
  ), by = Gene_name]
  
  # Calculate frequencies
  results[, `:=`(
    GC = ifelse(Total_count > 0, GC_count / Total_count, 0),
    GC1 = ifelse(Total1 > 0, GC1_count / Total1, 0),
    GC2 = ifelse(Total2 > 0, GC2_count / Total2, 0),
    GC3 = ifelse(Total3 > 0, GC3_count / Total3, 0),
    GC12 = ifelse(Total1 + Total2 > 0, (GC1_count + GC2_count) / (Total1 + Total2), 0),
    GC3s = ifelse(Total3s > 0, GC3s_count / Total3s, 0)
  )]
  
  # Return only the final columns
  return(results[, .(Gene_name, GC, GC1, GC2, GC3, GC12, GC3s)])
}
