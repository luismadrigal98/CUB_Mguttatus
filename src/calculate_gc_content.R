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
  
  codon_cols <- setdiff(names(codon_counts), "Gene_name")
  
  # Initialize results
  results <- data.frame(
    Gene_name = codon_counts$Gene_name,
    GC = numeric(nrow(codon_counts)),
    GC1 = numeric(nrow(codon_counts)),
    GC2 = numeric(nrow(codon_counts)),
    GC3 = numeric(nrow(codon_counts)),
    GC12 = numeric(nrow(codon_counts)),
    GC3s = numeric(nrow(codon_counts))
  )
  
  for(gene_idx in 1:nrow(codon_counts))
  {
    # Get counts for all codons
    total_gc <- 0
    total_at <- 0
    
    gc1 <- 0
    at1 <- 0
    gc2 <- 0
    at2 <- 0
    gc3 <- 0
    at3 <- 0
    gc3s <- 0
    at3s <- 0
    
    for(codon in codon_cols)
    {
      count <- as.numeric(codon_counts[gene_idx, codon, with = FALSE])
      
      if(count > 0)
      {
        bases <- strsplit(codon, "")[[1]]
        
        # Position 1
        if(bases[1] %in% c("G", "C"))
        {
          gc1 <- gc1 + count
          total_gc <- total_gc + count
        }
        else
        {
          at1 <- at1 + count
          total_at <- total_at + count
        }
        
        # Position 2
        if(bases[2] %in% c("G", "C"))
        {
          gc2 <- gc2 + count
          total_gc <- total_gc + count
        }
        else
        {
          at2 <- at2 + count
          total_at <- total_at + count
        }
        
        # Position 3
        if(bases[3] %in% c("G", "C"))
        {
          gc3 <- gc3 + count
          total_gc <- total_gc + count
        }
        else
        {
          at3 <- at3 + count
          total_at <- total_at + count
        }
        
        # Position 3 synonymous (exclude ATG and TGG)
        if(!(codon %in% c("ATG", "TGG")))
        {
          if(bases[3] %in% c("G", "C"))
          {
            gc3s <- gc3s + count
          }
          else
          {
            at3s <- at3s + count
          }
        }
      }
    }
    
    total_bases <- total_gc + total_at
    total1 <- gc1 + at1
    total2 <- gc2 + at2
    total3 <- gc3 + at3
    total12 <- total1 + total2
    total3s <- gc3s + at3s
    
    results$GC[gene_idx] <- ifelse(total_bases > 0, total_gc / total_bases, 0)
    results$GC1[gene_idx] <- ifelse(total1 > 0, gc1 / total1, 0)
    results$GC2[gene_idx] <- ifelse(total2 > 0, gc2 / total2, 0)
    results$GC3[gene_idx] <- ifelse(total3 > 0, gc3 / total3, 0)
    results$GC12[gene_idx] <- ifelse(total12 > 0, (gc1 + gc2) / total12, 0)
    results$GC3s[gene_idx] <- ifelse(total3s > 0, gc3s / total3s, 0)
  }
  
  return(results)
}
