PSPM <- function(syn_codons, rf) 
{
  #' Function to calculate the Position-Specific Probability Matrix
  #' 
  #' @param syn_codons Vector of synonymous codons
  #' @param rf Named vector of relative frequencies (must sum to 1)
  #' 
  #' @return PSPM, with nucleotides (rows) and positions (cols)
  #' ___________________________________________________________________________
  
  codons_m <- t(sapply(syn_codons, function(X) {
    unlist(strsplit(x = X, split = ""))
  }))
  
  # Set column names
  colnames(codons_m) <- 1:3 # Positions in triplet
  
  # Get the PSPM
  PSPM <- matrix(0, nrow = 4, ncol = 3) # Initialize with 0, not NA
  rownames(PSPM) <- c("A", "C", "G", "T")
  colnames(PSPM) <- 1:3
  
  # Bug 3 Fix (Replaced sapply with for loops):
  # This is the clear, standard R way to populate a matrix.
  for (pos in 1:3) {
    
    # Get all nucleotides at this position (e.g., c("C","C","C","C","A","A"))
    nuc_at_pos <- codons_m[, pos] 
    
    for (base in c("A", "C", "G", "T")) {
      
      # Find which codons have this 'base' at this 'pos'
      codons_with_base <- syn_codons[nuc_at_pos == base]
      
      # Sum the relative frequencies of those specific codons
      if (length(codons_with_base) > 0) {
        PSPM[base, pos] <- sum(rf[codons_with_base])
      }
      # (Else, it remains 0)
    }
  }
  
  return(PSPM)
}