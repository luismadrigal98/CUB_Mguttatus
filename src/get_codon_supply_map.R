get_codon_supply_map <- function(trna_counts)
{
  #' Create a mapping of codon-to-tRNA supply based on wobble rules
  #'
  #' @param trna_counts A data.table with 'Anticodon' and 'tRNA_count'
  #' @return A data.table with 'Codon' and 'tRNA_supply'
  #' ___________________________________________________________________________

  # Standard DNA-based wobble rules for anticodon 1st base -> codon 3rd base
  # (Anticodon T = RNA U)
  wobble_rules <- list(
    "G" = c("T", "C"),  # Anticodon G pairs with Codon T or C
    "C" = c("G"),       # Anticodon C pairs with Codon G
    "A" = c("T"),       # Anticodon A (often modified to I) pairs with T (U)
    "T" = c("A", "G")   # Anticodon T (U) pairs with Codon A or G
  )
  
  # Standard complement for bases 2 and 3
  complement <- c("A" = "T", "T" = "A", "G" = "C", "C" = "G")
  
  # This will store a list of all recognized codons and their tRNA source
  wobble_map_long <- list()
  
  for (i in 1:nrow(trna_counts)) {
    anticodon <- trna_counts$Anticodon[i]
    count <- trna_counts$tRNA_count[i]
    
    # Split anticodon
    ac_1 <- substr(anticodon, 1, 1) # Wobble base
    ac_2 <- substr(anticodon, 2, 2)
    ac_3 <- substr(anticodon, 3, 3)
    
    # Get codon bases (1 and 2 are reverse-complemented)
    codon_1 <- complement[ac_3]
    codon_2 <- complement[ac_2]
    
    # Get all possible codon 3rd bases from wobble rules
    codon_3_list <- wobble_rules[[ac_1]]
    
    if (is.null(codon_3_list)) next # Skip if anticodon is weird (e.g., 'N')
    
    # Create a row for each codon this anticodon can read
    for (c3 in codon_3_list) {
      codon <- paste0(codon_1, codon_2, c3)
      wobble_map_long[[length(wobble_map_long) + 1]] <- data.table(
        Codon = codon,
        tRNA_count = count
      )
    }
  }
  
  # Bind all rows
  wobble_map_dt <- rbindlist(wobble_map_long)
  
  # Sum the counts for codons that are read by multiple anticodons
  codon_supply <- wobble_map_dt[, .(tRNA_supply = sum(tRNA_count)), by = Codon]
  
  return(codon_supply)
}