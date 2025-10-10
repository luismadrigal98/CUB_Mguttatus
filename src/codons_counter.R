codons_counter <- function(sequence, gene, codons)
{
  #' Main function to count how many times a give codon appears in the seq input
  #' 
  #' @param seq Sequence for which the codons are going to be quantified. This is
  #' a string, that is ging to be split using `splitInPartsAx`.
  #' @param gene Gene name to be stored in the output data table.
  #' @param codons Codons to quantify.
  #' 
  #' @return data.table entry wiht the Gene name and the counts for each codon
  #' ___________________________________________________________________________
  
  sequence <- splitInPartsAux(sequence, 3)
  
  # Generate the counts
  seq_table <- table(sequence)
  
  # Store the results
  result <- data.frame(Gene_name = gene)
  
  # Create the other fields
  for (x in codons)
  {
    result[, x] <- ifelse(is.na(seq_table[x]), 0, seq_table[x])
  }
  
  return(result)
}