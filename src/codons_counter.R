codons_counter <- function(sequence, gene, codons)
{
  #' Main function to count how many times a given codon appears in the seq input
  #' 
  #' @param sequence Sequence for which the codons are going to be quantified.
  #' This is a string, that is going to be split using `splitInPartsAux`.
  #' @param gene Gene name to be stored in the output data table.
  #' @param codons Codons to quantify.
  #' 
  #' @return data.frame entry with the Gene name and the counts for each codon
  #' ___________________________________________________________________________
  #' 
  #' Counts are returned as INTEGER for every codon, including absent ones.
  #' The previous implementation used `ifelse(is.na(seq_table[x]), 0, ...)`,
  #' where the `0` literal is a double while a table lookup is an integer. A
  #' codon therefore came back double whenever it was absent from a gene, so
  #' after rbind across genes a column was integer only if every gene contained
  #' that codon. That mixed typing is what makes melt.data.table warn that
  #' 'measure.vars' are "not all of the same type". Values were correct either
  #' way; this keeps the storage type consistent (and is vectorised, so it also
  #' drops the per-codon loop).
  
  sequence <- splitInPartsAux(sequence, 3)
  seq_table <- table(sequence)
  
  counts <- stats::setNames(rep(0L, length(codons)), codons)
  present <- intersect(names(seq_table), codons)
  if (length(present)) counts[present] <- as.integer(seq_table[present])
  
  result <- data.frame(Gene_name = gene, stringsAsFactors = FALSE)
  result[codons] <- as.list(counts)
  
  return(result)
}
