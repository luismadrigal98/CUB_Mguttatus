trim_uninformative <- function(quant_table, genetic_code)
{
  #' Function to remove codons that code STOP, Trp, or Met from a quantification
  #' table.
  #' 
  #' @param quant_table Data table where each column represent a codon, and each
  #' row a gene. There are two accepted formats, data.table and matrix(for 
  #' compatibility with cubar)
  #' 
  #' @return Trimmed data table.
  #' ___________________________________________________________________________
  
  codons_for_trimming <- names(genetic_code[genetic_code %in% c("STOP", "Trp", "Met")])
  
  if(is.data.table(quant_table))
    {
      quant_table[, (codons_for_trimming) := NULL]
  }
  
  else
  {
    try(quant_table[, !colnames(quant_table) %in% codons_for_trimming])   
  }
}