splitInPartsAux <- function(string, size)
{
  #' Auxiliar function to split a string in sub-strings of fixed length
  #' 
  #' @description This function can be used to detect correct reading frames 
  #' (number of bp should be a multiple of 3) and to organize the transcripts
  #' in the constituting codons.
  #' 
  #' @param string String that is going to be split
  #' @param size Size of fragments after spliting
  #' 
  #' @return A vector of fragments
  #' _________________________________________________________________________
  
  pat <- paste0('.{1,', size, '}')
  unlist(stri_extract_all_regex(string, pat))
}