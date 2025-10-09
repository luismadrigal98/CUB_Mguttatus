gene_name_extractor <- function(name_string, 
                                field_sep = " ",
                                gene_name_pos = 1)
{
  #' This function is used to extract the gene names from a name string (names
  #' as stored in a DNAStringSet)
  #' 
  #' @param name_string String that contains the name of the genetic element
  #' @param field_sep Separator used to build the string
  #' @param gene_name_pos Position where the gene name appears (usually the first one)
  #' 
  #' @return Cleaned name
  #' ___________________________________________________________________________
  
  if(!is.character(name_string))
  {
    stop("Gene name string must be of class character")
  }
  
  return(unlist(strsplit(x = name_string, split = field_sep))[gene_name_pos])
}