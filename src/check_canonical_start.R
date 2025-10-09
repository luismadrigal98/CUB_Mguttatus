check_canonical_start <- function(transcript_set, 
                                  start = 'ATG')
{
  #' Function to check if the initial codon in a transcript is canonical (equal
  #' to start)
  #' 
  #' @param transcript_set Object of class DNAStringSet
  #' @param start Canonical start codon (ATG by default)
  #' _________________________________________________________________________
  
  assertthat::assert_that(class(transcript_set) == "DNAStringSet",
                          msg = "Input object (transcript_set) must be of class `DNAStringSet`")
  
  selector <- sapply(1:length(transcript_set), function(i)
  {
    as.character(trans[[i]][1:3]) == start
  })
  
  return(selector)
}