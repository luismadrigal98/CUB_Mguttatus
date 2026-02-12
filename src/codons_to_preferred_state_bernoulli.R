#' Function to decompose a transcript object of class DNAStringSet into a sequence
#' of codon states, indicating for each gene and codon position whether there is
#' a preferred codon or any of the unpreferred ones.
#' 
#' @param transcripts Object of DNAStringSet with main transcript per gene
#' @param preferred_codons Character vector of preferred codons
#' @param parallel Whether to run the analysis in parallel or not. Default is TRUE.
#' 
#' @return Data.table with three columns: GeneID, Position, Is_Preferred
#' _____________________________________________________________________________
codons_to_preferred_state_bernoulli <- function(transcripts, preferred_codons,
                                                parallel = TRUE)
{
  # Check the class of the input file using inherits for S4 safety
  assertthat::assert_that(inherits(transcripts, "DNAStringSet"), 
                          msg = "Input object (transcripts) must be of class `DNAStringSet`")
  
  # Check that preferred_codons is a vector
  assertthat::assert_that(is.character(preferred_codons), 
                          msg = "preferred_codons must be a character vector with the preferred codons")
  
  # Define the worker function to avoid code duplication and scoping issues
  process_transcript <- function(i) {
    # Extract sequence as character
    seq_char <- as.character(transcripts[[i]])
    
    # Split sequence in triplets
    seq_codons <- splitInPartsAux(seq_char, 3)
    
    # Extract gene name
    gene <- gene_name_extractor(names(transcripts)[i])
    
    # Determine preference
    is_preferred <- seq_codons %in% preferred_codons
    
    # Return a list (lighter than data.frame for iteration)
    list(GeneID = rep(gene, length(seq_codons)),
         Position = seq_along(seq_codons), 
         Is_Preferred = is_preferred)
  }
  
  # Run analysis
  if (parallel) {
    # Ensure parallel backend is registered by user before calling this function
    results_list <- foreach(i = seq_along(transcripts), 
                            .export = c("splitInPartsAux",
                                        "gene_name_extractor"),
                            .packages = c("data.table")) %dopar% {
                              process_transcript(i)
                            }
  } else {
    # lapply will use process_transcript which inherits preferred_codons from scope
    results_list <- lapply(seq_along(transcripts), process_transcript)    
  }
  
  # Efficiently combine results
  results <- data.table::rbindlist(results_list)
  data.table::setorder(results, GeneID)
  
  return(results)
}