codon_quant <- function(transcripts, codons, parallel = T,
                        check_canonical = T)
{
  #' This function will process the entries of a DNAStringSet an build a 
  #' quantification data.table, where each gene (row) will have associated a 
  #' count for each possible codon (64 columns).
  #' 
  #' @param transcripts DNAStringSet with the primary transcript per gene.
  #' @param codons Codons to quantify.
  #' @param parallel Whether to enable or not the parallel processing.
  #' @param check_canonical This flag enable the filtering out of genes that do
  #' not start with the canonical ATG. If TRUE is passed, ATG is assumed to be
  #' canonical.
  #' 
  #' @return Count data table of codon per gene
  #' ___________________________________________________________________________
  
  assertthat::assert_that(class(transcripts) == "DNAStringSet",
                          msg = "Input object (transcripts) must be of class `DNAStringSet`")
  
  # Check that the reading frame is correct (length of transcript is multiple of 3)
  transcripts <- transcripts[sapply(1:length(transcripts), function(i){
    length(splitInPartsAux(as.character(transcripts[[i]]), 1)) %% 3 == 0
  })]
  
  if(parallel)
  {
    results <- foreach(i = 1:length(transcripts), 
                       .export = c("splitInPartsAux",
                                   "codons_counter",
                                   "codons"),
                       .packages = c("data.table"),
                       .combine = 'rbind') %dopar%
      {
        codons_counter(sequence = as.character(transcripts[[i]]),
                       gene = gene_name_extractor(names(transcripts[i])),
                       codons = codons)
      }
    
    results <- results |> as.data.table() |> setorder(Gene_name)
  }
  
  else # Sequential approach (better for debugging)
  {
    results <- lapply(X = 1:length(transcripts), 
                      FUN = function(i)
                      {
                        codons_counter(sequence = as.character(transcripts[[i]]),
                                       gene = gene_name_extractor(names(transcripts[i])),
                                       codons = codons)
                      }
    )
    
    results <- do.call("rbind", results)
    
    results <- results |> as.data.table() |> setorder(Gene_name)
  }
  
  return(results)
}