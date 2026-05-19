check_cds <- function(transcript_set,
                      stops = c("TAA", "TAG", "TGA"))
{
  #' Validate coding sequences.
  #'
  #' Keeps transcripts whose length is a multiple of 3, that end on a stop
  #' codon, and that contain no internal stop codons. Used immediately after
  #' check_canonical_start() to produce a clean CDS DNAStringSet.
  #'
  #' @param transcript_set DNAStringSet of (canonical-start) CDS sequences.
  #' @param stops Stop codons (default standard genetic code).
  #' @return DNAStringSet of validated CDS.

  assertthat::assert_that(
    inherits(transcript_set, "DNAStringSet"),
    msg = "Input object (transcript_set) must be of class `DNAStringSet`"
  )

  if (length(transcript_set) == 0L) return(transcript_set)

  widths <- BiocGenerics::width(transcript_set)
  ok_len <- (widths %% 3L) == 0L

  last_codon <- as.character(Biostrings::subseq(
    transcript_set,
    start = pmax(widths - 2L, 1L),
    end   = widths
  ))
  ok_stop <- ok_len & (last_codon %in% stops)

  cds_no_stop <- Biostrings::subseq(transcript_set,
                                    start = 1L,
                                    end   = pmax(widths - 3L, 0L))
  no_internal_stop <- vapply(seq_along(cds_no_stop), function(i) {
    s <- cds_no_stop[[i]]
    if (length(s) < 3L) return(TRUE)
    codons <- as.character(Biostrings::codons(s))
    !any(codons %in% stops)
  }, logical(1))

  keep <- ok_len & ok_stop & no_internal_stop
  n_drop <- sum(!keep)
  if (n_drop > 0L) {
    message(sprintf(
      "check_cds(): dropped %d / %d transcripts (length %% 3, terminal stop, no internal stop).",
      n_drop, length(keep)
    ))
  }
  transcript_set[keep]
}
