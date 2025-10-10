##' @title CUB in Mimulus guttaus
##' 
##' @author Luis Javier Madrigal-Roca & John K. Kelly
##' ____________________________________________________________________________

## *****************************************************************************
## 1) Load required libraries ----
## _____________________________________________________________________________

library(coRdon)
library(Biostrings)
library(data.table)
library(assertthat)
library(stringi)
library(data.table)

# 1.1) Definition of globals ----

genetic_code_dna_long <- c(
  "TTT"="Phe", "TTC"="Phe", "TTA"="Leu", "TTG"="Leu",
  "TCT"="Ser", "TCC"="Ser", "TCA"="Ser", "TCG"="Ser",
  "TAT"="Tyr", "TAC"="Tyr", "TAA"="STOP", "TAG"="STOP",
  "TGT"="Cys", "TGC"="Cys", "TGA"="STOP", "TGG"="Trp",
  "CTT"="Leu", "CTC"="Leu", "CTA"="Leu", "CTG"="Leu",
  "CCT"="Pro", "CCC"="Pro", "CCA"="Pro", "CCG"="Pro",
  "CAT"="His", "CAC"="His", "CAA"="Gln", "CAG"="Gln",
  "CGT"="Arg", "CGC"="Arg", "CGA"="Arg", "CGG"="Arg",
  "ATT"="Ile", "ATC"="Ile", "ATA"="Ile", "ATG"="Met",
  "ACT"="Thr", "ACC"="Thr", "ACA"="Thr", "ACG"="Thr",
  "AAT"="Asn", "AAC"="Asn", "AAA"="Lys", "AAG"="Lys",
  "AGT"="Ser", "AGC"="Ser", "AGA"="Arg", "AGG"="Arg",
  "GTT"="Val", "GTC"="Val", "GTA"="Val", "GTG"="Val",
  "GCT"="Ala", "GCC"="Ala", "GCA"="Ala", "GCG"="Ala",
  "GAT"="Asp", "GAC"="Asp", "GAA"="Glu", "GAG"="Glu",
  "GGT"="Gly", "GGC"="Gly", "GGA"="Gly", "GGG"="Gly"
)

## *****************************************************************************
## 2) Set work directory ----
## _____________________________________________________________________________

setwd(".")

## *****************************************************************************
## 3) Load the data ----
## _____________________________________________________________________________

## 3.1) Analysis from transcript (if available is a shortcut) ----

trans <- Biostrings::readDNAStringSet(filepath = "./data/Mguttatusvar_IM767_887_v2.1.cds_primaryTranscriptOnly.fa", 
                                      format = 'fasta')

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
  
  # Check that the gene has a canonical start ATG
  filter <- check_canonical(transcripts)
  transcripts <- transcripts[filter] # Filter out non-canonical genes
  
  # Check that the reading frame is correct (length of transcript is multiple of 3)
  transcripts <- transcripts[sapply(1:length(transcripts), function(i){
    length(splitInPartsAux(as.character(transcripts[[i]]), 1)) %% 3 == 0
  })]
  
  if(parallel)
  {
    results <- foreach(i = 1:length(transcripts), 
                       .export = c("splitInPartsAux"),
                       .packages = c("data.table")) %dopar%
      {
        
      }
  }
  
  else # Sequential approach (better for debugging)
  {
    
  }
}

codons_counter <- function(seq)
{
  #' Main function to count how many times a give codon appears in the seq input
  #' 
  #' @param seq Sequence for which the codons are going to be quantified
  #' 
  #' @return data.table entry wiht the Gene name and the counts for each codon
  #' ___________________________________________________________________________
  
  
}

## 3.2) Analysis from fasta and gff3 ----