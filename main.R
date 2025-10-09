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
  #' not start with the canonical ATG.
  #' 
  #' @return Count data table of codon per gene
  #' ___________________________________________________________________________
  
  assertthat::assert_that(class(transcripts) == "DNAStringSet",
                          msg = "Input object (transcripts) must be of class `DNAStringSet`")
  
  # Check that the gene has a canonical start ATG
  
  splitInPartsAux <- function(string, size)
  {
    #' Auxiliar function to split a string in sub-strings of fixed length
    #' _________________________________________________________________________
    
    pat <- paste0('.{1,', size, '}')
    stri_extract_all_regex(string, pat)
  }
}

## 3.2) Analysis from fasta and gff3 ----