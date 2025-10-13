##' @title CUB in Mimulus guttaus
##' 
##' @author Luis Javier Madrigal-Roca & John K. Kelly
##' ____________________________________________________________________________

## *****************************************************************************
## 1) Load required libraries and set up environment ----
## _____________________________________________________________________________

required_libraries <- c('data.table', 'Biostrings', 'assertthat', 
                        'stringi', 'foreach', 'doParallel', 'coRdon',
                        'doFuture')

set_environment(required_pckgs = required_libraries, personal_seed = 1998, 
                parallel_backend = T, n_cores = 10)

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

codon_usage <- codon_quant(trans, codons = names(genetic_code_dna_long), 
                           parallel = T)

## 3.2) Analysis from fasta and gff3 ----