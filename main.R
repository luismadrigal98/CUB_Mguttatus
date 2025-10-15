##' @title CUB in Mimulus guttaus
##' 
##' @author Luis Javier Madrigal-Roca & John K. Kelly
##' ____________________________________________________________________________

## *****************************************************************************
## 1) Load required libraries and set up environment ----
## _____________________________________________________________________________

# Source the set_environment function first
source("./src/set_environment.R")

required_libraries <- c('data.table', 'Biostrings', 'assertthat', 
                        'stringi', 'foreach', 'doParallel',
                        'doFuture', 'ggplot2')

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
## 2) Source all helper functions ----
## _____________________________________________________________________________

source("./src/splitInPartsAux.R")
source("./src/gene_name_extractor.R")
source("./src/check_canonical_start.R")
source("./src/codons_counter.R")
source("./src/codon_quant.R")
source("./src/calculate_rscu.R")
source("./src/calculate_enc.R")
source("./src/calculate_gc_content.R")
source("./src/visualize_codon_usage.R")
source("./src/neutrality_analysis.R")
source("./src/cub_summary.R")

## *****************************************************************************
## 3) Set work directory ----
## _____________________________________________________________________________

setwd(".")

## *****************************************************************************
## 4) Load the data ----
## _____________________________________________________________________________

## 4.1) Analysis from transcript (if available is a shortcut) ----

trans <- Biostrings::readDNAStringSet(filepath = "./data/Mguttatusvar_IM767_887_v2.1.cds_primaryTranscriptOnly.fa", 
                                      format = 'fasta')

codon_usage <- codon_quant(trans, codons = names(genetic_code_dna_long), 
                           parallel = T)

## *****************************************************************************
## 5) Comprehensive CUB Analysis ----
## _____________________________________________________________________________

message("Performing comprehensive codon usage bias analysis...")

# Run complete analysis and generate all outputs
cub_results <- cub_summary(codon_usage, genetic_code_dna_long, 
                          output_dir = "./results")

# Create amino acid-specific codon logos
create_aa_specific_logos(codon_usage, genetic_code_dna_long,
                         output_dir = "./results/codon_logos")

## *****************************************************************************
## 6) Additional analyses (optional) ----
## _____________________________________________________________________________

## Individual metric calculations (if needed separately):

# Calculate RSCU
# rscu_values <- calculate_rscu(codon_usage, genetic_code_dna_long)

# Calculate ENC
# enc_values <- calculate_enc(codon_usage, genetic_code_dna_long)

# Calculate GC content
# gc_content <- calculate_gc_content(codon_usage)

# Create specific visualizations
# visualize_codon_usage(codon_usage, genetic_code_dna_long, 
#                      "custom_heatmap.pdf", type = "heatmap")
# neutrality_plot(gc_content, "custom_neutrality.pdf")
# enc_plot(enc_values, gc_content, "custom_enc.pdf")
# pr2_bias_plot(codon_usage, "custom_pr2.pdf")

# Create logo for specific amino acid
# create_codon_logo(codon_usage, genetic_code_dna_long, "Leu", "leucine_logo.pdf")

message("\nAnalysis complete! Check the './results' directory for all outputs.")

## 4.2) Analysis from fasta and gff3 ----