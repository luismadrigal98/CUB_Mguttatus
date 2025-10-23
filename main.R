##' @title CUB in Mimulus guttaus
##' 
##' @author Luis Javier Madrigal-Roca & John K. Kelly
##' ____________________________________________________________________________

## *****************************************************************************
## 2) Set work directory ----
## _____________________________________________________________________________

setwd(".")

## *****************************************************************************
## 2) Load required libraries and set up environment ----
## _____________________________________________________________________________

# Source the set_environment function first
source("./src/set_environment.R")

required_libraries <- c('data.table', 'Biostrings', 'assertthat', 
                        'stringi', 'foreach', 'doParallel',
                        'doFuture', 'ggplot2', 'grid', 'gridExtra',
                        'ggseqlogo', 'FactoMineR',
                        'factoextra')

set_environment(required_pckgs = required_libraries, personal_seed = 1998, 
                parallel_backend = T, n_cores = 10)

# 1.1) Definition of globals ----
# Look-up table

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
## 3) Load the data ----
## _____________________________________________________________________________

## 3.1) Analysis from transcript (if available is a shortcut) ----

trans <- Biostrings::readDNAStringSet(filepath = "./data/Mguttatusvar_IM767_887_v2.1.cds_primaryTranscriptOnly.fa", 
                                      format = 'fasta')

codon_usage <- codon_quant(trans, codons = names(genetic_code_dna_long), 
                           parallel = T)

## 3.2) Clean the codon usage object (remove the STOP codon)

codon_usage <- codon_usage |>
  trim_uninformative(genetic_code = genetic_code_dna_long)

## *****************************************************************************
## 4) Comprehensive CUB Analysis ----
## _____________________________________________________________________________

message("Performing comprehensive codon usage bias analysis...")

# Run complete analysis and generate all outputs
cub_results <- cub_summary(codon_usage, genetic_code_dna_long, 
                          output_dir = "./results")

# Create amino acid-specific codon logos
create_aa_specific_logos(codon_usage, genetic_code_dna_long,
                         output_dir = "./results/codon_logos")

# Create comprehensive codon logo
create_comprehensive_codon_logo(codon_usage, genetic_code_dna_long,
                               output_file = "./results/comprehensive_codon_logo.pdf")

## *****************************************************************************
## 5) tRNA abundance correlation analysis ----
## _____________________________________________________________________________

message("Analyzing correlation between codon usage and tRNA abundance...")

# Perform tRNA-codon correlation analysis using filtered tRNA data
tRNA_correlation_results <- tRNA_codon_correlation(
  codon_counts = codon_usage,
  tRNA_file = "./data/Mguttatusvar_IM767_887_v2.0_tRNA_filtered.txt",
  genetic_code = genetic_code_dna_long,
  output_dir = "./results/tRNA_analysis",
  test_method = "spearman"  # Can also use "pearson" or "kendall"
)

message("tRNA correlation analysis complete!")

## *****************************************************************************
## 6) Additional analyses (optional) ----
## _____________________________________________________________________________

## Individual metric calculations (if needed separately):

# Calculate RSCU
rscu_values <- calculate_rscu(codon_usage, genetic_code_dna_long)

# Calculate ENC
enc_values <- calculate_enc(codon_usage, genetic_code_dna_long)

# Calculate RF
rf_values <- calculate_rf(codon_usage, genetic_code_dna_long)

# Get he PSPM
pspm_overall <- calculate_overall_PSPM(rf_values, genetic_code_dna_long)

# Create logos
create_aa_logo(pspm_overall)

# Calculate GC content
gc_content <- calculate_gc_content(codon_usage)

# Create specific visualizations
visualize_codon_usage(codon_usage, genetic_code_dna_long,
                     "custom_heatmap.pdf", type = "heatmap")
neutrality_plot(gc_content, "custom_neutrality.pdf")
enc_plot(enc_values, gc_content, "custom_enc.pdf")
pr2_bias_plot(codon_usage, "custom_pr2.pdf")

message("\nAnalysis complete! Check the './results' directory for all outputs.")