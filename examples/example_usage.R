##' @title Example Usage of CUB Analysis Pipeline
##' 
##' @description This script demonstrates how to use the CUB analysis functions
##' with step-by-step examples.
##' 
##' @author Luis Javier Madrigal-Roca & John K. Kelly
##' ____________________________________________________________________________

## Load required libraries
required_libraries <- c('data.table', 'Biostrings', 'assertthat', 
                        'stringi', 'foreach', 'doParallel',
                        'doFuture', 'ggplot2')

# Source set_environment function
source("./src/set_environment.R")
set_environment(required_pckgs = required_libraries, personal_seed = 1998, 
                parallel_backend = TRUE, n_cores = 10)

# Source all analysis functions
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

# Define genetic code
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

## =============================================================================
## EXAMPLE 1: Complete automated analysis
## =============================================================================

# Load transcript data
trans <- Biostrings::readDNAStringSet(
  filepath = "./data/Mguttatusvar_IM767_887_v2.1.cds_primaryTranscriptOnly.fa", 
  format = 'fasta'
)

# Count codons
codon_usage <- codon_quant(trans, codons = names(genetic_code_dna_long), 
                           parallel = TRUE)

# Run comprehensive analysis (generates all outputs)
cub_results <- cub_summary(codon_usage, genetic_code_dna_long, 
                          output_dir = "./results")

# Create amino acid-specific logos
create_aa_specific_logos(codon_usage, genetic_code_dna_long,
                         output_dir = "./results/codon_logos")

## =============================================================================
## EXAMPLE 2: Step-by-step analysis with individual functions
## =============================================================================

# Calculate RSCU
rscu_values <- calculate_rscu(codon_usage, genetic_code_dna_long)
print(head(rscu_values))

# Calculate ENC
enc_values <- calculate_enc(codon_usage, genetic_code_dna_long)
print(head(enc_values))
print(paste("Mean ENC:", mean(enc_values$ENC, na.rm = TRUE)))

# Calculate GC content
gc_content <- calculate_gc_content(codon_usage)
print(head(gc_content))

## =============================================================================
## EXAMPLE 3: Creating individual visualizations
## =============================================================================

# Create codon usage heatmap
visualize_codon_usage(codon_usage, genetic_code_dna_long, 
                     "./results/my_heatmap.pdf", type = "heatmap")

# Create codon usage barplot
visualize_codon_usage(codon_usage, genetic_code_dna_long, 
                     "./results/my_barplot.pdf", type = "barplot")

# Create neutrality plot (mutation vs selection)
neutrality_plot(gc_content, "./results/my_neutrality.pdf")

# Create ENC plot (identify genes under selection)
enc_plot(enc_values, gc_content, "./results/my_enc.pdf")

# Create PR2 bias plot (purine/pyrimidine bias)
pr2_bias_plot(codon_usage, "./results/my_pr2.pdf")

## =============================================================================
## EXAMPLE 4: Creating logos for specific amino acids
## =============================================================================

# Create logo for Leucine (most degenerate, 6 codons)
create_codon_logo(codon_usage, genetic_code_dna_long, "Leu", 
                 "./results/leucine_logo.pdf")

# Create logo for Serine (6 codons)
create_codon_logo(codon_usage, genetic_code_dna_long, "Ser", 
                 "./results/serine_logo.pdf")

# Create logo for Phenylalanine (2 codons)
create_codon_logo(codon_usage, genetic_code_dna_long, "Phe", 
                 "./results/phenylalanine_logo.pdf")

## =============================================================================
## EXAMPLE 5: Exploring results
## =============================================================================

# Access different components of results
summary_table <- cub_results$summary_table  # Complete data with all metrics
enc_data <- cub_results$enc_results         # ENC values only
gc_data <- cub_results$gc_results           # GC metrics only
statistics <- cub_results$statistics        # Summary statistics

# Print summary statistics
cat("\n=== Summary Statistics ===\n")
print(statistics)

# Identify genes with extreme codon bias (low ENC)
high_bias_genes <- enc_values[enc_values$ENC < 30, ]
cat("\nGenes with high codon bias (ENC < 30):\n")
print(head(high_bias_genes))

# Identify genes with low codon bias (high ENC)
low_bias_genes <- enc_values[enc_values$ENC > 55, ]
cat("\nGenes with low codon bias (ENC > 55):\n")
print(head(low_bias_genes))

# Examine GC content distribution
cat("\nGC content distribution:\n")
cat(sprintf("Mean GC: %.2f%%\n", mean(gc_content$GC) * 100))
cat(sprintf("Mean GC3s: %.2f%%\n", mean(gc_content$GC3s) * 100))
cat(sprintf("GC12 vs GC3 correlation: %.3f\n", 
            cor(gc_content$GC12, gc_content$GC3, use = "complete.obs")))

message("\n=== Analysis complete! ===")
message("All results saved to ./results/")
