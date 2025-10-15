##' @title Test CUB Functions with Synthetic Data
##' 
##' @description This script tests the CUB analysis functions with synthetic data
##' to ensure they work correctly without requiring the full dataset.
##' ____________________________________________________________________________

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(Biostrings)
  library(ggplot2)
})

# Source all functions
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

cat("Testing CUB analysis functions...\n\n")

# Test 1: splitInPartsAux
cat("Test 1: splitInPartsAux\n")
test_seq <- "ATGGCATAA"
result <- splitInPartsAux(test_seq, 3)
cat("Input:", test_seq, "\n")
cat("Output:", paste(result, collapse = ", "), "\n")
stopifnot(length(result) == 3)
cat("PASSED\n\n")

# Test 2: gene_name_extractor
cat("Test 2: gene_name_extractor\n")
test_name <- "GENE001 protein coding"
result <- gene_name_extractor(test_name)
cat("Input:", test_name, "\n")
cat("Output:", result, "\n")
stopifnot(result == "GENE001")
cat("PASSED\n\n")

# Test 3: Create synthetic transcript data
cat("Test 3: Creating synthetic transcript data\n")
# Create some synthetic sequences with codon bias
seq1 <- "ATGTTCTTCGCCGCCGCCTAA"  # Met-Phe-Phe-Ala-Ala-Ala-STOP
seq2 <- "ATGCTGCTGCTGGGTGGTTAA"  # Met-Leu-Leu-Leu-Gly-Gly-STOP
seq3 <- "ATGAAAAAAGACTCTTCTTAA"  # Met-Lys-Lys-Thr-Ser-Ser-STOP

synthetic_trans <- DNAStringSet(c(seq1, seq2, seq3))
names(synthetic_trans) <- c("GENE001", "GENE002", "GENE003")
cat("Created", length(synthetic_trans), "synthetic transcripts\n")
cat("PASSED\n\n")

# Test 4: check_canonical_start
cat("Test 4: check_canonical_start\n")
result <- check_canonical_start(synthetic_trans)
cat("Canonical starts:", sum(result), "/", length(result), "\n")
stopifnot(sum(result) == 3)  # All should start with ATG
cat("PASSED\n\n")

# Test 5: codons_counter
cat("Test 5: codons_counter\n")
result <- codons_counter(as.character(synthetic_trans[[1]]), "GENE001", 
                         names(genetic_code_dna_long))
cat("Gene:", result$Gene_name, "\n")
cat("Total codons counted:", sum(result[, -1]), "\n")
stopifnot(result$Gene_name == "GENE001")
cat("PASSED\n\n")

# Test 6: codon_quant
cat("Test 6: codon_quant\n")
codon_counts <- codon_quant(synthetic_trans, names(genetic_code_dna_long), 
                            parallel = FALSE)
cat("Dimensions:", nrow(codon_counts), "genes x", ncol(codon_counts), "columns\n")
cat("Genes:", paste(codon_counts$Gene_name, collapse = ", "), "\n")
stopifnot(nrow(codon_counts) == 3)
cat("PASSED\n\n")

# Test 7: calculate_rscu
cat("Test 7: calculate_rscu\n")
rscu_values <- calculate_rscu(codon_counts, genetic_code_dna_long)
cat("Dimensions:", nrow(rscu_values), "genes x", ncol(rscu_values), "columns\n")
stopifnot(nrow(rscu_values) == 3)
cat("PASSED\n\n")

# Test 8: calculate_enc
cat("Test 8: calculate_enc\n")
enc_values <- calculate_enc(codon_counts, genetic_code_dna_long)
cat("Dimensions:", nrow(enc_values), "genes x", ncol(enc_values), "columns\n")
cat("ENC values:", paste(round(enc_values$ENC, 2), collapse = ", "), "\n")
stopifnot(nrow(enc_values) == 3)
stopifnot(all(enc_values$ENC >= 20 & enc_values$ENC <= 61))
cat("PASSED\n\n")

# Test 9: calculate_gc_content
cat("Test 9: calculate_gc_content\n")
gc_content <- calculate_gc_content(codon_counts)
cat("Dimensions:", nrow(gc_content), "genes x", ncol(gc_content), "columns\n")
cat("Mean GC:", round(mean(gc_content$GC), 3), "\n")
stopifnot(nrow(gc_content) == 3)
stopifnot(all(gc_content$GC >= 0 & gc_content$GC <= 1))
cat("PASSED\n\n")

# Test 10: Test with larger synthetic dataset for plotting
cat("Test 10: Creating larger synthetic dataset for visualization tests\n")
set.seed(42)
n_genes <- 100

# Generate random sequences with varying codon bias
generate_random_cds <- function(n_codons = 50) {
  codons <- names(genetic_code_dna_long)
  codons <- codons[genetic_code_dna_long != "STOP"]
  
  # Start with ATG
  seq <- "ATG"
  
  # Add random codons with some bias
  for(i in 1:(n_codons - 2)) {
    seq <- paste0(seq, sample(codons, 1))
  }
  
  # End with stop
  seq <- paste0(seq, "TAA")
  
  return(seq)
}

large_synthetic_trans <- DNAStringSet(sapply(1:n_genes, function(i) {
  generate_random_cds(sample(20:100, 1))
}))
names(large_synthetic_trans) <- paste0("GENE", sprintf("%04d", 1:n_genes))

# Run codon quantification
large_codon_counts <- codon_quant(large_synthetic_trans, 
                                   names(genetic_code_dna_long), 
                                   parallel = FALSE)
cat("Created dataset with", nrow(large_codon_counts), "genes\n")
cat("PASSED\n\n")

# Test 11: Full analysis pipeline
cat("Test 11: Running full analysis pipeline\n")
dir.create("./test_results", showWarnings = FALSE, recursive = TRUE)

tryCatch({
  # Calculate all metrics
  large_enc <- calculate_enc(large_codon_counts, genetic_code_dna_long)
  large_gc <- calculate_gc_content(large_codon_counts)
  
  cat("  - ENC calculation: OK\n")
  cat("  - GC content calculation: OK\n")
  
  # Create visualizations
  visualize_codon_usage(large_codon_counts, genetic_code_dna_long, 
                       "./test_results/test_heatmap.pdf", type = "heatmap")
  cat("  - Heatmap creation: OK\n")
  
  visualize_codon_usage(large_codon_counts, genetic_code_dna_long, 
                       "./test_results/test_barplot.pdf", type = "barplot")
  cat("  - Barplot creation: OK\n")
  
  neutrality_plot(large_gc, "./test_results/test_neutrality.pdf")
  cat("  - Neutrality plot: OK\n")
  
  enc_plot(large_enc, large_gc, "./test_results/test_enc.pdf")
  cat("  - ENC plot: OK\n")
  
  pr2_bias_plot(large_codon_counts, "./test_results/test_pr2.pdf")
  cat("  - PR2 plot: OK\n")
  
  # Create a codon logo
  create_codon_logo(large_codon_counts, genetic_code_dna_long, "Leu", 
                   "./test_results/test_leucine_logo.pdf")
  cat("  - Codon logo: OK\n")
  
  cat("PASSED\n\n")
  
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  stop(e)
})

cat("\n======================\n")
cat("ALL TESTS PASSED!\n")
cat("======================\n")
cat("\nTest outputs saved to ./test_results/\n")
