#!/usr/bin/env Rscript
# Test CAI calculation

library(dplyr)

# Source required functions
source("./src/calculate_cai.R")

# Load genetic code
load("./data/genetic_code_dna_long.rda")

# Create simple test data
# 3 genes with different codon usage patterns
test_codon_usage <- data.frame(
  Gene_name = c("Gene1", "Gene2", "Gene3"),
  # Alanine codons (GCT, GCC, GCA, GCG)
  GCT = c(10, 2, 5),   # Gene1 prefers GCT
  GCC = c(1, 8, 5),    # Gene2 prefers GCC
  GCA = c(1, 1, 5),
  GCG = c(1, 1, 5),
  # Leucine codons (TTA, TTG, CTT, CTC, CTA, CTG)
  TTA = c(15, 1, 5),   # Gene1 prefers TTA
  TTG = c(1, 1, 5),
  CTT = c(1, 1, 5),
  CTC = c(1, 15, 5),   # Gene2 prefers CTC
  CTA = c(1, 1, 5),
  CTG = c(1, 1, 5),
  stringsAsFactors = FALSE
)

cat("Test codon usage data:\n")
print(test_codon_usage)

# Use Gene1 as reference (highly expressed)
reference_genes <- c("Gene1")

cat("\n\nCalculating CAI with Gene1 as reference...\n")
cat("Expected: Gene1 should have CAI = 1.0 (it IS the reference)\n")
cat("          Gene2 should have lower CAI (uses different codons)\n")
cat("          Gene3 should have intermediate CAI (uniform usage)\n\n")

# Calculate CAI
cai_results <- calculate_cai(
  codon_counts = test_codon_usage,
  reference_genes = reference_genes,
  genetic_code = genetic_code_dna_long
)

cat("\n\nFinal CAI values:\n")
print(cai_results$cai_values)

cat("\n\nRelative adaptiveness (w) values:\n")
w_relevant <- cai_results$w_table |>
  filter(codon %in% c("GCT", "GCC", "GCA", "GCG", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG"))
print(w_relevant)

cat("\n\nTest completed successfully!\n")
