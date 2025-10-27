# Test Suite for calculate_enc.R
#
# Tests the Effective Number of Codons (ENC) calculation
# ENC ranges from 20 (extreme bias) to 61 (no bias)

library(testthat)
library(data.table)

# Source the function
source("../src/calculate_enc.R")

context("ENC Calculation Tests")

# Create standard genetic code for testing
genetic_code <- c(
  # Leucine (6-fold)
  "TTA" = "Leu", "TTG" = "Leu", "CTT" = "Leu", "CTC" = "Leu", "CTA" = "Leu", "CTG" = "Leu",
  # Serine (6-fold)
  "TCT" = "Ser", "TCC" = "Ser", "TCA" = "Ser", "TCG" = "Ser", "AGT" = "Ser", "AGC" = "Ser",
  # Arginine (6-fold)
  "CGT" = "Arg", "CGC" = "Arg", "CGA" = "Arg", "CGG" = "Arg", "AGA" = "Arg", "AGG" = "Arg",
  # Isoleucine (3-fold)
  "ATT" = "Ile", "ATC" = "Ile", "ATA" = "Ile",
  # Valine (4-fold)
  "GTT" = "Val", "GTC" = "Val", "GTA" = "Val", "GTG" = "Val",
  # Alanine (4-fold)
  "GCT" = "Ala", "GCC" = "Ala", "GCA" = "Ala", "GCG" = "Ala",
  # Glycine (4-fold)
  "GGT" = "Gly", "GGC" = "Gly", "GGA" = "Gly", "GGG" = "Gly",
  # Proline (4-fold)
  "CCT" = "Pro", "CCC" = "Pro", "CCA" = "Pro", "CCG" = "Pro",
  # Threonine (4-fold)
  "ACT" = "Thr", "ACC" = "Thr", "ACA" = "Thr", "ACG" = "Thr",
  # Lysine (2-fold)
  "AAA" = "Lys", "AAG" = "Lys",
  # Asparagine (2-fold)
  "AAC" = "Asn", "AAT" = "Asn",
  # Glutamic acid (2-fold)
  "GAA" = "Glu", "GAG" = "Glu",
  # Aspartic acid (2-fold)
  "GAC" = "Asp", "GAT" = "Asp",
  # Glutamine (2-fold)
  "CAA" = "Gln", "CAG" = "Gln",
  # Histidine (2-fold)
  "CAC" = "His", "CAT" = "His",
  # Tyrosine (2-fold)
  "TAC" = "Tyr", "TAT" = "Tyr",
  # Cysteine (2-fold)
  "TGC" = "Cys", "TGT" = "Cys",
  # Phenylalanine (2-fold)
  "TTC" = "Phe", "TTT" = "Phe",
  # Met, Trp (1-fold each - excluded from ENC calculation body)
  "ATG" = "Met", "TGG" = "Trp",
  # Stop codons
  "TAA" = "STOP", "TAG" = "STOP", "TGA" = "STOP"
)

# Test 1: ENC is within theoretical bounds
test_that("ENC values are within theoretical bounds [20, 61]", {
  # Diverse codon usage across multiple amino acid families
  codon_counts <- data.table(
    Gene_name = "test_gene",
    # 2-fold degenerate
    AAA = 10, AAG = 15,
    GAA = 12, GAG = 9,
    # 3-fold degenerate (Ile)
    ATT = 8, ATC = 10, ATA = 5,
    # 4-fold degenerate
    GCT = 5, GCC = 7, GCA = 6, GCG = 4,
    # 6-fold degenerate
    TTA = 3, TTG = 5, CTT = 4, CTC = 6, CTA = 3, CTG = 7
  )
  
  result <- calculate_enc(codon_counts, genetic_code)
  
  expect_true(result$ENC >= 20,
              info = "ENC should be >= 20 (maximum bias)")
  expect_true(result$ENC <= 61,
              info = "ENC should be <= 61 (no bias)")
})

# Test 2: No bias case (ENC should approach 61)
test_that("ENC approaches 61 for uniform codon usage", {
  # Equal usage of all synonymous codons
  codon_counts <- data.table(
    Gene_name = "no_bias_gene",
    # 6-fold degenerate - all equal
    TTA = 10, TTG = 10, CTT = 10, CTC = 10, CTA = 10, CTG = 10,
    TCT = 10, TCC = 10, TCA = 10, TCG = 10, AGT = 10, AGC = 10,
    CGT = 10, CGC = 10, CGA = 10, CGG = 10, AGA = 10, AGG = 10,
    # 4-fold degenerate - all equal
    GTT = 10, GTC = 10, GTA = 10, GTG = 10,
    GCT = 10, GCC = 10, GCA = 10, GCG = 10,
    # 3-fold degenerate
    ATT = 10, ATC = 10, ATA = 10,
    # 2-fold degenerate - all equal
    AAA = 10, AAG = 10,
    GAA = 10, GAG = 10
  )
  
  result <- calculate_enc(codon_counts, genetic_code)
  
  expect_true(result$ENC > 55,
              info = "ENC should be close to 61 for uniform usage")
})

# Test 3: Extreme bias case (ENC should approach 20)
test_that("ENC approaches 20 for extreme codon bias", {
  # Only one codon per amino acid family
  codon_counts <- data.table(
    Gene_name = "biased_gene",
    # 6-fold - only one used
    TTA = 100, TTG = 0, CTT = 0, CTC = 0, CTA = 0, CTG = 0,
    TCT = 100, TCC = 0, TCA = 0, TCG = 0, AGT = 0, AGC = 0,
    CGT = 100, CGC = 0, CGA = 0, CGG = 0, AGA = 0, AGG = 0,
    # 4-fold - only one used
    GTT = 100, GTC = 0, GTA = 0, GTG = 0,
    GCT = 100, GCC = 0, GCA = 0, GCG = 0,
    # 3-fold - only one used
    ATT = 100, ATC = 0, ATA = 0,
    # 2-fold - only one used
    AAA = 100, AAG = 0,
    GAA = 100, GAG = 0
  )
  
  result <- calculate_enc(codon_counts, genetic_code)
  
  expect_true(result$ENC < 35,
              info = "ENC should be low for extreme bias")
})

# Test 4: Output structure
test_that("ENC output has correct structure", {
  codon_counts <- data.table(
    Gene_name = "test_gene",
    AAA = 10, AAG = 15,
    GAA = 12, GAG = 9,
    ATT = 8, ATC = 10, ATA = 5,  # 3-fold degenerate
    GCT = 5, GCC = 7, GCA = 6, GCG = 4,  # 4-fold degenerate
    TTA = 3, TTG = 5, CTT = 4, CTC = 6, CTA = 3, CTG = 7  # 6-fold degenerate
  )
  
  result <- calculate_enc(codon_counts, genetic_code)
  
  expect_true("Gene_name" %in% names(result))
  expect_true("ENC" %in% names(result))
  expect_equal(nrow(result), 1)
  expect_equal(result$Gene_name, "test_gene")
  expect_true(is.numeric(result$ENC))
})

# Test 5: Multiple genes
test_that("ENC works with multiple genes", {
  codon_counts <- data.table(
    Gene_name = c("gene1", "gene2"),
    AAA = c(10, 20),
    AAG = c(15, 5),
    GAA = c(12, 18),
    GAG = c(9, 7),
    ATT = c(8, 12),
    ATC = c(10, 8),
    ATA = c(5, 10),
    GCT = c(5, 7),
    GCC = c(7, 5),
    GCA = c(6, 8),
    GCG = c(4, 6),
    TTA = c(3, 4),
    TTG = c(5, 6),
    CTT = c(4, 5),
    CTC = c(6, 7),
    CTA = c(3, 4),
    CTG = c(7, 8)
  )
  
  result <- calculate_enc(codon_counts, genetic_code)
  
  expect_equal(nrow(result), 2)
  expect_equal(result$Gene_name, c("gene1", "gene2"))
  expect_true(all(result$ENC >= 20))
  expect_true(all(result$ENC <= 61))
})

# Test 6: ENC is numeric and finite
test_that("ENC values are numeric and finite", {
  codon_counts <- data.table(
    Gene_name = "test_gene",
    AAA = 10, AAG = 15,
    ATT = 8, ATC = 10, ATA = 5,
    GCT = 5, GCC = 7, GCA = 6, GCG = 4,
    TTA = 3, TTG = 5, CTT = 4, CTC = 6, CTA = 3, CTG = 7  # 6-fold
  )
  
  result <- calculate_enc(codon_counts, genetic_code)
  
  expect_true(is.numeric(result$ENC))
  expect_true(is.finite(result$ENC))
  expect_false(is.na(result$ENC))
})

cat("✓ All ENC calculation tests passed!\n\n")
