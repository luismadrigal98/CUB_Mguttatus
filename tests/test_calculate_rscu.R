##' Unit tests for RSCU calculation
##' 
##' Tests the calculate_rscu function for mathematical correctness
##' ____________________________________________________________________________

library(testthat)
library(data.table)

# Source the function
source("../src/calculate_rscu.R")

context("RSCU Calculation Tests")

# Test 1: RSCU definition - mean should be 1.0
test_that("Mean RSCU within amino acid family equals 1.0", {
  codon_counts <- data.table(
    Gene_name = "test_gene",
    # Leucine with bias
    TTA = 10, TTG = 20, CTT = 30, CTC = 20, CTA = 10, CTG = 10
  )
  
  genetic_code <- c(
    "TTA" = "Leu", "TTG" = "Leu", "CTT" = "Leu", 
    "CTC" = "Leu", "CTA" = "Leu", "CTG" = "Leu"
  )
  
  rscu_result <- calculate_rscu(codon_counts, genetic_code)
  
  # Get RSCU values for Leucine codons
  leu_codons <- c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG")
  rscu_values <- as.numeric(rscu_result[1, leu_codons, with = FALSE])
  
  # Mean RSCU should be 1.0 (or very close due to floating point)
  expect_equal(mean(rscu_values), 1.0, tolerance = 1e-10)
})

# Test 2: Uniform usage should give RSCU = 1.0 for all codons
test_that("RSCU equals 1.0 for uniform codon usage", {
  codon_counts <- data.table(
    Gene_name = "uniform_gene",
    # All Leu codons used equally
    TTA = 10, TTG = 10, CTT = 10, CTC = 10, CTA = 10, CTG = 10,
    # All Ser codons used equally  
    TCT = 5, TCC = 5, TCA = 5, TCG = 5, AGT = 5, AGC = 5
  )
  
  genetic_code <- c(
    "TTA" = "Leu", "TTG" = "Leu", "CTT" = "Leu", 
    "CTC" = "Leu", "CTA" = "Leu", "CTG" = "Leu",
    "TCT" = "Ser", "TCC" = "Ser", "TCA" = "Ser", 
    "TCG" = "Ser", "AGT" = "Ser", "AGC" = "Ser"
  )
  
  rscu_result <- calculate_rscu(codon_counts, genetic_code)
  
  all_codons <- names(genetic_code)
  rscu_values <- as.numeric(rscu_result[1, all_codons, with = FALSE])
  
  # All RSCU values should be 1.0
  expect_true(all(abs(rscu_values - 1.0) < 1e-10))
})

# Test 3: Extreme bias - one codon used exclusively
test_that("RSCU handles extreme codon bias correctly", {
  codon_counts <- data.table(
    Gene_name = "biased_gene",
    # Only TTA used for Leucine (6 codons)
    TTA = 60, TTG = 0, CTT = 0, CTC = 0, CTA = 0, CTG = 0
  )
  
  genetic_code <- c(
    "TTA" = "Leu", "TTG" = "Leu", "CTT" = "Leu", 
    "CTC" = "Leu", "CTA" = "Leu", "CTG" = "Leu"
  )
  
  rscu_result <- calculate_rscu(codon_counts, genetic_code)
  
  # TTA should have RSCU = 6.0 (6 codons for Leu, only 1 used)
  expect_equal(rscu_result$TTA, 6.0, tolerance = 1e-10)
  
  # Other codons should have RSCU = 0
  expect_equal(rscu_result$TTG, 0.0)
  expect_equal(rscu_result$CTT, 0.0)
})

# Test 4: RSCU with 2-codon family
test_that("RSCU correct for 2-codon amino acids", {
  codon_counts <- data.table(
    Gene_name = "test_gene",
    # Lysine: 2 codons with 3:1 ratio
    AAA = 30, AAG = 10
  )
  
  genetic_code <- c("AAA" = "Lys", "AAG" = "Lys")
  
  rscu_result <- calculate_rscu(codon_counts, genetic_code)
  
  # AAA used 3x more, so RSCU should be 1.5
  # AAG used 1x, so RSCU should be 0.5
  expect_equal(rscu_result$AAA, 1.5, tolerance = 1e-10)
  expect_equal(rscu_result$AAG, 0.5, tolerance = 1e-10)
  
  # Mean should still be 1.0
  expect_equal(mean(c(rscu_result$AAA, rscu_result$AAG)), 1.0, tolerance = 1e-10)
})

# Test 5: Multiple genes
test_that("RSCU works with multiple genes independently", {
  codon_counts <- data.table(
    Gene_name = c("gene1", "gene2"),
    TTA = c(10, 20), TTG = c(10, 0)  # Leu codons
  )
  
  genetic_code <- c("TTA" = "Leu", "TTG" = "Leu")
  
  rscu_result <- calculate_rscu(codon_counts, genetic_code)
  
  # Gene 1: equal usage -> both RSCU = 1.0
  expect_equal(rscu_result$TTA[1], 1.0, tolerance = 1e-10)
  expect_equal(rscu_result$TTG[1], 1.0, tolerance = 1e-10)
  
  # Gene 2: only TTA used -> RSCU = 2.0, 0.0
  expect_equal(rscu_result$TTA[2], 2.0, tolerance = 1e-10)
  expect_equal(rscu_result$TTG[2], 0.0)
})

# Test 6: Zero counts handling
test_that("RSCU handles zero counts correctly", {
  codon_counts <- data.table(
    Gene_name = "gene_with_zeros",
    TTA = 0, TTG = 0, CTT = 10
  )
  
  genetic_code <- c("TTA" = "Leu", "TTG" = "Leu", "CTT" = "Leu")
  
  # Should not produce NaN or Inf
  rscu_result <- calculate_rscu(codon_counts, genetic_code)
  
  expect_true(all(is.finite(as.numeric(rscu_result[, c("TTA", "TTG", "CTT"), with = FALSE]))))
})

# Test 7: Output structure
test_that("RSCU output has correct structure", {
  codon_counts <- data.table(
    Gene_name = c("gene1", "gene2"),
    TTA = c(10, 20), ATG = c(5, 10)
  )
  
  genetic_code <- c("TTA" = "Leu", "ATG" = "Met")
  
  rscu_result <- calculate_rscu(codon_counts, genetic_code)
  
  expect_true("Gene_name" %in% names(rscu_result))
  expect_true("TTA" %in% names(rscu_result))
  expect_true("ATG" %in% names(rscu_result))
  expect_equal(nrow(rscu_result), 2)
})

# Test 8: RSCU values are non-negative
test_that("RSCU values are always non-negative", {
  codon_counts <- data.table(
    Gene_name = "test_gene",
    TTA = 5, TTG = 10, CTT = 15, CTC = 20
  )
  
  genetic_code <- c("TTA" = "Leu", "TTG" = "Leu", 
                   "CTT" = "Leu", "CTC" = "Leu")
  
  rscu_result <- calculate_rscu(codon_counts, genetic_code)
  
  codon_cols <- c("TTA", "TTG", "CTT", "CTC")
  rscu_values <- as.numeric(rscu_result[1, codon_cols, with = FALSE])
  
  expect_true(all(rscu_values >= 0))
})

cat("✓ All RSCU calculation tests passed!\n")
