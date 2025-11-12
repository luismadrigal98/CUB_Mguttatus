# Test Suite for CUB_g_test.R
#
# Tests the G-test for codon usage bias
# Tests all 4 modes: by_gene, by_genome, by_aminoacid, heterogeneity_per_aa

library(testthat)
library(data.table)

# Source the function
source("../src/CUB_g_test.R")

context("G-test for Codon Usage Bias")

# Create standard genetic code for testing
genetic_code <- c(
  "AAA" = "Lys", "AAG" = "Lys",
  "GAA" = "Glu", "GAG" = "Glu",
  "CAA" = "Gln", "CAG" = "Gln",
  "AAC" = "Asn", "AAT" = "Asn",
  "GCT" = "Ala", "GCC" = "Ala", "GCA" = "Ala", "GCG" = "Ala",
  "TTA" = "Leu", "TTG" = "Leu", "CTT" = "Leu", "CTC" = "Leu", "CTA" = "Leu", "CTG" = "Leu"
)

# Test 1: Mode by_gene works and detects bias
test_that("G-test by_gene mode detects extreme bias", {
  # Extreme bias: only one codon per amino acid
  codon_counts <- data.table(
    Gene_name = "biased_gene",
    AAA = 100, AAG = 0,  # Only AAA used
    GAA = 100, GAG = 0,  # Only GAA used
    CAA = 100, CAG = 0   # Only CAA used
  )
  
  result <- CUB_g_test(codon_counts, genetic_code, mode = "by_gene")
  
  expect_true("Gene_name" %in% names(result))
  expect_true("G_statistic" %in% names(result))
  expect_true("df" %in% names(result))
  expect_true("p_value" %in% names(result))
  expect_true("significant" %in% names(result))
  
  # Extreme bias should have high G-statistic and low p-value
  expect_true(result$G_statistic > 100)
  expect_true(result$p_value < 0.001)
  expect_true(result$significant)
})

# Test 2: Mode by_gene detects no bias
test_that("G-test by_gene mode detects no bias for equal usage", {
  # Equal usage: no bias
  codon_counts <- data.table(
    Gene_name = "unbiased_gene",
    AAA = 50, AAG = 50,
    GAA = 50, GAG = 50,
    CAA = 50, CAG = 50
  )
  
  result <- CUB_g_test(codon_counts, genetic_code, mode = "by_gene")
  
  # Equal usage should have low G-statistic and high p-value
  expect_true(result$G_statistic < 1)
  expect_true(result$p_value > 0.05)
  expect_false(result$significant)
})

# Test 3: Mode by_gene handles multiple genes with FDR correction
test_that("G-test by_gene applies FDR correction for multiple genes", {
  codon_counts <- data.table(
    Gene_name = c("gene1", "gene2", "gene3"),
    AAA = c(100, 50, 75), 
    AAG = c(0, 50, 25),
    GAA = c(100, 50, 80), 
    GAG = c(0, 50, 20)
  )
  
  result <- CUB_g_test(codon_counts, genetic_code, mode = "by_gene", correct.p_values = TRUE)
  
  expect_equal(nrow(result), 3)
  expect_true("p_value_adj" %in% names(result))
  
  # gene1 is biased, gene2 is not
  expect_true(result[Gene_name == "gene1", significant])
  expect_false(result[Gene_name == "gene2", significant])
})

# Test 4: Mode by_genome pools all genes
test_that("G-test by_genome mode pools all genes correctly", {
  codon_counts <- data.table(
    Gene_name = c("gene1", "gene2", "gene3"),
    AAA = c(30, 40, 30),  # Total: 100
    AAG = c(10, 5, 5),    # Total: 20
    GAA = c(25, 35, 20),  # Total: 80
    GAG = c(15, 10, 10)   # Total: 35
  )
  
  result <- CUB_g_test(codon_counts, genetic_code, mode = "by_genome")
  
  expect_equal(nrow(result), 1)
  expect_true("G_statistic" %in% names(result))
  expect_true(is.numeric(result$G_statistic))
  expect_true(result$df > 0)
})

# Test 5: Mode by_aminoacid tests each AA separately
test_that("G-test by_aminoacid mode tests each amino acid", {
  codon_counts <- data.table(
    Gene_name = c("gene1", "gene2"),
    AAA = c(100, 100),  # Lys - biased
    AAG = c(10, 10),
    GAA = c(50, 50),    # Glu - not biased
    GAG = c(50, 50),
    GCT = c(25, 25),    # Ala - 4-fold degenerate
    GCC = c(25, 25),
    GCA = c(25, 25),
    GCG = c(25, 25)
  )
  
  result <- CUB_g_test(codon_counts, genetic_code, mode = "by_aminoacid", correct.p_values = TRUE)
  
  expect_true("Amino_acid" %in% names(result))
  expect_true(nrow(result) >= 2)  # At least Lys, Glu, maybe Ala
  
  # Lysine should be significant (biased)
  lys_row <- result[Amino_acid == "Lys"]
  expect_true(nrow(lys_row) > 0)
  expect_true(lys_row$significant)
  
  # Glutamic acid should not be significant (equal usage)
  glu_row <- result[Amino_acid == "Glu"]
  expect_true(nrow(glu_row) > 0)
  expect_false(glu_row$significant)
})

# Test 6: Mode heterogeneity_per_aa detects consistent vs variable usage
test_that("G-test heterogeneity_per_aa detects variation across genes", {
  # Lysine: gene1 prefers AAA, gene2 prefers AAG (heterogeneous)
  # Glutamic acid: both genes prefer GAA equally (homogeneous)
  codon_counts <- data.table(
    Gene_name = c("gene1", "gene2", "gene3"),
    AAA = c(90, 10, 90),   # Heterogeneous across genes
    AAG = c(10, 90, 10),
    GAA = c(80, 80, 80),   # Homogeneous across genes
    GAG = c(20, 20, 20)
  )
  
  result <- CUB_g_test(codon_counts, genetic_code, mode = "heterogeneity_per_aa", correct.p_values = TRUE)
  
  expect_true("Amino_acid" %in% names(result))
  expect_true("G_heterogeneity" %in% names(result))
  
  # Lysine should show significant heterogeneity
  lys_row <- result[Amino_acid == "Lys"]
  if (nrow(lys_row) > 0) {
    expect_true(lys_row$G_heterogeneity > result[Amino_acid == "Glu", G_heterogeneity])
  }
})

# Test 7: G-statistic is always non-negative
test_that("G-statistic is always non-negative", {
  codon_counts <- data.table(
    Gene_name = c("gene1", "gene2"),
    AAA = c(60, 40),
    AAG = c(40, 60),
    GAA = c(70, 30),
    GAG = c(30, 70)
  )
  
  for (mode in c("by_gene", "by_genome", "by_aminoacid")) {
    result <- CUB_g_test(codon_counts, genetic_code, mode = mode)
    expect_true(all(result$G_statistic >= 0),
                info = paste("Mode:", mode))
  }
})

# Test 8: Degrees of freedom calculation is correct
test_that("Degrees of freedom are calculated correctly", {
  # 2-fold degenerate: df = 1
  codon_counts <- data.table(
    Gene_name = "test",
    AAA = 50, AAG = 50
  )
  
  result <- CUB_g_test(codon_counts, genetic_code, mode = "by_genome")
  expect_equal(result$df, 1)
  
  # 4-fold degenerate: df = 3
  codon_counts <- data.table(
    Gene_name = "test",
    GCT = 25, GCC = 25, GCA = 25, GCG = 25
  )
  
  result <- CUB_g_test(codon_counts, genetic_code, mode = "by_genome")
  expect_equal(result$df, 3)
})

# Test 9: Function handles missing codons gracefully
test_that("G-test handles missing codons (zero counts)", {
  codon_counts <- data.table(
    Gene_name = "sparse_gene",
    AAA = 100, AAG = 0,  # One codon not used
    GAA = 0, GAG = 0     # Entire amino acid not used
  )
  
  # Should not throw error
  expect_error({
    result <- CUB_g_test(codon_counts, genetic_code, mode = "by_gene")
  }, NA)
})

# Test 10: Output is data.table
test_that("G-test returns data.table object", {
  codon_counts <- data.table(
    Gene_name = "test",
    AAA = 60, AAG = 40
  )
  
  result <- CUB_g_test(codon_counts, genetic_code, mode = "by_genome")
  
  expect_true(inherits(result, "data.table"))
})

cat("✓ All G-test tests passed!\n\n")
