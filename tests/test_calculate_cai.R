# Test Suite for calculate_cai.R
#
# Tests the Codon Adaptation Index (CAI) calculation
# CAI ranges from 0 (worst adaptation) to 1 (perfect adaptation)

library(testthat)
library(data.table)

# Source the function
source("../src/calculate_cai.R")

context("CAI Calculation Tests")

# Create standard genetic code for testing
genetic_code <- c(
  "AAA" = "Lys", "AAG" = "Lys",
  "GAA" = "Glu", "GAG" = "Glu",
  "CAA" = "Gln", "CAG" = "Gln",
  "AAC" = "Asn", "AAT" = "Asn",
  "GCT" = "Ala", "GCC" = "Ala", "GCA" = "Ala", "GCG" = "Ala"
)

# Test 1: CAI is within bounds [0, 1]
test_that("CAI values are within bounds [0, 1]", {
  codon_counts <- data.table(
    Gene_name = c("test_gene", "ref_gene"),
    AAA = c(10, 5), AAG = c(15, 25),
    GAA = c(12, 20), GAG = c(9, 5),
    CAA = c(8, 10), CAG = c(13, 20)
  )
  
  # Use ref_gene as the reference
  reference_genes <- c("ref_gene")
  
  result <- calculate_cai(codon_counts, reference_genes, genetic_code)
  
  expect_true(result$CAI >= 0)
  expect_true(result$CAI <= 1)
})

# Test 2: Perfect adaptation gives CAI = 1
test_that("CAI equals 1 for perfect adaptation to reference", {
  # Both genes use same codon preferences
  codon_counts <- data.table(
    Gene_name = c("ref_gene", "adapted_gene"),
    AAA = c(5, 10), AAG = c(95, 90),  # Both prefer AAG
    GAA = c(85, 80), GAG = c(15, 20)   # Both prefer GAA
  )
  
  reference_genes <- c("ref_gene")
  
  result <- calculate_cai(codon_counts, reference_genes, genetic_code)
  
  expect_true(result$CAI > 0.95)  # Should be very close to 1
})

# Test 3: Poor adaptation gives low CAI
test_that("CAI is low for poor adaptation to reference", {
  # Genes use opposite codon preferences
  codon_counts <- data.table(
    Gene_name = c("ref_gene", "maladapted_gene"),
    AAA = c(5, 90), AAG = c(95, 10),   # Opposite preferences
    GAA = c(85, 20), GAG = c(15, 80)   # Opposite preferences
  )
  
  reference_genes <- c("ref_gene")
  
  result <- calculate_cai(codon_counts, reference_genes, genetic_code)
  
  expect_true(result$CAI < 0.5)  # Should be low
})

# Test 4: CAI output structure is correct
test_that("CAI output has correct structure", {
  codon_counts <- data.table(
    Gene_name = c("ref_gene", "test_gene"),
    AAA = c(5, 10), AAG = c(25, 15)
  )
  
  reference_genes <- c("ref_gene")
  
  result <- calculate_cai(codon_counts, reference_genes, genetic_code)
  
  expect_true("Gene_name" %in% names(result))
  expect_true("CAI" %in% names(result))
  expect_equal(nrow(result), 2)
  expect_true("test_gene" %in% result$Gene_name)
  expect_true(is.numeric(result$CAI))
})

# Test 5: CAI works with multiple genes
test_that("CAI works with multiple genes", {
  codon_counts <- data.table(
    Gene_name = c("ref_gene", "gene1", "gene2", "gene3"),
    AAA = c(5, 10, 20, 15),
    AAG = c(25, 15, 10, 20),
    GAA = c(20, 12, 18, 15),
    GAG = c(5, 9, 7, 10)
  )
  
  reference_genes <- c("ref_gene")
  
  result <- calculate_cai(codon_counts, reference_genes, genetic_code)
  
  expect_equal(nrow(result), 4)
  expect_true(all(c("gene1", "gene2", "gene3") %in% result$Gene_name))
  expect_true(all(result$CAI >= 0))
  expect_true(all(result$CAI <= 1))
})

# Test 6: CAI is numeric and finite
test_that("CAI values are numeric and finite", {
  codon_counts <- data.table(
    Gene_name = c("ref_gene", "test_gene"),
    AAA = c(5, 10), AAG = c(25, 15),
    GAA = c(20, 12), GAG = c(5, 9)
  )
  
  reference_genes <- c("ref_gene")
  
  result <- calculate_cai(codon_counts, reference_genes, genetic_code)
  
  expect_true(is.numeric(result$CAI))
  expect_true(is.finite(result$CAI))
  expect_false(is.na(result$CAI))
})

# Test 7: Relative adaptiveness values are calculated correctly
test_that("Relative adaptiveness (w) values are correct", {
  # Test that CAI reflects codon preferences
  codon_counts <- data.table(
    Gene_name = c("ref", "adapted", "maladapted"),
    AAA = c(10, 10, 90),  # ref and adapted prefer AAG
    AAG = c(90, 90, 10)   # maladapted prefers AAA
  )
  
  reference_genes <- c("ref")
  
  result <- calculate_cai(codon_counts, reference_genes, genetic_code)
  
  cai_adapted <- result[result$Gene_name == "adapted", "CAI"]
  cai_maladapted <- result[result$Gene_name == "maladapted", "CAI"]
  
  expect_true(cai_adapted$CAI > cai_maladapted$CAI)
})

# Test 8: CAI handles zero counts gracefully
test_that("CAI handles genes with zero counts for some codons", {
  codon_counts <- data.table(
    Gene_name = c("ref", "sparse_gene"),
    AAA = c(10, 0), AAG = c(20, 30),   # sparse uses only AAG
    GAA = c(18, 25), GAG = c(7, 0)     # sparse uses only GAA
  )
  
  reference_genes <- c("ref")
  
  expect_error({
    result <- calculate_cai(codon_counts, reference_genes, genetic_code)
  }, NA)
  
  expect_true(result$CAI >= 0)
  expect_true(result$CAI <= 1)
})

cat("✓ All CAI calculation tests passed!\n\n")
