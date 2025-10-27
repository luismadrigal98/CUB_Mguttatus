# Test Suite for calculate_gc_content.R
#
# Tests GC content calculation including overall GC, GC1, GC2, GC3, GC12, and GC3s
# GC3s excludes Met (ATG) and Trp (TGG) as they are not degenerate at 3rd position

library(testthat)
library(data.table)

# Source the function
source("../src/calculate_gc_content.R")

context("GC Content Calculation Tests")

# Test 1: GC content is properly bounded [0, 1]
test_that("GC content is between 0 and 1", {
  codon_counts <- data.table(
    Gene_name = "test_gene",
    AAA = 10, AAC = 5, AAG = 15, AAT = 8
  )
  
  gc_result <- calculate_gc_content(codon_counts)
  
  expect_true(gc_result$GC >= 0,
              info = "GC content should be >= 0.0")
  expect_true(gc_result$GC <= 1,
              info = "GC content should be <= 1.0")
  expect_true(gc_result$GC3s >= 0)
  expect_true(gc_result$GC3s <= 1)
})

# Test 2: 100% GC codons
test_that("100% GC codons gives GC content = 1.0", {
  # GGG = Gly (all G)
  codon_counts <- data.table(
    Gene_name = "all_gc_gene",
    GGG = 100, GGC = 50, GCG = 30, CGC = 40
  )
  
  gc_result <- calculate_gc_content(codon_counts)
  expect_equal(gc_result$GC, 1, tolerance = 0.01)
})

# Test 3: 100% AT codons
test_that("100% AT codons gives GC content = 0.0", {
  codon_counts <- data.table(
    Gene_name = "all_at_gene",
    AAA = 50, AAT = 40, TAT = 30, TAA = 10  # all A/T
  )
  
  gc_result <- calculate_gc_content(codon_counts)
  expect_equal(gc_result$GC, 0, tolerance = 0.01)
})

# Test 4: GC3s excludes Met and Trp
test_that("GC3s correctly excludes Met (ATG) and Trp (TGG)", {
  # Mix ATG/TGG with 4-fold degenerate sites
  codon_counts <- data.table(
    Gene_name = "test_gene",
    ATG = 10,  # Met - excluded from GC3s
    TGG = 10,  # Trp - excluded from GC3s
    # Alanine (4-fold degenerate) - 50% GC at position 3
    GCA = 5, GCC = 5, GCG = 5, GCT = 5
  )
  
  gc_result <- calculate_gc_content(codon_counts)
  
  # GC3s should only consider the Ala codons (half GC at pos 3)
  expect_equal(gc_result$GC3s, 0.5, tolerance = 0.01,
               info = "GC3s should exclude Met and Trp")
})

# Test 5: Correct calculation for known sequence
test_that("GC content calculated correctly for known sequence", {
  # Create a balanced codon set with 50% GC
  codon_counts <- data.table(
    Gene_name = "balanced_gene",
    GCA = 10,  # G,C,A = 2/3 GC
    TAT = 10,  # T,A,T = 0/3 GC
    CGC = 10   # C,G,C = 3/3 GC
  )
  
  gc_result <- calculate_gc_content(codon_counts)
  
  # Total: 10 codons * 3 bases = 30 bases
  # GC bases: GCA(2) + TAT(0) + CGC(3) = 5*10 = 50 bases out of 90
  # Expected: 50/90 = 0.556
  expect_equal(gc_result$GC, 5/9, tolerance = 0.01)
})

# Test 6: GC1, GC2, GC3 are calculated correctly
test_that("Position-specific GC content is correct", {
  codon_counts <- data.table(
    Gene_name = "test_gene",
    GCA = 10  # G(pos1)=GC, C(pos2)=GC, A(pos3)=AT
  )
  
  gc_result <- calculate_gc_content(codon_counts)
  
  expect_equal(gc_result$GC1, 1.0, tolerance = 0.01,
               info = "First position should be 100% GC")
  expect_equal(gc_result$GC2, 1.0, tolerance = 0.01,
               info = "Second position should be 100% GC")
  expect_equal(gc_result$GC3, 0.0, tolerance = 0.01,
               info = "Third position should be 0% GC")
})

# Test 7: GC12 is correctly calculated
test_that("GC12 combines first two positions correctly", {
  codon_counts <- data.table(
    Gene_name = "test_gene",
    GCA = 10,  # GC, GC, AT
    ATG = 10   # AT, AT, GC
  )
  
  gc_result <- calculate_gc_content(codon_counts)
  
  # Position 1+2: GCA contributes 2 GC / 2 total, ATG contributes 0 GC / 2 total
  # Total: 20 GC / 40 bases = 0.5
  expect_equal(gc_result$GC12, 0.5, tolerance = 0.01)
})

# Test 8: Output structure
test_that("GC content output has correct structure", {
  codon_counts <- data.table(
    Gene_name = "test_gene",
    AAA = 10
  )
  
  gc_result <- calculate_gc_content(codon_counts)
  
  expect_true("GC" %in% names(gc_result))
  expect_true("GC1" %in% names(gc_result))
  expect_true("GC2" %in% names(gc_result))
  expect_true("GC3" %in% names(gc_result))
  expect_true("GC12" %in% names(gc_result))
  expect_true("GC3s" %in% names(gc_result))
  expect_equal(nrow(gc_result), 1)
})

# Test 9: Handles single codon type
test_that("GC content handles single codon type", {
  codon_counts <- data.table(
    Gene_name = "single_codon",
    GGG = 100  # All G
  )
  
  gc_result <- calculate_gc_content(codon_counts)
  expect_equal(gc_result$GC, 1.0, tolerance = 0.01)
})

# Test 10: All GC metrics are consistent
test_that("GC content and GC3s are internally consistent", {
  codon_counts <- data.table(
    Gene_name = "test_gene",
    AAA = 10, GGG = 10, CCC = 10, TTT = 10
  )
  
  gc_result <- calculate_gc_content(codon_counts)
  
  # All metrics should be bounded [0, 1]
  expect_true(gc_result$GC >= 0 && gc_result$GC <= 1)
  expect_true(gc_result$GC1 >= 0 && gc_result$GC1 <= 1)
  expect_true(gc_result$GC2 >= 0 && gc_result$GC2 <= 1)
  expect_true(gc_result$GC3 >= 0 && gc_result$GC3 <= 1)
  expect_true(gc_result$GC12 >= 0 && gc_result$GC12 <= 1)
  expect_true(gc_result$GC3s >= 0 && gc_result$GC3s <= 1)
  
  # All should be numeric
  expect_true(is.numeric(gc_result$GC))
  expect_true(is.numeric(gc_result$GC3s))
})

cat("✓ All GC content calculation tests passed!\n\n")
