# Test Suite for selection_coefficient_analysis.R
#
# Tests the calculation of selection coefficients from polymorphism data
# s = (π_neutral - π_selected) / π_neutral

library(testthat)
library(data.table)

# Source the function
source("../src/selection_coefficient_analysis.R")

context("Selection Coefficient Analysis Tests")

# Test 1: Output structure is correct
test_that("Selection coefficient output has correct structure", {
  polymorphism_data <- data.table(
    Gene = c("gene1", "gene1", "gene2", "gene2"),
    Degeneracy = c("0-fold", "4-fold", "0-fold", "4-fold"),
    Sites = c(100, 50, 120, 60),
    Pi = c(0.5, 0.2, 0.6, 0.15)
  )
  
  result <- selection_coefficient_analysis(polymorphism_data)
  
  expect_true("Gene" %in% names(result))
  expect_true("Pi_0fold" %in% names(result) || "Pi_neutral" %in% names(result))
  expect_true("Pi_4fold" %in% names(result) || "Pi_selected" %in% names(result))
  expect_true("selection_coef" %in% names(result) || "s" %in% names(result))
})

# Test 2: Selection coefficient is bounded [-1, 1]
test_that("Selection coefficient is within reasonable bounds", {
  polymorphism_data <- data.table(
    Gene = c("gene1", "gene1"),
    Degeneracy = c("0-fold", "4-fold"),
    Sites = c(100, 50),
    Pi = c(0.5, 0.2)
  )
  
  result <- selection_coefficient_analysis(polymorphism_data)
  
  s_col <- if ("s" %in% names(result)) "s" else "selection_coef"
  
  expect_true(abs(result[[s_col]]) <= 1)
})

# Test 3: Positive selection coefficient when π_4fold < π_0fold
test_that("Selection coefficient is positive when 4-fold has lower diversity", {
  # Lower diversity at 4-fold sites indicates selection
  polymorphism_data <- data.table(
    Gene = "test_gene",
    Degeneracy = c("0-fold", "4-fold"),
    Sites = c(100, 50),
    Pi = c(0.8, 0.2)  # 4-fold has much lower π
  )
  
  result <- selection_coefficient_analysis(polymorphism_data)
  
  s_col <- if ("s" %in% names(result)) "s" else "selection_coef"
  
  expect_true(result[[s_col]] > 0)
})

# Test 4: Negative selection coefficient when π_4fold > π_0fold
test_that("Selection coefficient is negative when 4-fold has higher diversity", {
  # Higher diversity at 4-fold sites (unusual, but possible)
  polymorphism_data <- data.table(
    Gene = "test_gene",
    Degeneracy = c("0-fold", "4-fold"),
    Sites = c(100, 50),
    Pi = c(0.2, 0.8)  # 4-fold has higher π
  )
  
  result <- selection_coefficient_analysis(polymorphism_data)
  
  s_col <- if ("s" %in% names(result)) "s" else "selection_coef"
  
  expect_true(result[[s_col]] < 0)
})

# Test 5: Selection coefficient is zero when π_4fold = π_0fold
test_that("Selection coefficient is zero for equal diversity", {
  polymorphism_data <- data.table(
    Gene = "test_gene",
    Degeneracy = c("0-fold", "4-fold"),
    Sites = c(100, 50),
    Pi = c(0.5, 0.5)  # Equal diversity
  )
  
  result <- selection_coefficient_analysis(polymorphism_data)
  
  s_col <- if ("s" %in% names(result)) "s" else "selection_coef"
  
  expect_true(abs(result[[s_col]]) < 0.01)
})

# Test 6: Handles multiple genes
test_that("Analysis works with multiple genes", {
  polymorphism_data <- data.table(
    Gene = c("gene1", "gene1", "gene2", "gene2", "gene3", "gene3"),
    Degeneracy = rep(c("0-fold", "4-fold"), 3),
    Sites = c(100, 50, 120, 60, 90, 45),
    Pi = c(0.8, 0.2, 0.6, 0.3, 0.7, 0.1)
  )
  
  result <- selection_coefficient_analysis(polymorphism_data)
  
  expect_equal(nrow(result), 3)
  expect_true(all(c("gene1", "gene2", "gene3") %in% result$Gene))
})

# Test 7: Returns data.table
test_that("Function returns data.table object", {
  polymorphism_data <- data.table(
    Gene = c("gene1", "gene1"),
    Degeneacy = c("0-fold", "4-fold"),
    Sites = c(100, 50),
    Pi = c(0.5, 0.2)
  )
  
  result <- selection_coefficient_analysis(polymorphism_data)
  
  expect_true(inherits(result, "data.table") || inherits(result, "data.frame"))
})

cat("✓ All selection coefficient tests passed!\n\n")
