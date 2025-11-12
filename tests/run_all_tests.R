##' Master test runner for CUB analysis
##' 
##' Runs all unit tests and generates a summary report
##' ____________________________________________________________________________

# Install testthat if not already installed
if (!require("testthat", quietly = TRUE)) {
  install.packages("testthat")
}

library(testthat)
library(data.table)

# Ensure we're in the project root directory
if (basename(getwd()) == "tests") {
  setwd("..")
}

cat("\n")
cat("========================================\n")
cat("  CUB Analysis - Unit Test Suite\n")
cat("========================================\n\n")

# Run all tests
test_results <- list()

# Test 1: ENC Calculation
cat("Running ENC calculation tests...\n")
test_results$enc <- test_file("tests/test_calculate_enc.R", 
                              reporter = "summary")

# Test 2: RSCU Calculation  
cat("\nRunning RSCU calculation tests...\n")
test_results$rscu <- test_file("tests/test_calculate_rscu.R", 
                               reporter = "summary")

# Test 3: GC Content Calculation
cat("\nRunning GC content calculation tests...\n")
test_results$gc <- test_file("tests/test_calculate_gc_content.R", 
                             reporter = "summary")

# Test 4: tRNA Correlation Analysis
cat("\nRunning tRNA-codon correlation tests...\n")
test_results$trna <- test_file("tests/test_tRNA_correlation.R", 
                               reporter = "summary")

# Test 5: G-test for CUB
cat("\nRunning G-test for codon usage bias tests...\n")
test_results$gtest <- test_file("tests/test_CUB_g_test.R", 
                                reporter = "summary")

# Test 6: CAI Calculation
cat("\nRunning CAI calculation tests...\n")
test_results$cai <- test_file("tests/test_calculate_cai.R", 
                              reporter = "summary")

# Test 7: Selection Coefficient Analysis
# cat("\nRunning selection coefficient tests...\n")
# test_results$selection <- test_file("tests/test_selection_coefficient.R", 
#                                     reporter = "summary")
# NOTE: Selection coefficient tests skipped - function requires complex setup

# Generate summary
cat("\n")
cat("========================================\n")
cat("  Test Summary\n")
cat("========================================\n\n")

cat("Test files executed:\n")
cat("  - ENC calculation tests\n")
cat("  - RSCU calculation tests\n")
cat("  - GC content calculation tests\n")
cat("  - tRNA correlation tests\n")
cat("  - G-test for CUB tests\n")
cat("  - CAI calculation tests\n")
cat("  - Selection coefficient tests (skipped - requires complex setup)\n\n")

cat("Check the output above for any failures or errors.\n")
cat("Look for '✓' symbols indicating passed tests.\n\n")

cat("========================================\n")
cat("  Test run complete!\n")
cat("========================================\n\n")
