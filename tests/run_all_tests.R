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

# Set working directory to project root
setwd("..")

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

# Generate summary
cat("\n")
cat("========================================\n")
cat("  Test Summary\n")
cat("========================================\n\n")

cat("Test files executed:\n")
cat("  - ENC calculation tests\n")
cat("  - RSCU calculation tests\n")
cat("  - GC content calculation tests\n")
cat("  - tRNA correlation tests\n\n")

cat("Check the output above for any failures or errors.\n")
cat("Look for '✓' symbols indicating passed tests.\n\n")

cat("========================================\n")
cat("  Test run complete!\n")
cat("========================================\n\n")
