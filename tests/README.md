# Unit Tests for CUB Analysis

This directory contains comprehensive unit tests for the Codon Usage Bias (CUB) analysis pipeline.

## Test Files

### Core Calculation Tests
- **`test_calculate_enc.R`**: Tests for Effective Number of Codons (ENC) calculation
  - Validates ENC stays within theoretical bounds [20, 61]
  - Tests extreme bias cases (one codon per AA)
  - Tests uniform usage (no bias)
  - Validates handling of missing codons

- **`test_calculate_rscu.R`**: Tests for Relative Synonymous Codon Usage (RSCU)
  - Validates mean RSCU = 1.0 within amino acid families
  - Tests uniform usage (RSCU = 1.0 for all codons)
  - Tests extreme bias scenarios
  - Validates 2-codon and multi-codon families

- **`test_calculate_gc_content.R`**: Tests for GC content calculations
  - Validates GC content bounds [0, 1]
  - Tests GC3s calculation (excludes Met and Trp)
  - Validates 4-fold degenerate sites
  - Tests edge cases (100% GC, 100% AT)

- **`test_tRNA_correlation.R`**: Tests for tRNA-codon correlation analysis
  - Validates anticodon to codon conversion (reverse complement)
  - Tests correlation with different methods (Spearman, Pearson, Kendall)
  - Validates RSCU calculation consistency
  - Tests handling of missing tRNA data

## Running the Tests

### Run All Tests
```r
# From the tests directory
source("run_all_tests.R")
```

### Run Individual Test Files
```r
# From the tests directory
library(testthat)
test_file("test_calculate_enc.R")
test_file("test_calculate_rscu.R")
test_file("test_calculate_gc_content.R")
test_file("test_tRNA_correlation.R")
```

### Run from Command Line
```bash
cd tests
Rscript run_all_tests.R
```

## Test Coverage

The test suite covers:
- ✅ **Mathematical correctness** of all CUB metrics
- ✅ **Edge cases** (extreme bias, uniform usage, zero counts)
- ✅ **Boundary conditions** (valid ranges for all metrics)
- ✅ **Multiple genes** (batch processing)
- ✅ **Data structure** (correct output format)
- ✅ **Error handling** (missing data, invalid inputs)

## Expected Output

When all tests pass, you should see:
```
========================================
  CUB Analysis - Unit Test Suite
========================================

Running ENC calculation tests...
✓ All ENC calculation tests passed!

Running RSCU calculation tests...
✓ All RSCU calculation tests passed!

Running GC content calculation tests...
✓ All GC content calculation tests passed!

Running tRNA-codon correlation tests...
✓ All tRNA-codon correlation tests passed!

========================================
  Test Summary
========================================

ENC tests:                     ✓ PASS (6/6 passed)
RSCU tests:                    ✓ PASS (8/8 passed)
GC tests:                      ✓ PASS (10/10 passed)
TRNA tests:                    ✓ PASS (6/6 passed)

Total: 30 tests, 30 passed, 0 failed

========================================
  ✓ ALL TESTS PASSED!
========================================
```

## Adding New Tests

To add new tests:

1. Create a new test file: `test_<function_name>.R`
2. Follow the structure:
```r
library(testthat)
library(data.table)
source("../src/<function_name>.R")

context("<Function Name> Tests")

test_that("<description of test>", {
  # Test code
  expect_equal(result, expected_value)
})
```
3. Add the test file to `run_all_tests.R`

## Dependencies

- `testthat`: R package for unit testing
- `data.table`: For efficient data manipulation

Install with:
```r
install.packages("testthat")
install.packages("data.table")
```

## Continuous Integration

These tests can be integrated into a CI/CD pipeline:
```bash
# Run tests and exit with error code if any fail
Rscript tests/run_all_tests.R || exit 1
```

## Test Philosophy

These tests follow best practices:
- **Independent**: Each test can run standalone
- **Reproducible**: Same inputs always produce same outputs
- **Fast**: All tests complete in seconds
- **Comprehensive**: Cover normal cases, edge cases, and error conditions
- **Documented**: Clear descriptions of what each test validates

## Troubleshooting

If tests fail:
1. Check that all source files are in the `../src/` directory
2. Ensure all required packages are installed
3. Verify the working directory is set correctly
4. Review the error messages for specific failures

## Contact

For questions about the test suite, refer to the main README or contact the repository maintainer.
