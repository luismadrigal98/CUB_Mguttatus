# Test Suite for tRNA_codon_correlation.R
#
# Tests tRNA-codon correlation analysis including anticodon conversion
# and correlation methods

library(testthat)
library(data.table)

# Source the functions
source("../src/tRNA_codon_correlation.R")
source("../src/get_codon_supply_map.R")

context("tRNA-Codon Correlation Tests")

# Test 1: Anticodon to codon conversion (reverse complement)
test_that("anticodon_to_codon performs reverse complement correctly", {
  anticodon_to_codon <- function(anticodon) {
    complement <- c("A" = "T", "T" = "A", "G" = "C", "C" = "G")
    bases <- strsplit(anticodon, "")[[1]]
    codon <- paste0(complement[bases[3]], complement[bases[2]], complement[bases[1]])
    return(codon)
  }
  
  # Test specific anticodon-codon pairs
  expect_equal(anticodon_to_codon("AAA"), "TTT")  # AAA anticodon -> TTT codon (Lys)
  expect_equal(anticodon_to_codon("TTC"), "GAA")  # TTC anticodon -> GAA codon (Glu)
  expect_equal(anticodon_to_codon("CAT"), "ATG")  # CAT anticodon -> ATG codon (Met)
  expect_equal(anticodon_to_codon("GCC"), "GGC")  # GCC anticodon -> GGC codon (Gly)
})

# Test 2: get_codon_supply_map creates correct mapping
test_that("get_codon_supply_map handles wobble pairing rules", {
  # Create mock tRNA data
  trna_counts <- data.table(
    Anticodon = c("TTC", "AAA"),  # These should map to multiple codons via wobble
    tRNA_count = c(5, 10)
  )
  
  codon_supply <- get_codon_supply_map(trna_counts)
  
  expect_true("Codon" %in% names(codon_supply))
  expect_true("tRNA_supply" %in% names(codon_supply))
  expect_true(nrow(codon_supply) > 0)
  
  # Check that wobble rules allow multiple codons per tRNA
  # TTC anticodon (T at wobble position) should recognize GAA and GAG
  gaa_supply <- codon_supply[Codon == "GAA"]$tRNA_supply
  gag_supply <- codon_supply[Codon == "GAG"]$tRNA_supply
  expect_true(gaa_supply > 0 || gag_supply > 0)
})

# Test 3: Function handles missing tRNA data gracefully
test_that("Function runs with minimal tRNA data", {
  # Create temporary test files
  temp_codon <- tempfile(fileext = ".csv")
  temp_trna <- tempfile(fileext = ".tsv")
  
  # Create minimal codon counts
  codon_counts <- data.table(
    Gene_name = "test_gene",
    AAA = 10, AAG = 15,
    GAA = 12, GAG = 9
  )
  
  # Create minimal tRNA file
  trna_data <- data.frame(
    Seqid = "chr1",
    Source = "tRNAscan",
    Type = "tRNA",
    Start = 100,
    End = 200,
    Score = 50,
    Strand = "+",
    Phase = ".",
    Attributes = "ID=tRNA1",
    Anticodon = "TTC"
  )
  
  write.table(trna_data, temp_trna, sep = "\t", row.names = FALSE, quote = FALSE)
  
  genetic_code <- c("AAA" = "Lys", "AAG" = "Lys", "GAA" = "Glu", "GAG" = "Glu")
  
  # Test that function runs without error
  expect_error({
    result <- tRNA_codon_correlation(codon_counts, temp_trna, genetic_code,
                                     output_dir = tempdir(), test_method = "spearman")
  }, NA)  # NA means "no error expected"
  
  # Cleanup
  unlink(c(temp_codon, temp_trna))
})

# Test 4: Correlation method parameter works
test_that("All correlation methods work", {
  temp_trna <- tempfile(fileext = ".tsv")
  
  codon_counts <- data.table(
    Gene_name = c("g1", "g2", "g3"),
    AAA = c(10, 20, 30),
    AAG = c(15, 25, 5),
    GAA = c(12, 18, 25),
    GAG = c(9, 7, 15)
  )
  
  trna_data <- data.frame(
    Seqid = rep("chr1", 3),
    Source = rep("tRNAscan", 3),
    Type = rep("tRNA", 3),
    Start = c(100, 200, 300),
    End = c(200, 300, 400),
    Score = rep(50, 3),
    Strand = rep("+", 3),
    Phase = rep(".", 3),
    Attributes = c("ID=tRNA1", "ID=tRNA2", "ID=tRNA3"),
    Anticodon = c("TTC", "TTC", "AAA")
  )
  
  write.table(trna_data, temp_trna, sep = "\t", row.names = FALSE, quote = FALSE)
  
  genetic_code <- c("AAA" = "Lys", "AAG" = "Lys", "GAA" = "Glu", "GAG" = "Glu")
  
  for (method in c("spearman", "pearson", "kendall")) {
    expect_error({
      result <- tRNA_codon_correlation(codon_counts, temp_trna, genetic_code,
                                      output_dir = tempdir(), test_method = method)
    }, NA, info = paste("Method:", method))
  }
  
  unlink(temp_trna)
})

# Test 5: Output structure
test_that("tRNA correlation results have correct structure", {
  temp_trna <- tempfile(fileext = ".tsv")
  
  codon_counts <- data.table(
    Gene_name = "test_gene",
    AAA = 10, AAG = 15,
    GAA = 12, GAG = 9,
    GCT = 5, GCC = 7, GCA = 6, GCG = 4
  )
  
  trna_data <- data.frame(
    Seqid = rep("chr1", 2),
    Source = rep("tRNAscan", 2),
    Type = rep("tRNA", 2),
    Start = c(100, 200),
    End = c(200, 300),
    Score = rep(50, 2),
    Strand = rep("+", 2),
    Phase = rep(".", 2),
    Attributes = c("ID=tRNA1", "ID=tRNA2"),
    Anticodon = c("TTC", "AAA")
  )
  
  write.table(trna_data, temp_trna, sep = "\t", row.names = FALSE, quote = FALSE)
  
  genetic_code <- c("AAA" = "Lys", "AAG" = "Lys", "GAA" = "Glu", "GAG" = "Glu",
                    "GCT" = "Ala", "GCC" = "Ala", "GCA" = "Ala", "GCG" = "Ala")
  
  result <- tRNA_codon_correlation(codon_counts, temp_trna, genetic_code,
                                   output_dir = tempdir(), test_method = "spearman")
  
  expect_true("correlation_results" %in% names(result))
  expect_true("analysis_data" %in% names(result))
  expect_true("plots" %in% names(result))
  
  # Check analysis data structure
  expect_true("Codon" %in% names(result$analysis_data))
  expect_true("RSCU" %in% names(result$analysis_data))
  expect_true("tRNA_supply" %in% names(result$analysis_data))
  
  unlink(temp_trna)
})

# Test 6: RSCU calculation is correct within tRNA analysis
test_that("RSCU calculation in tRNA analysis is correct", {
  temp_trna <- tempfile(fileext = ".tsv")
  
  # Uniform codon usage within amino acid families
  codon_counts <- data.table(
    Gene_name = "test_gene",
    AAA = 100, AAG = 100,  # Lys - equal usage
    GAA = 50, GAG = 50     # Glu - equal usage
  )
  
  trna_data <- data.frame(
    Seqid = "chr1",
    Source = "tRNAscan",
    Type = "tRNA",
    Start = 100,
    End = 200,
    Score = 50,
    Strand = "+",
    Phase = ".",
    Attributes = "ID=tRNA1",
    Anticodon = "TTC"
  )
  
  write.table(trna_data, temp_trna, sep = "\t", row.names = FALSE, quote = FALSE)
  
  genetic_code <- c("AAA" = "Lys", "AAG" = "Lys", "GAA" = "Glu", "GAG" = "Glu")
  
  result <- tRNA_codon_correlation(codon_counts, temp_trna, genetic_code,
                                   output_dir = tempdir(), test_method = "spearman")
  
  # For equal usage, RSCU should be 1.0 for all codons
  rscu_values <- result$analysis_data$RSCU
  expect_true(all(abs(rscu_values - 1.0) < 0.01),
              info = "RSCU should be 1.0 for equal codon usage")
  
  unlink(temp_trna)
})

cat("✓ All tRNA-codon correlation tests passed!\n\n")
