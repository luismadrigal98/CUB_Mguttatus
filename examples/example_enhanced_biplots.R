#!/usr/bin/env Rscript

##' Example usage of enhanced biplot functions
##' 
##' This script demonstrates how to use the new statistical testing
##' and enhanced biplot functions independently of main.R
##' 
##' @author Luis J. Madrigal-Roca
##' @date November 12, 2025

# Setup ----
library(data.table)
library(dplyr)
library(ggplot2)

# Source functions
source("./src/test_codon_proportions.R")
source("./src/enhanced_biplot.R")

# Load required data ----

cat("Loading data...\n")

# 1. Load genetic code
genetic_code_dna <- c(
  TTT = "Phe", TTC = "Phe", TTA = "Leu", TTG = "Leu",
  TCT = "Ser", TCC = "Ser", TCA = "Ser", TCG = "Ser",
  TAT = "Tyr", TAC = "Tyr", TAA = "STOP", TAG = "STOP",
  TGT = "Cys", TGC = "Cys", TGA = "STOP", TGG = "Trp",
  CTT = "Leu", CTC = "Leu", CTA = "Leu", CTG = "Leu",
  CCT = "Pro", CCC = "Pro", CCA = "Pro", CCG = "Pro",
  CAT = "His", CAC = "His", CAA = "Gln", CAG = "Gln",
  CGT = "Arg", CGC = "Arg", CGA = "Arg", CGG = "Arg",
  ATT = "Ile", ATC = "Ile", ATA = "Ile", ATG = "Met",
  ACT = "Thr", ACC = "Thr", ACA = "Thr", ACG = "Thr",
  AAT = "Asn", AAC = "Asn", AAA = "Lys", AAG = "Lys",
  AGT = "Ser", AGC = "Ser", AGA = "Arg", AGG = "Arg",
  GTT = "Val", GTC = "Val", GTA = "Val", GTG = "Val",
  GCT = "Ala", GCC = "Ala", GCA = "Ala", GCG = "Ala",
  GAT = "Asp", GAC = "Asp", GAA = "Glu", GAG = "Glu",
  GGT = "Gly", GGC = "Gly", GGA = "Gly", GGG = "Gly"
)

# 2. Load codon usage data (assumes main.R has been run)
if (!file.exists("./results/cub_analysis_complete.csv")) {
  stop("Error: Run main.R first to generate required data files")
}

cub_data <- read.csv("./results/cub_analysis_complete.csv")

# 3. Load expression data and create groups
if (file.exists("./results/expression_enc_cai.csv")) {
  expression_data <- read.csv("./results/expression_enc_cai.csv")
  
  # Create expression groups (Top 5%, Bottom 5%, Middle 90%)
  top5_threshold <- quantile(expression_data$log10_mean_cpm, 0.95, na.rm = TRUE)
  bottom5_threshold <- quantile(expression_data$log10_mean_cpm, 0.05, na.rm = TRUE)
  
  expression_data <- expression_data %>%
    mutate(
      Expression_Group = case_when(
        log10_mean_cpm >= top5_threshold ~ "Top 5%",
        log10_mean_cpm <= bottom5_threshold ~ "Bottom 5%",
        TRUE ~ "Middle 90%"
      )
    )
} else {
  stop("Error: Expression data not found. Run main.R through section 5.")
}

# 4. Load codon counts (from codon_usage object in main.R)
# For this example, we'll simulate - in practice, load from main.R output
cat("Note: This example requires running main.R first to generate all data\n")
cat("See main.R sections 6.3.1 and 7 for full implementation\n\n")

# Example: Test codon proportions ----

cat("=== Example 1: Testing Codon Proportions ===\n\n")

# This assumes you have codon_usage data.table from main.R
# For demonstration purposes only - replace with actual data

# Load CAI weights if available
if (file.exists("./results/optimal_codons_weights.csv")) {
  w_table <- read.csv("./results/optimal_codons_weights.csv")
  cat("✓ Loaded CAI weights\n")
} else {
  cat("⚠ CAI weights not found - run main.R through section 6.2\n")
}

# Example workflow (requires actual codon usage data):
cat("\nTo run statistical tests:\n")
cat("1. Prepare codon usage data for Top 5% and Bottom 5% genes\n")
cat("2. Call test_codon_proportions():\n\n")

cat("  test_results <- test_codon_proportions(\n")
cat("    selected_usage = top5_codon_usage,\n")
cat("    neutral_usage = bottom95_codon_usage,\n")
cat("    genetic_code = genetic_code_dna,\n")
cat("    method = 'chisq',\n")
cat("    fdr_correction = TRUE\n")
cat("  )\n\n")

cat("3. Create summary plots:\n\n")
cat("  plot_codon_selection_summary(test_results, \n")
cat("    output_file = 'codon_selection.pdf')\n\n")

cat("  plot_codon_classification_heatmap(test_results, w_table,\n")
cat("    output_file = 'codon_classification.pdf')\n\n")

# Example: Enhanced biplots ----

cat("\n=== Example 2: Creating Enhanced Biplots ===\n\n")

cat("To create enhanced biplots:\n")
cat("1. Run CA or PCA analysis (see main.R section 7)\n")
cat("2. Prepare gene data with expression groups\n")
cat("3. Call create_enhanced_biplot():\n\n")

cat("  # For CA results\n")
cat("  enhanced_ca <- create_enhanced_biplot(\n")
cat("    ordination_result = ca_result,\n")
cat("    gene_data = gene_expr_groups,\n")
cat("    codon_test_results = test_results,\n")
cat("    w_table = w_table,\n")
cat("    dims = c(1, 2),\n")
cat("    color_by = 'combined',\n")
cat("    show_only_significant = FALSE,\n")
cat("    title = 'Enhanced CA Biplot',\n")
cat("    output_file = 'enhanced_ca_biplot.pdf'\n")
cat("  )\n\n")

cat("  # For PCA results\n")
cat("  enhanced_pca <- create_enhanced_biplot(\n")
cat("    ordination_result = pca_result,\n")
cat("    gene_data = gene_expr_groups,\n")
cat("    codon_test_results = test_results,\n")
cat("    w_table = w_table,\n")
cat("    dims = c(1, 2),\n")
cat("    color_by = 'selection',\n")
cat("    show_only_significant = TRUE,\n")
cat("    title = 'Enhanced PCA Biplot',\n")
cat("    output_file = 'enhanced_pca_biplot.pdf'\n")
cat("  )\n\n")

# Example: Loading analysis ----

cat("\n=== Example 3: Analyzing Codon Loading Patterns ===\n\n")

cat("To analyze which codons drive separation:\n\n")

cat("  loading_analysis <- analyze_codon_loading_patterns(\n")
cat("    ordination_result = ca_result,\n")
cat("    codon_test_results = test_results,\n")
cat("    dims = c(1, 2)\n")
cat("  )\n\n")

cat("  # Save results\n")
cat("  write.csv(loading_analysis, 'loading_analysis.csv')\n\n")

cat("  # Check if significant codons have higher loadings\n")
cat("  sig_loadings <- loading_analysis %>% filter(Significant)\n")
cat("  cat(sprintf('Mean loading magnitude (sig): %.3f\\n',\n")
cat("              mean(sig_loadings$Loading_Magnitude)))\n\n")

# Color scheme reference ----

cat("\n=== Color Scheme Reference ===\n\n")

cat("Selection Status (color_by = 'selection'):\n")
cat("  🔴 Red    : Under Selection + Preferred (w=1)\n")
cat("  🟠 Orange : Under Selection (not preferred)\n")
cat("  🔵 Blue   : Avoided in high expression\n")
cat("  ⚪ Gray   : Neutral (no significant difference)\n\n")

cat("CAI Preference (color_by = 'preference'):\n")
cat("  🔴 Red  : Preferred codons (w = 1.0)\n")
cat("  ⚪ Gray : Non-preferred codons (w < 1.0)\n\n")

cat("AT vs GC Ending (color_by = 'ending'):\n")
cat("  🟠 Orange : AT-ending codons\n")
cat("  🔵 Blue   : GC-ending codons\n\n")

cat("Combined (color_by = 'combined'):\n")
cat("  🔴 Dark Red  : Selection + Preferred + GC-ending\n")
cat("  🟠 Orange    : Selection + Preferred + AT-ending\n")
cat("  🟡 Yellow    : Selection + Non-pref + GC-ending\n")
cat("  🟨 Lt Yellow : Selection + Non-pref + AT-ending\n")
cat("  🔵 Blue      : Avoided in high expression\n")
cat("  ⚪ Gray      : Neutral\n\n")

# Summary ----

cat("=== Quick Start ===\n\n")

cat("1. Run main.R through section 7 to generate all required data\n")
cat("2. Check results/codon_proportion_test_results.csv for test results\n")
cat("3. View enhanced biplots in results/ directory:\n")
cat("   - CA_enhanced_biplot_*.pdf\n")
cat("   - PCA_enhanced_biplot_*.pdf\n")
cat("   - *_biplot_panel_4plots.pdf (multi-panel figures)\n")
cat("4. Read loading analysis CSV files for quantitative results\n")
cat("5. Consult doc/ENHANCED_BIPLOT_GUIDE.md for interpretation\n\n")

cat("✓ Example usage demonstration complete!\n")
cat("See main.R for full working implementation.\n")
