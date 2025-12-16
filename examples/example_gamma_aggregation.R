#' Example: Aggregating Gamma to Gene Level and Comparing with AnaCoDa
#' 
#' This script demonstrates how to:
#' 1. Aggregate per-AA gamma estimates to gene-level selection intensity
#' 2. Compare with AnaCoDa-derived selection coefficients
#' 3. Validate comparability of the two approaches
#' 
#' @author Luis Javier Madrigal-Roca
#' @date 2024-12-15

# Load required functions
source("./src/integrate_intronic_polymorphism.R")
source("./src/derivation_gamma_from_polymorphism.R")

# ==============================================================================
# STEP 1: Load gamma results (from polymorphism analysis)
# ==============================================================================

# Assuming you've already run estimate_gamma_by_gene_with_neutral_params()
# and have gamma_results object
#
# gamma_results should have columns:
#   - Gene
#   - AA
#   - Gamma (selection coefficient)
#   - Total_Alleles (sample size)
#   - Significant (logical)

# Example: Load from saved file
# gamma_results <- fread("./results/gamma_estimates_by_gene_AA.csv")

# ==============================================================================
# STEP 2: Aggregate gamma to gene level
# ==============================================================================

# This function creates gene-level metrics comparable to AnaCoDa
gamma_gene_level <- aggregate_gamma_per_gene(
  gamma_results = gamma_results,
  codon_usage_df = codon_usage,  # Your codon usage matrix
  genetic_code = genetic_code_dna_long  # Codon to AA mapping
)

# Output columns:
#   - Gene_name
#   - Selection_Intensity: mean(|gamma|) - DIRECTLY COMPARABLE to AnaCoDa S_coeff
#   - Gamma_Weighted_Mean: weighted mean of gamma values
#   - Gamma_Mean: simple mean of gamma values
#   - Total_Selection_Load: sum(|gamma|)
#   - N_AA_Analyzed: number of amino acids with estimates
#   - Prop_Significant: proportion of AAs with significant selection

# ==============================================================================
# STEP 3: Compare with AnaCoDa selection coefficients
# ==============================================================================

# Load AnaCoDa selection intensity (from main.R analysis)
# This is the output from section 8.3 where you calculate:
#   selection_coeff_intensity <- data.frame(
#     Gene_name = names(sel_intensity),
#     S_coeff = as.vector(sel_intensity)
#   )

# Direct comparison
comparison <- compare_gamma_with_anacoda(
  gamma_gene_level = gamma_gene_level,
  anacoda_intensity = selection_coeff_intensity
)

# This will print:
#   - Spearman correlation between Selection_Intensity and S_coeff
#   - Interpretation (strong/moderate/weak agreement)
#   - Quartile distributions for both metrics

# ==============================================================================
# STEP 4: Visualize comparison
# ==============================================================================

library(ggplot2)

# Scatter plot: Gamma vs AnaCoDa
p1 <- ggplot(comparison, aes(x = S_coeff, y = Selection_Intensity)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(
    title = "Comparison: Gamma vs AnaCoDa Selection Coefficients",
    subtitle = sprintf("Spearman ρ = %.3f (n = %d genes)",
                       cor(comparison$S_coeff, comparison$Selection_Intensity, 
                           method = "spearman"),
                       nrow(comparison)),
    x = "AnaCoDa Selection Intensity (S_coeff)",
    y = "Gamma Selection Intensity (mean |γ|)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11)
  )

# Log-log plot for better dynamic range
p2 <- ggplot(comparison, aes(x = log10(S_coeff + 0.01), 
                              y = log10(Selection_Intensity + 0.01))) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  labs(
    title = "Log-Scale Comparison",
    x = "log10(AnaCoDa S_coeff + 0.01)",
    y = "log10(Gamma Intensity + 0.01)"
  ) +
  theme_bw()

# Combine plots
library(patchwork)
combined_plot <- p1 / p2

ggsave("./results/gamma_vs_anacoda_comparison.pdf", 
       plot = combined_plot, 
       width = 10, height = 12)

# ==============================================================================
# STEP 5: Validate with CAI/CDC
# ==============================================================================

# Both gamma and AnaCoDa should correlate with CAI/CDC
validation <- validate_against_cai(
  gamma_results = gamma_results,
  integrated_data = integrated_data  # Your main data with CAI, CDC
)

# ==============================================================================
# INTERPRETATION GUIDE
# ==============================================================================

cat("\n=== Interpretation Guide ===\n\n")

cat("1. EXPECTED RESULTS:\n")
cat("   - Positive correlation between Gamma and AnaCoDa (ρ > 0.3)\n")
cat("   - Both correlate positively with CAI and CDC\n")
cat("   - Both correlate positively with expression\n\n")

cat("2. WHY MAGNITUDES DIFFER:\n")
cat("   - AnaCoDa: Based on codon frequencies within genes\n")
cat("   - Gamma: Based on population polymorphism data\n")
cat("   - Different scales (deltaEta vs 4Nes)\n")
cat("   - Different units (per codon vs per AA)\n\n")

cat("3. WHAT MATTERS:\n")
cat("   - RANK ORDER should be consistent\n")
cat("   - Genes with high AnaCoDa S_coeff should have high Gamma\n")
cat("   - Both identify the same biological signal (selection for codon bias)\n\n")

cat("4. BIOLOGICAL VALIDATION:\n")
cat("   - High-expression genes should have higher selection intensity\n")
cat("   - Ribosomal proteins should show strong positive selection\n")
cat("   - Genes with low CDC should have low/negative gamma\n\n")

# ==============================================================================
# EXAMPLE: Check ribosomal proteins
# ==============================================================================

# Identify ribosomal proteins (example gene set)
ribosomal_genes <- comparison$Gene_name[grep("^Rp", comparison$Gene_name)]

if (length(ribosomal_genes) > 0) {
  cat(sprintf("\n=== Ribosomal Proteins (n = %d) ===\n", length(ribosomal_genes)))
  
  ribo_data <- comparison[comparison$Gene_name %in% ribosomal_genes, ]
  
  cat(sprintf("Mean Gamma Intensity: %.3f\n", mean(ribo_data$Selection_Intensity)))
  cat(sprintf("Mean AnaCoDa S_coeff: %.3f\n", mean(ribo_data$S_coeff)))
  cat(sprintf("Proportion with Gamma > 1: %.1f%%\n\n", 
              100 * mean(ribo_data$Gamma_Mean > 1)))
}

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

# Save aggregated gamma data
fwrite(gamma_gene_level, "./results/gamma_gene_level_selection.csv")

# Save comparison table
fwrite(comparison, "./results/gamma_vs_anacoda_comparison.csv")

cat("\n✓ Analysis complete! Results saved to ./results/\n\n")
