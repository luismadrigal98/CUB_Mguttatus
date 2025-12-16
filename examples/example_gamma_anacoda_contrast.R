#' Mathematical Contrast: Polymorphism vs AnaCoDa Selection Estimates
#' 
#' This script demonstrates the rigorous validation framework comparing
#' polymorphism-based gamma estimates with AnaCoDa's mechanistic model.
#' 
#' Key Innovation: Uses mathematically precise aggregation formula
#'   S_poly = (1/L) * Sum_AA(Count_Unpref_AA * gamma_AA)
#' 
#' @author Luis Javier Madrigal-Roca
#' @date 2024-12-15

# Load required functions
source("./src/integrate_intronic_polymorphism.R")
source("./src/derivation_gamma_from_polymorphism.R")

# ==============================================================================
# PREREQUISITE DATA
# ==============================================================================

# You should have already run:
# 1. estimate_gamma_by_gene_with_neutral_params() → gamma_results
# 2. AnaCoDa analysis → selection_coeff_intensity
# 3. Codon usage quantification → codon_usage
# 4. Preferred codon identification → preferred_codons

# Example loading (adjust paths as needed):
# gamma_results <- fread("./results/gamma_estimates_by_gene_AA.csv")
# codon_usage <- fread("./results/codon_usage_matrix.csv")
# preferred_codons <- fread("./data/preferred_codons_roc.csv")

# ==============================================================================
# MATHEMATICAL CONTRAST ANALYSIS
# ==============================================================================

cat("\n" %+% "=" %+% rep("=", 78) %+% "\n")
cat("VALIDATION FRAMEWORK: Polymorphism vs AnaCoDa\n")
cat(rep("=", 80) %+% "\n\n")

cat("This analysis implements the rigorous mathematical contrast:\n\n")
cat("  S_poly = (1/L) * Σ_AA [Count_Unpref(AA) * γ_AA]\n\n")
cat("Where:\n")
cat("  • L = Total gene length (codons)\n")
cat("  • Count_Unpref(AA) = Number of unpreferred codons for each AA\n")
cat("  • γ_AA = Selection coefficient (4Nes) favoring preferred codon\n\n")
cat("This metric is DIRECTLY comparable to AnaCoDa's S_coeff.\n")
cat(rep("-", 80) %+% "\n\n")

# Run the contrast analysis
contrast_results <- contrast_gamma_anacoda(
  gamma_results = gamma_results,
  codon_usage = codon_usage,
  preferred_codons = preferred_codons,
  anacoda_intensity = selection_coeff_intensity,
  genetic_code = genetic_code_dna_long  # From main.R
)

# ==============================================================================
# RESULTS INTERPRETATION
# ==============================================================================

cat("\n=== INTERPRETATION GUIDE ===\n\n")

cat("Expected Results for Valid Model:\n")
cat("  1. Spearman ρ > 0.5: Strong concordance\n")
cat("  2. p-value < 0.01: Highly significant\n")
cat("  3. Both methods rank genes similarly\n\n")

cat("What High Correlation Means:\n")
cat("  ✓ Polymorphism data validates mechanistic model\n")
cat("  ✓ Both approaches detect same biological signal\n")
cat("  ✓ Selection on codon bias is real, not artifact\n\n")

cat("What Low Correlation Suggests:\n")
cat("  • Different timescales (contemporary vs historical)\n")
cat("  • Model violations (demography, recombination)\n")
cat("  • Measurement error in one or both methods\n\n")

# ==============================================================================
# BIOLOGICAL VALIDATION
# ==============================================================================

cat("=== BIOLOGICAL VALIDATION ===\n\n")

# Check if high-expression genes show higher selection
if ("Geom_Exp" %in% names(contrast_results) || 
    exists("integrated_data")) {
  
  # Merge with expression if not already present
  if (!"Geom_Exp" %in% names(contrast_results) && exists("integrated_data")) {
    contrast_results <- merge(contrast_results, 
                              integrated_data[, .(Gene_name, Geom_Exp)],
                              by = "Gene_name", all.x = TRUE)
  }
  
  # Quartile analysis
  q75_exp <- quantile(contrast_results$Geom_Exp, 0.75, na.rm = TRUE)
  q25_exp <- quantile(contrast_results$Geom_Exp, 0.25, na.rm = TRUE)
  
  high_exp <- contrast_results[Geom_Exp >= q75_exp]
  low_exp <- contrast_results[Geom_Exp <= q25_exp]
  
  cat("High Expression (Q4) vs Low Expression (Q1):\n\n")
  cat(sprintf("  S_poly (High Exp): Mean = %.4f, SD = %.4f\n",
              mean(high_exp$S_poly, na.rm = TRUE),
              sd(high_exp$S_poly, na.rm = TRUE)))
  cat(sprintf("  S_poly (Low Exp):  Mean = %.4f, SD = %.4f\n\n",
              mean(low_exp$S_poly, na.rm = TRUE),
              sd(low_exp$S_poly, na.rm = TRUE)))
  
  # Wilcoxon test
  wilcox_poly <- wilcox.test(high_exp$S_poly, low_exp$S_poly)
  wilcox_anacoda <- wilcox.test(high_exp$S_coeff, low_exp$S_coeff)
  
  cat(sprintf("  Wilcoxon p-value (S_poly):   %.2e\n", wilcox_poly$p.value))
  cat(sprintf("  Wilcoxon p-value (S_coeff):  %.2e\n\n", wilcox_anacoda$p.value))
  
  if (wilcox_poly$p.value < 0.01 && wilcox_anacoda$p.value < 0.01) {
    cat("✓ BOTH methods detect expression-selection correlation!\n")
    cat("  Strong evidence for translational selection hypothesis.\n\n")
  }
}

# ==============================================================================
# DIAGNOSTIC PLOTS
# ==============================================================================

cat("=== GENERATING DIAGNOSTIC PLOTS ===\n\n")

library(ggplot2)
library(patchwork)

# Plot 1: Log-log scale for better dynamic range
p_log <- ggplot(contrast_results, 
                aes(x = log10(S_coeff + 0.001), 
                    y = log10(S_poly + 0.001))) +
  geom_point(alpha = 0.4, size = 1.5, color = "steelblue") +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(
    title = "Log-Scale Comparison",
    x = expression(log[10](S[AnaCoDa])),
    y = expression(log[10](S[poly]))
  ) +
  theme_bw()

# Plot 2: Residuals
contrast_results[, Residual := S_poly - S_coeff]

p_resid <- ggplot(contrast_results, 
                  aes(x = S_coeff, y = Residual)) +
  geom_point(alpha = 0.4, size = 1.5, color = "darkgreen") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", color = "blue", se = TRUE) +
  labs(
    title = "Residual Analysis",
    x = expression(S[AnaCoDa]),
    y = expression(S[poly] - S[AnaCoDa])
  ) +
  theme_bw()

# Combine and save
combined <- p_log / p_resid

ggsave("./results/gamma_anacoda_diagnostics.pdf",
       plot = combined, width = 8, height = 10)

cat("✓ Diagnostic plots saved: ./results/gamma_anacoda_diagnostics.pdf\n\n")

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

# Save contrast table
fwrite(contrast_results, "./results/gamma_anacoda_contrast.csv")

cat("✓ Results saved: ./results/gamma_anacoda_contrast.csv\n\n")

# ==============================================================================
# SUMMARY REPORT
# ==============================================================================

cat(rep("=", 80) %+% "\n")
cat("VALIDATION SUMMARY\n")
cat(rep("=", 80) %+% "\n\n")

cat(sprintf("Genes analyzed: %d\n", nrow(contrast_results)))
cat(sprintf("Spearman correlation: ρ = %.4f\n", 
            unique(contrast_results$Spearman_rho)[1]))
cat(sprintf("Statistical significance: p = %.2e\n\n", 
            unique(contrast_results$Spearman_p)[1]))

if (unique(contrast_results$Spearman_rho)[1] > 0.5) {
  cat("CONCLUSION: ✓ VALIDATION SUCCESSFUL\n")
  cat("  • Strong concordance between independent methods\n")
  cat("  • Polymorphism data supports mechanistic model\n")
  cat("  • Selection on codon bias is robustly detected\n\n")
} else if (unique(contrast_results$Spearman_rho)[1] > 0.3) {
  cat("CONCLUSION: ⚠ MODERATE VALIDATION\n")
  cat("  • Methods show general agreement\n")
  cat("  • Some discrepancies may warrant investigation\n")
  cat("  • Selection signal present but method-dependent\n\n")
} else {
  cat("CONCLUSION: ✗ VALIDATION INCONCLUSIVE\n")
  cat("  • Weak concordance requires further investigation\n")
  cat("  • Check for systematic biases or model violations\n")
  cat("  • Consider alternative explanations\n\n")
}

cat(rep("=", 80) %+% "\n")
cat("Analysis complete!\n")
cat(rep("=", 80) %+% "\n\n")
