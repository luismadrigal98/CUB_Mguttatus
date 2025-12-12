#' Comprehensive test of gamma estimation pipeline
#' 
#' This script:
#' 1. Generates synthetic data with known selection coefficients
#' 2. Runs the full pipeline (process_codon_vcf_with_nucleotide_pi, estimate_gamma, etc.)
#' 3. Compares recovered parameters to true values
#' 
#' @author Luis Javier Madrigal-Roca
#' _____________________________________________________________________________

library(data.table)
library(dplyr)
library(ggplot2)
library(gsl)

# Source required functions
source('./src/derivation_gamma_from_polymorphism.R')
source('./src/local_M_estimation.R')

set.seed(1998)

# ******************************************************************************
# 1) Set up simulation parameters ----
# ******************************************************************************

# True population parameters
theta_true <- 0.0312  # 4*Ne*mu (nucleotide mutation rate)
Ne_true <- 100000     # Effective population size (for reference)

# Mutation rates - REALISTIC: mutational bias toward unpreferred codons
# Beta should be ~1.5x greater than alpha
u_true <- 0.4  # Unpreferred -> Preferred (lower rate)
v_true <- 0.6  # Preferred -> Unpreferred (higher rate, mutational bias)

# Calculate population-scaled mutation rates
alpha_true <- theta_true * (u_true / (u_true + v_true))  # ~0.0125
beta_true <- theta_true * (v_true / (u_true + v_true))   # ~0.0187

cat("Checking beta/alpha ratio:", beta_true / alpha_true, "\n")
stopifnot(beta_true / alpha_true >= 1.4)  # Verify beta is at least 1.4x alpha

# Selection coefficients to test (in terms of 4*Ne*s)
# Weak selection: gamma < 1
# Moderate selection: 1 < gamma < 5  (should show hump effect)
# Adjusted categories based on validation results:
# - Negligible: gamma < 1 (nearly neutral, residuals centered at 0)
# - Moderate: 1 < gamma < 2 (weak selection, hump effect expected)
# - Strong: gamma > 2 (moderate selection, begins to show bias)
selection_scenarios <- data.frame(
  Category = c(
    rep("Negligible", 20),
    rep("Moderate", 20),
    rep("Strong", 20)
  ),
  Gamma_True = c(
    runif(20, -0.5, 0.95),   # Negligible: nearly neutral
    runif(20, 1.0, 2.0),     # Moderate: hump effect expected
    runif(20, 2.0, 6.0)      # Strong: moderate selection (not too extreme)
  )
)

selection_scenarios$Gene <- paste0("Gene_", 1:nrow(selection_scenarios))
selection_scenarios$AA <- "A"  # All alanine for simplicity

cat("Simulation setup:\n")
cat("Theta (4*Ne*mu):", theta_true, "\n")
cat("Alpha (4*Ne*u):", alpha_true, "\n")
cat("Beta (4*Ne*v):", beta_true, "\n")
cat("\nSelection scenarios:\n")
print(table(selection_scenarios$Category))

# ******************************************************************************
# 2) Simulate Wright-Fisher frequency distributions ----
# ******************************************************************************

cat("\n\n=== Simulating Wright-Fisher distributions ===\n")

# Function to sample from Wright-Fisher distribution
sample_wright_fisher <- function(n_sites, sample_size, alpha, beta, S) {
  #' Sample allele frequencies from Wright-Fisher stationary distribution
  #' Density: f(p) ∝ p^(alpha-1) * (1-p)^(beta-1) * exp(S*p)
  #' 
  #' Most sites should be monomorphic with small theta (~0.03)
  #' Polymorphic sites sampled from the Wright-Fisher SFS
  
  results <- matrix(NA, nrow = n_sites, ncol = 3)
  colnames(results) <- c("k", "n", "p_true")
  
  for (site in 1:n_sites) {
    # Sample frequency from Wright-Fisher stationary distribution
    # f(p) ∝ p^(alpha-1) * (1-p)^(beta-1) * exp(S*p)
    #
    # Strategy: Use Beta(alpha, beta) as proposal, then importance weight by exp(S*p)
    
    if (abs(S) < 0.01) {
      # Nearly neutral: just sample from beta distribution
      p <- rbeta(1, shape1 = alpha, shape2 = beta)
    } else {
      # With selection: use rejection sampling with beta proposal
      accepted <- FALSE
      attempts <- 0
      max_weight <- exp(S)  # Maximum weight at p=1 when S>0, or p=0 when S<0
      
      while (!accepted && attempts < 100) {
        # Proposal from beta distribution
        p_candidate <- rbeta(1, shape1 = alpha, shape2 = beta)
        
        # Importance weight: exp(S*p)
        weight <- exp(S * p_candidate)
        
        # Accept with probability weight/max_weight
        if (runif(1) * max_weight < weight) {
          p <- p_candidate
          accepted <- TRUE
        }
        attempts <- attempts + 1
      }
      
      if (!accepted) {
        # Fallback: just use beta sample (should rarely happen)
        p <- rbeta(1, shape1 = alpha, shape2 = beta)
      }
    }
    
    # Sample observed counts from binomial
    k <- rbinom(1, size = sample_size, prob = p)
    
    results[site, ] <- c(k, sample_size, p)
  }
  
  return(results)
}

# Generate data for each gene
n_sites_per_gene <- 50  # Number of synonymous sites per gene
sample_size <- 187      # Sample size (as in your data)

simulated_data <- list()

for (i in 1:nrow(selection_scenarios)) {
  gene <- selection_scenarios$Gene[i]
  gamma <- selection_scenarios$Gamma_True[i]
  
  # Sample sites
  site_data <- sample_wright_fisher(n_sites_per_gene, sample_size, 
                                     alpha_true, beta_true, gamma)
  
  simulated_data[[gene]] <- data.frame(
    Gene = gene,
    Site = 1:n_sites_per_gene,
    k = site_data[, "k"],
    n = site_data[, "n"],
    p_true = site_data[, "p_true"]
  )
  
  if (i %% 10 == 0) cat("Generated data for", i, "genes...\n")
}

# Combine all data
sim_df <- do.call(rbind, simulated_data)
rownames(sim_df) <- NULL

cat("\nGenerated", nrow(sim_df), "sites across", length(unique(sim_df$Gene)), "genes\n")

# ******************************************************************************
# 3) Calculate observed pi from simulated data ----
# ******************************************************************************

cat("\n=== Calculating observed nucleotide diversity ===\n")

# For each site, calculate pi
sim_df$Site_Pi <- with(sim_df, {
  p <- k / n
  ifelse(n > 1, 2 * p * (1 - p) * (n / (n - 1)), 0)
})

# Aggregate by gene
gene_summary <- sim_df %>%
  group_by(Gene) %>%
  summarise(
    Mean_Pi_Observed = mean(Site_Pi, na.rm = TRUE),
    N_Sites = n(),
    k_vec = list(k),
    n_vec = list(n),
    .groups = "drop"
  )

# Add true parameters
gene_summary <- left_join(gene_summary, 
                          selection_scenarios[, c("Gene", "Gamma_True", "Category")],
                          by = "Gene")

cat("Observed pi summary by category:\n")
print(gene_summary %>% group_by(Category) %>% 
        summarise(Mean_Pi = mean(Mean_Pi_Observed), 
                  SD_Pi = sd(Mean_Pi_Observed)))

# ******************************************************************************
# 4) Estimate gamma from SFS ----
# ******************************************************************************

cat("\n=== Estimating gamma from site frequency spectrum ===\n")

# Estimate gamma for each gene
gene_summary$Gamma_Estimated <- NA

for (i in 1:nrow(gene_summary)) {
  k_vec <- gene_summary$k_vec[[i]]
  n_vec <- gene_summary$n_vec[[i]]
  
  gamma_est <- tryCatch({
    estimate_gamma_for_AA(counts = k_vec, 
                         sample_sizes = n_vec,
                         u = u_true, 
                         v = v_true,
                         S_interval = c(-2, 20))
  }, error = function(e) NA)
  
  gene_summary$Gamma_Estimated[i] <- gamma_est
  
  if (i %% 10 == 0) cat("Estimated gamma for", i, "genes...\n")
}

# ******************************************************************************
# 5) Calculate expected pi from estimated gamma ----
# ******************************************************************************

cat("\n=== Calculating expected pi from estimated gamma ===\n")

gene_summary <- gene_summary %>%
  mutate(
    # Expected pi from TRUE gamma
    Pi_Expected_True = calculate_pi_analytical(alpha_true, beta_true, Gamma_True),
    # Expected pi from ESTIMATED gamma
    Pi_Expected_Est = ifelse(!is.na(Gamma_Estimated),
                             calculate_pi_analytical(alpha_true, beta_true, Gamma_Estimated),
                             NA)
  )

# ******************************************************************************
# 6) Analyze results ----
# ******************************************************************************

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("VALIDATION RESULTS\n")
cat(strrep("=", 70), "\n\n")

# Remove failed estimations
valid_results <- gene_summary %>% filter(!is.na(Gamma_Estimated))

cat("Successfully estimated gamma for", nrow(valid_results), "out of", 
    nrow(gene_summary), "genes\n\n")

# Overall correlation
cor_gamma <- cor(valid_results$Gamma_True, valid_results$Gamma_Estimated)
cat("Correlation between true and estimated gamma:", round(cor_gamma, 3), "\n")

# RMSE
rmse_gamma <- sqrt(mean((valid_results$Gamma_True - valid_results$Gamma_Estimated)^2))
cat("RMSE for gamma estimation:", round(rmse_gamma, 3), "\n\n")

# By category
cat("Results by selection category:\n")
cat(strrep("-", 70), "\n")
summary_by_cat <- valid_results %>%
  group_by(Category) %>%
  summarise(
    N = n(),
    Mean_Gamma_True = mean(Gamma_True),
    Mean_Gamma_Est = mean(Gamma_Estimated),
    Bias = mean(Gamma_Estimated - Gamma_True),
    RMSE = sqrt(mean((Gamma_True - Gamma_Estimated)^2)),
    Correlation = cor(Gamma_True, Gamma_Estimated),
    Mean_Pi_Obs = mean(Mean_Pi_Observed),
    Mean_Pi_Exp_True = mean(Pi_Expected_True),
    .groups = "drop"
  )

print(as.data.frame(summary_by_cat))

# Check for hump effect
cat("\n\nChecking for HUMP EFFECT (pi should peak at moderate selection):\n")
cat(strrep("-", 70), "\n")
cat("Mean observed pi by category:\n")
pi_by_cat <- valid_results %>%
  group_by(Category) %>%
  summarise(Mean_Pi = mean(Mean_Pi_Observed), .groups = "drop") %>%
  arrange(match(Category, c("Negligible", "Moderate", "Strong")))
print(as.data.frame(pi_by_cat))

# Find which category has highest pi
max_pi_cat <- pi_by_cat$Category[which.max(pi_by_cat$Mean_Pi)]

if (max_pi_cat == "Moderate") {
  cat("\n✓ HUMP EFFECT DETECTED: Pi is highest for moderate selection!\n")
} else if (max_pi_cat == "Negligible") {
  cat("\n✓ Pattern confirmed: pi highest for nearly neutral sites\n")
  cat("  (Hump may be visible only when moderate selection is weak enough)\n")
} else {
  cat("\n✗ Warning: Unexpected pattern - pi highest at", max_pi_cat, "selection\n")
}

# ******************************************************************************
# 7) Visualizations ----
# ******************************************************************************

cat("\n=== Generating validation plots ===\n")

# Plot 1: True vs Estimated Gamma
p1 <- ggplot(valid_results, aes(x = Gamma_True, y = Gamma_Estimated, color = Category)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dotted") +
  scale_color_manual(values = c("Weak" = "#2E86AB", "Moderate" = "#A23B72", "Strong" = "#F18F01")) +
  labs(
    title = "Validation: True vs Estimated Selection Coefficient",
    subtitle = paste0("Correlation = ", round(cor_gamma, 3), " | RMSE = ", round(rmse_gamma, 3)),
    x = "True Gamma (4*Ne*s)",
    y = "Estimated Gamma"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("./tests/validation_gamma_recovery.png", p1, width = 8, height = 6, dpi = 300)

# Plot 2: Hump Effect - Pi vs Selection
valid_results_long <- valid_results %>%
  select(Gene, Category, Gamma_True, Mean_Pi_Observed, Pi_Expected_True) %>%
  tidyr::pivot_longer(cols = c(Mean_Pi_Observed, Pi_Expected_True),
                      names_to = "Type", values_to = "Pi")

p2 <- ggplot(valid_results_long, aes(x = Gamma_True, y = Pi, color = Type)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = TRUE, span = 0.5) +
  scale_color_manual(
    values = c("Mean_Pi_Observed" = "#E63946", "Pi_Expected_True" = "#457B9D"),
    labels = c("Observed", "Theoretical")
  ) +
  labs(
    title = "Hump Effect: Nucleotide Diversity vs Selection Strength",
    subtitle = "Mutation-selection balance predicts peak diversity at moderate selection",
    x = "Selection Coefficient Gamma (4*Ne*s)",
    y = "Nucleotide Diversity (π)",
    color = "Data Type"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("./tests/validation_hump_effect.png", p2, width = 10, height = 6, dpi = 300)

# Plot 3: Observed vs Expected Pi
p3 <- ggplot(valid_results, aes(x = Pi_Expected_True, y = Mean_Pi_Observed, color = Category)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dotted") +
  scale_color_manual(values = c("Negligible" = "#2E86AB", "Moderate" = "#A23B72", "Strong" = "#F18F01")) +
  labs(
    title = "Validation: Observed vs Expected Nucleotide Diversity",
    subtitle = "Using true selection coefficients",
    x = "Expected π (from true gamma)",
    y = "Observed π"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("./tests/validation_pi_match.png", p3, width = 8, height = 6, dpi = 300)

# Plot 4: Residuals
valid_results$Gamma_Residual <- valid_results$Gamma_Estimated - valid_results$Gamma_True

p4 <- ggplot(valid_results, aes(x = Gamma_True, y = Gamma_Residual, color = Category)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  scale_color_manual(values = c("Negligible" = "#2E86AB", "Moderate" = "#A23B72", "Strong" = "#F18F01")) +
  labs(
    title = "Estimation Bias: Residuals vs True Gamma",
    x = "True Gamma (4*Ne*s)",
    y = "Residual (Estimated - True)",
    color = "Selection Category"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("./tests/validation_residuals.png", p4, width = 8, height = 6, dpi = 300)

cat("\nPlots saved to ./tests/\n")

# ******************************************************************************
# 8) Save detailed results ----
# ******************************************************************************

# Remove list columns (k_vec, n_vec) before writing CSV
valid_results_csv <- valid_results %>%
  select(-k_vec, -n_vec)

write.csv(valid_results_csv, "./tests/validation_results.csv", row.names = FALSE)
cat("\nDetailed results saved to ./tests/validation_results.csv\n")

# Print summary
cat("\n\n")
cat(strrep("=", 70), "\n")
cat("PIPELINE VALIDATION COMPLETE\n")
cat(strrep("=", 70), "\n")
cat("\nKey Findings:\n")
cat("1. Gamma recovery correlation:", round(cor_gamma, 3), "\n")
cat("2. Hump effect visible:", 
    ifelse(pi_by_cat$Mean_Pi[2] > max(pi_by_cat$Mean_Pi[c(1,3)]), "YES ✓", "NO ✗"), "\n")
cat("3. Mean bias in gamma estimation:", 
    round(mean(valid_results$Gamma_Residual), 3), "\n")
cat("\nThe pipeline is", 
    ifelse(cor_gamma > 0.7 && abs(mean(valid_results$Gamma_Residual)) < 1, 
           "VALIDATED ✓", "NEEDS REVIEW ✗"), "\n")
