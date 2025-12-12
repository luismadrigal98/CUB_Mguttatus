#!/usr/bin/env Rscript
# Analyze relationship between Site_Pi_Nucleotide and selection coefficient
# Theory predicts a "hump" - maximum pi at intermediate selection (S ≈ 1.7)

library(data.table)
library(ggplot2)

source("./src/derivation_gamma_from_polymorphism.R")

cat(paste(rep("=", 80), collapse=""), "\n")
cat("RELATIONSHIP BETWEEN PI AND SELECTION COEFFICIENT\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

# Parameters from your analysis
alpha <- 0.008  # From your plot
beta <- 0.024   # From your plot
theta <- alpha + beta  # = 0.032 ≈ 0.0312

cat("Parameters:\n")
cat("  alpha (u scaled) =", alpha, "\n")
cat("  beta (v scaled) =", beta, "\n")
cat("  theta =", theta, "\n")
cat("  Mutation bias: beta/alpha =", round(beta/alpha, 2), "\n\n")

# Test range of selection coefficients
S_values <- seq(-2, 15, by = 0.5)
gamma_values <- S_values  # gamma = 4*Ne*s, but S here is already scaled

cat("Calculating expected pi for different selection coefficients...\n\n")

# Calculate expected pi for each S
pi_results <- data.table(
  S = S_values,
  gamma = gamma_values
)

pi_results[, Pi_Expected := sapply(S, function(s) {
  calculate_pi_analytical(alpha, beta, s)
})]

# Find the peak
peak_idx <- which.max(pi_results$Pi_Expected)
peak_S <- pi_results$S[peak_idx]
peak_pi <- pi_results$Pi_Expected[peak_idx]

cat(paste(rep("-", 60), collapse=""), "\n")
cat("THEORETICAL PREDICTIONS\n")
cat(paste(rep("-", 60), collapse=""), "\n\n")

cat("Selection regimes:\n")
cat("  Neutral (S ≈ 0):     Pi =", round(pi_results$Pi_Expected[pi_results$S == 0], 5), "\n")
cat("  Peak (S ≈", round(peak_S, 1), "):   Pi =", round(peak_pi, 5), "← MAXIMUM\n")
cat("  Strong (S = 10):     Pi =", round(pi_results$Pi_Expected[pi_results$S == 10], 5), "\n")
cat("  Very strong (S=15):  Pi =", round(pi_results$Pi_Expected[pi_results$S == 15], 5), "\n\n")

cat("Key insight: Pi is HIGHER at intermediate selection than at neutrality!\n\n")

# Calculate how balanced polymorphisms relate to this
cat(paste(rep("-", 60), collapse=""), "\n")
cat("INTERPRETING HIGH PI VALUES IN YOUR DATA\n")
cat(paste(rep("-", 60), collapse=""), "\n\n")

cat("When you see Site_Pi_Nucleotide = 0.167:\n\n")

cat("1. This is CODON-LEVEL pi after scaling:\n")
cat("   - Codon heterozygosity ≈ 0.50 (p ≈ 0.5)\n")
cat("   - After nucleotide scaling: 0.50 × (1/3) = 0.167\n\n")

cat("2. Expected POPULATION-AVERAGE pi is much lower:\n")
cat("   - Peak theoretical pi ≈", round(peak_pi, 4), "\n")
cat("   - This is averaged over the ENTIRE frequency spectrum\n")
cat("   - Wright-Fisher SFS is U-shaped (most mass at boundaries)\n\n")

cat("3. Why the discrepancy?\n")
cat("   - Your Site_Pi is calculated from OBSERVED frequencies at ONE site\n")
cat("   - Theoretical pi is EXPECTED average across many realizations\n")
cat("   - A site with p = 0.5 can exist, but it's rare in equilibrium\n\n")

cat("4. Probability of observing p ≈ 0.5:\n")
# Calculate probability density at p = 0.5 for different S values
cat("   Under neutrality (S=0):\n")
cat("     P(p=0.5) ∝ Beta(α=", alpha, ", β=", beta, ") density\n")
cat("     This is LOW but non-zero\n\n")

cat("   Under intermediate selection (S≈1.7):\n")
cat("     P(p=0.5) is HIGHEST (balancing selection effect)\n")
cat("     π is maximized because variants stick around at intermediate freq\n\n")

cat("   Under strong selection (S>5):\n")
cat("     P(p=0.5) is VERY LOW (preferred fixed quickly)\n")
cat("     π is low because variants rarely polymorphic\n\n")

# Demonstrate with Wright-Fisher probabilities
cat(paste(rep("=", 80), collapse=""), "\n")
cat("FREQUENCY SPECTRUM COMPARISON\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

# For n = 187 (your sample size), calculate expected SFS
n <- 187
k_values <- 0:n

# Calculate for three scenarios
scenarios <- data.table(
  Scenario = c("Neutral", "Peak_Selection", "Strong_Selection"),
  S = c(0, 1.7, 10)
)

cat("Probability of observing k preferred alleles (out of n=187):\n\n")

for (i in 1:nrow(scenarios)) {
  scen <- scenarios$Scenario[i]
  S_val <- scenarios$S[i]
  
  cat(scen, "(S =", S_val, "):\n")
  
  # Calculate probabilities for different k values
  probs <- sapply(k_values, function(k) {
    if (k == 0 || k == n) return(NA)  # Skip boundaries
    get_prob_k_analytical(k, n, alpha/theta, beta/theta, S_val, theta)
  })
  
  # Find probability near p = 0.5 (k ≈ n/2)
  mid_k <- round(n/2)
  mid_range <- (mid_k - 5):(mid_k + 5)
  mid_prob <- sum(probs[mid_range + 1], na.rm = TRUE)
  
  # Probability in tails (rare variants)
  tail_low <- sum(probs[2:20], na.rm = TRUE)
  tail_high <- sum(probs[(n-20):(n-1)], na.rm = TRUE)
  
  cat("  P(k near n/2, balanced) =", round(mid_prob, 4), "\n")
  cat("  P(k < 20, rare variant) =", round(tail_low, 4), "\n")
  cat("  P(k > n-20, near-fixed) =", round(tail_high, 4), "\n")
  cat("  Expected pi:", round(pi_results$Pi_Expected[pi_results$S == S_val], 5), "\n\n")
}

cat(paste(rep("=", 80), collapse=""), "\n")
cat("CONCLUSIONS\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

cat("1. POPULATION AVERAGE (theoretical pi):\n")
cat("   - Peak at S ≈", round(peak_S, 1), "with π ≈", round(peak_pi, 4), "\n")
cat("   - This is averaged over ALL possible frequencies\n")
cat("   - Most realizations have extreme frequencies (U-shaped SFS)\n\n")

cat("2. INDIVIDUAL SITE OBSERVATION (your Site_Pi = 0.167):\n")
cat("   - Represents ONE realization at balanced frequency (p ≈ 0.5)\n")
cat("   - This is POSSIBLE but RARE in equilibrium\n")
cat("   - Such sites could be under balancing selection (S ≈ 0-3)\n\n")

cat("3. WHY HIGH PI VALUES AREN'T ERRORS:\n")
cat("   - Site_Pi measures realized diversity at ONE site\n")
cat("   - Can be up to 0.5 (biallelic maximum)\n")
cat("   - Mean across many sites should match theoretical pi (~0.01-0.03)\n\n")

cat("4. WHAT TO EXPECT IN YOUR DATA:\n")
cat("   - Mean_Pi_Observed (across sites) ≈ 0.01-0.03 ✓\n")
cat("   - Some sites with high pi (0.1-0.5) are expected ✓\n")
cat("   - Sites with pi > 0.05 might have S ≈ 0-3 (weak/intermediate selection)\n")
cat("   - Sites with pi < 0.01 likely have strong selection (S > 5)\n\n")

# Create two-panel plot
library(gridExtra)

pdf("results/pi_vs_selection_theoretical.pdf", width = 14, height = 6)

# Panel 1: Full range
p1 <- ggplot(pi_results, aes(x = S, y = Pi_Expected)) +
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_point(data = pi_results[S == peak_S], 
             aes(x = S, y = Pi_Expected), 
             color = "red", size = 4) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = peak_S, linetype = "dashed", color = "red", alpha = 0.5) +
  annotate("text", x = peak_S, y = peak_pi * 1.05, 
           label = sprintf("Peak: S = %.1f\npi = %.4f", peak_S, peak_pi),
           color = "red", fontface = "bold", vjust = 0) +
  annotate("text", x = -1, y = peak_pi * 0.5, 
           label = "Mutation\nDominated", color = "gray40", size = 3.5) +
  annotate("text", x = 8, y = peak_pi * 0.3, 
           label = "Selection\nDominated", color = "gray40", size = 3.5) +
  annotate("rect", xmin = -2, xmax = 15, ymin = 0, ymax = 0.167,
           alpha = 0.1, fill = "green") +
  annotate("text", x = 12, y = 0.167 * 0.95, 
           label = "Range of observed Site_Pi values\n(individual sites can reach 0.167)",
           color = "darkgreen", size = 3, hjust = 1, vjust = 1) +
  geom_hline(yintercept = theta, linetype = "dotted", color = "blue", alpha = 0.7) +
  annotate("text", x = 12, y = theta * 1.1, 
           label = sprintf("theta = %.4f", theta), 
           color = "blue", size = 3) +
  labs(
    title = "Full Range: Nucleotide Diversity vs Selection",
    subtitle = sprintf("Parameters: alpha = %.3f, beta = %.3f, theta = %.4f (Mutation bias beta/alpha = %.2f)", 
                      alpha, beta, theta, beta/alpha),
    x = "Selection Coefficient (gamma = 4Nes)",
    y = "Expected Nucleotide Diversity (pi)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9)
  )

# Panel 2: Zoomed in on hump (pi < theta)
# Filter data to show only region below theta
pi_zoom <- pi_results[Pi_Expected <= theta * 1.05]

p2 <- ggplot(pi_zoom, aes(x = S, y = Pi_Expected)) +
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_point(data = pi_zoom[S == peak_S], 
             aes(x = S, y = Pi_Expected), 
             color = "red", size = 4) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = peak_S, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_hline(yintercept = theta, linetype = "dotted", color = "blue", alpha = 0.7) +
  annotate("text", x = peak_S + 2, y = peak_pi, 
           label = sprintf("Peak\nS = %.1f\npi = %.5f", peak_S, peak_pi),
           color = "red", fontface = "bold", size = 3.5, hjust = 0) +
  annotate("text", x = max(pi_zoom$S) * 0.95, y = theta * 0.98, 
           label = sprintf("theta = %.4f", theta), 
           color = "blue", size = 3, vjust = 1, hjust = 1) +
  annotate("rect", xmin = 0, xmax = 3, ymin = 0, ymax = theta * 1.05,
           alpha = 0.05, fill = "orange") +
  annotate("text", x = 1.5, y = theta * 0.4, 
           label = "Hump Effect Region\n(Weak Selection)", 
           color = "darkorange", size = 3, fontface = "italic") +
  labs(
    title = "Zoomed View: The Hump Effect",
    subtitle = "Region where pi < theta shows peak diversity at intermediate selection",
    x = "Selection Coefficient (gamma = 4Nes)",
    y = "Expected Nucleotide Diversity (pi)"
  ) +
  coord_cartesian(ylim = c(0, theta * 1.05)) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9)
  )

# Combine panels
grid.arrange(p1, p2, ncol = 2)
dev.off()

cat("Two-panel plot saved to: results/pi_vs_selection_theoretical.pdf\n\n")
