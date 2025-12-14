# --- Simulation Parameters ---
# Total Mutation Rate (theta = 4Ne*mu) approx 0.015 (typical for plants)
theta_total <- 0.032

# Bias: Mutation favors Unpreferred (v > u) by 3x
# This matches the Mimulus AT-bias
bias_factor <- 3 
alpha_sim <- theta_total / (1 + bias_factor) # u (Unpref -> Pref)
beta_sim  <- theta_total - alpha_sim         # v (Pref -> Unpref)

# Create a sequence of Selection Coefficients (S)
# From Neutral (0) to Strong (15)
S_values <- seq(0, 15, by = 0.1)

# Calculate Pi for each S
pi_values <- sapply(S_values, function(s) {
  calculate_pi_analytical(alpha_sim, beta_sim, s)
})

# Create Dataframe for Plotting
plot_df <- data.frame(S = S_values, Pi = pi_values)

# Find the Peak for annotation
peak_S <- plot_df$S[which.max(plot_df$Pi)]
peak_Pi <- max(plot_df$Pi)

# --- Generate the Plot ---
ggplot(plot_df, aes(x = S, y = Pi)) +
  # The Curve
  geom_line(color = "#2c3e50", linewidth = 1.2) +
  
  # Area Shading
  geom_area(fill = "#3498db", alpha = 0.2) +
  
  # Vertical line at Peak
  geom_vline(xintercept = peak_S, linetype = "dashed", color = "red") +
  
  # Annotations
  annotate("text", x = 0.45, y = min(pi_values), label = "Mutation\nDominated", 
           vjust = 0, hjust = 0.5, color = "gray40", fontface = "italic") +
  
  annotate("text", x = peak_S + 2, y = peak_Pi + 0.00005, 
           label = paste0("Peak Diversity\n(S = ", peak_S, ")"), 
           color = "red", fontface = "bold", size = 3.5) +
  
  annotate("text", x = 12, y = min(pi_values), label = "Selection\nDominated", 
           vjust = 0.6, hjust = 0.5, color = "gray40", fontface = "italic") +
  
  # Labels
  labs(title = "Theoretical increase of Nucleotide Diversity (Pi) under Weak Selection",
       subtitle = paste0("Parameters: Alpha (u) = ", round(alpha_sim, 4), 
                         ", Beta (v) = ", round(beta_sim, 4), " (Mutation Bias v > u)"),
       x = expression(paste("Selection Coefficient (", gamma == 4*N[e]*s, ")")),
       y = expression(paste("Expected Nucleotide Diversity (", pi, ")"))) +
  
  theme_bw() +
  theme(plot.title = element_text(face = "bold"),
        panel.grid = element_blank())
