#' ROC Model Codon Trajectory Plot
#' 
#' Functions to visualize predicted vs observed codon frequencies across
#' expression levels using the AnaCoDa ROC model.
#'
#' ROC model: P(codon_i | phi) = exp(-dM_i - dEta_i * phi) / Z
#'
#' @author Luis Javier Madrigal-Roca
#' @date 2024-12-04

# ==============================================================================
# LOAD CSP PARAMETERS
# ==============================================================================

#' Load CSP parameters from AnaCoDa output files
#'
#' @param mutation_file Path to Cluster_X_Mutation.csv
#' @param selection_file Path to Cluster_X_Selection.csv
#' @return Data frame with AA, Codon, dM, dEta
load_csp_parameters <- function(mutation_file, selection_file) {
  
  dM_df <- read.csv(mutation_file, stringsAsFactors = FALSE)
  dEta_df <- read.csv(selection_file, stringsAsFactors = FALSE)
  
  csp <- merge(
    dM_df[, c("AA", "Codon", "Mean")],
    dEta_df[, c("AA", "Codon", "Mean")],
    by = c("AA", "Codon"),
    suffixes = c("_dM", "_dEta")
  )
  
  names(csp)[names(csp) == "Mean_dM"] <- "dM"
  names(csp)[names(csp) == "Mean_dEta"] <- "dEta"
  
  # Identify reference codons (dM = 0 and dEta = 0)
  csp$is_reference <- (csp$dM == 0 & csp$dEta == 0)
  
  message(sprintf("Loaded CSP for %d codons, %d amino acids", 
                  nrow(csp), length(unique(csp$AA))))
  
  return(csp)
}

# ==============================================================================
# PREDICT CODON PROBABILITIES
# ==============================================================================

#' Calculate predicted codon probabilities for a range of phi values
#'
#' @param phi_values Numeric vector of phi (expression) values
#' @param csp_df CSP parameters data frame
#' @return Data frame with phi, AA, Codon, predicted_prob
predict_across_phi_range <- function(phi_values, csp_df) {
  
  results <- list()
  
  for (phi in phi_values) {
    # For each AA, calculate multinomial probabilities
    aa_list <- split(csp_df, csp_df$AA)
    
    pred_list <- lapply(aa_list, function(aa_df) {
      # log P(codon) = -dM - dEta * phi - log(Z)
      log_unnorm <- -aa_df$dM - aa_df$dEta * phi
      
      # Log-sum-exp for numerical stability
      max_log <- max(log_unnorm)
      log_Z <- max_log + log(sum(exp(log_unnorm - max_log)))
      
      aa_df$predicted_prob <- exp(log_unnorm - log_Z)
      aa_df$phi <- phi
      
      aa_df[, c("AA", "Codon", "phi", "predicted_prob", "is_reference")]
    })
    
    results[[length(results) + 1]] <- do.call(rbind, pred_list)
  }
  
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  return(out)
}

# ==============================================================================
# MAP AMINO ACIDS TO ANACODA CONVENTION
# ==============================================================================

#' Map observed codon data to AnaCoDa AA convention
#' 
#' AnaCoDa uses "Z" for Serine AGN codons (AGC, AGT)
#'
#' @param codon_df Data frame with AA and Codon columns
#' @return Data frame with AA_anacoda column added
map_aa_to_anacoda <- function(codon_df) {
  
  codon_df$AA_anacoda <- codon_df$AA
  
  # AGC and AGT are coded as "Z" in AnaCoDa
  codon_df$AA_anacoda[codon_df$Codon %in% c("AGC", "AGT")] <- "Z"
  
  return(codon_df)
}

# ==============================================================================
# CALCULATE OBSERVED FREQUENCIES BY EXPRESSION BIN
# ==============================================================================

#' Calculate observed codon frequencies in expression bins
#'
#' @param obs_data Data frame with Gene, Codon, AA_anacoda, Observed_freq, Exp_log10
#' @param n_bins Number of expression bins (default: 10)
#' @return Data frame with binned observed frequencies and SDs
calculate_observed_by_bin <- function(obs_data, n_bins = 10) {
  
  # Create expression bins
  breaks <- quantile(obs_data$Exp_log10, probs = seq(0, 1, length.out = n_bins + 1), 
                     na.rm = TRUE)
  
  # Ensure unique breaks
  breaks <- unique(breaks)
  n_actual_bins <- length(breaks) - 1
  
  obs_data$Exp_bin <- cut(obs_data$Exp_log10, breaks = breaks, 
                          include.lowest = TRUE, labels = FALSE)
  
  # Calculate mean and SD per codon per bin
  obs_summary <- obs_data |>
    dplyr::filter(!is.na(Exp_bin)) |>
    dplyr::group_by(AA_anacoda, Codon, Exp_bin) |>
    dplyr::summarize(
      Observed_mean = mean(Observed_freq, na.rm = TRUE),
      Observed_sd = sd(Observed_freq, na.rm = TRUE),
      Exp_mean = mean(Exp_log10, na.rm = TRUE),
      n_genes = dplyr::n(),
      .groups = "drop"
    )
  
  return(obs_summary)
}

# ==============================================================================
# MAIN PLOTTING FUNCTION
# ==============================================================================

#' Plot codon frequency trajectories across expression
#'
#' Creates a faceted plot showing predicted (lines) vs observed (points + error bars)
#' codon frequencies across expression levels, similar to AnaCoDa diagnostic plots.
#'
#' @param csp_df CSP parameters (from load_csp_parameters)
#' @param obs_binned Observed frequencies by bin (from calculate_observed_by_bin)
#' @param phi_range Range of phi values for prediction curves (default: c(-4, 5))
#' @param phi_n Number of points for smooth prediction curves (default: 100)
#' @param output_file Path to save PDF (optional)
#' @return ggplot object
plot_codon_trajectories <- function(csp_df, obs_binned, 
                                     phi_range = c(-4, 5), 
                                     phi_n = 100,
                                     output_file = NULL) {
  
  require(ggplot2)
  
  # Generate prediction curves
  phi_seq <- seq(phi_range[1], phi_range[2], length.out = phi_n)
  pred_curves <- predict_across_phi_range(phi_seq, csp_df)
  
  # Mark reference codons with asterisk in legend
  pred_curves$Codon_label <- ifelse(pred_curves$is_reference, 
                                     paste0(pred_curves$Codon, "*"),
                                     pred_curves$Codon)
  
  # Same for observed data
  obs_binned <- merge(obs_binned, 
                      unique(csp_df[, c("AA", "Codon", "is_reference")]),
                      by.x = c("AA_anacoda", "Codon"),
                      by.y = c("AA", "Codon"),
                      all.x = TRUE)
  
  obs_binned$Codon_label <- ifelse(obs_binned$is_reference, 
                                    paste0(obs_binned$Codon, "*"),
                                    obs_binned$Codon)
  
  # Create plot
  p <- ggplot() +
    # Prediction curves (lines)
    geom_line(data = pred_curves,
              aes(x = phi, y = predicted_prob, color = Codon_label),
              linewidth = 1) +
    # Observed points
    geom_point(data = obs_binned,
               aes(x = Exp_mean, y = Observed_mean, color = Codon_label),
               size = 2) +
    # Error bars (SD)
    geom_errorbar(data = obs_binned,
                  aes(x = Exp_mean, 
                      ymin = pmax(0, Observed_mean - Observed_sd),
                      ymax = pmin(1, Observed_mean + Observed_sd),
                      color = Codon_label),
                  width = 0.1, alpha = 0.7) +
    # Facet by amino acid
    facet_wrap(~ AA_anacoda, scales = "free_y", ncol = 5) +
    # Axis limits
    scale_y_continuous(limits = c(0, 1)) +
    # Labels
    labs(
      x = "Expression (log10 scale)",
      y = "Codon Frequency",
      title = "ROC Model: Predicted vs Observed Codon Usage",
      subtitle = "Lines = Model prediction | Points = Observed (mean ± SD) | * = Reference codon",
      color = "Codon"
    ) +
    # Theme
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      strip.text = element_text(size = 12, face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(nrow = 2, override.aes = list(linewidth = 2)))
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 16, height = 14)
    message(sprintf("Saved: %s", output_file))
  }
  
  return(p)
}

# ==============================================================================
# CONVENIENCE WRAPPER
# ==============================================================================

#' Run complete codon trajectory analysis
#'
#' @param mutation_file Path to AnaCoDa mutation CSV
#' @param selection_file Path to AnaCoDa selection CSV
#' @param codon_freq_df Data frame with Gene, Codon, AA, Observed_freq columns
#' @param expression_df Data frame with Gene, Exp_log10 columns
#' @param output_file Path to save plot (optional)
#' @param n_bins Number of expression bins for observed data
#' @return List with plot, csp_params, and trajectory data
run_trajectory_analysis <- function(mutation_file, selection_file,
                                     codon_freq_df, expression_df,
                                     output_file = NULL, n_bins = 10) {
  
  # 1. Load CSP parameters
  csp <- load_csp_parameters(mutation_file, selection_file)
  
  # 2. Merge codon frequencies with expression
  obs_data <- merge(codon_freq_df, expression_df, by = "Gene")
  
  # 3. Map to AnaCoDa AA convention
  obs_data <- map_aa_to_anacoda(obs_data)
  
  # 4. Calculate observed frequencies by expression bin
  obs_binned <- calculate_observed_by_bin(obs_data, n_bins = n_bins)
  
  message(sprintf("Binned observations: %d AA-codon-bin combinations", nrow(obs_binned)))
  
  # 5. Create trajectory plot
  p <- plot_codon_trajectories(csp, obs_binned, output_file = output_file)
  
  return(list(
    plot = p,
    csp = csp,
    observed_binned = obs_binned
  ))
}
