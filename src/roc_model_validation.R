#' ROC Model Codon Trajectory Plot - AnaCoDa Style
#' 
#' Functions to visualize predicted vs observed codon frequencies across
#' expression levels using the AnaCoDa ROC model.
#'
#' ROC model: P(codon_i | phi) = exp(-dM_i - dEta_i * phi) / Z
#' where phi is on linear scale (not log)
#'
#' @author Luis Javier Madrigal-Roca
#' @date 2024-12-04

# ==============================================================================
# AnaCoDa COLOR SCHEME
# ==============================================================================

# Exact colors from AnaCoDa's colorSchemes.R for each codon
.codonColors <- list(
  GCA = "blue", GCC = "darkorange", GCG = "purple", GCT = "green4",
  TGC = "darkorange", TGT = "green4",
  GAC = "darkorange", GAT = "green4",
  GAA = "blue", GAG = "purple",

  TTC = "darkorange", TTT = "green4",
  GGA = "blue", GGC = "darkorange", GGG = "purple", GGT = "green4",
  CAC = "darkorange", CAT = "green4",
  ATA = "blue", ATC = "darkorange", ATT = "green4",
  AAA = "blue", AAG = "purple",
  CTA = "blue", CTC = "darkorange", CTG = "purple", CTT = "green4",
  TTA = "darkturquoise", TTG = "deeppink3",
  AAC = "darkorange", AAT = "green4",
  CCA = "blue", CCC = "darkorange", CCG = "purple", CCT = "green4",
  CAA = "blue", CAG = "purple",
  CGA = "blue", CGC = "darkorange", CGG = "purple", CGT = "green4",
  AGA = "darkturquoise", AGG = "deeppink3",
  TCA = "blue", TCC = "darkorange", TCG = "purple", TCT = "green4",
  ACA = "blue", ACC = "darkorange", ACG = "purple", ACT = "green4",
  GTA = "blue", GTC = "darkorange", GTG = "purple", GTT = "green4",
  TAC = "darkorange", TAT = "green4",
  AGC = "darkorange", AGT = "green4",
  TGG = "blue"
)

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
  
  # Identify OPTIMAL/PREFERRED codon: most negative dEta (strongest selection)
  # This is what AnaCoDa marks with asterisk
  csp$is_optimal <- FALSE
  for (aa in unique(csp$AA)) {
    aa_idx <- which(csp$AA == aa)
    min_dEta_idx <- aa_idx[which.min(csp$dEta[aa_idx])]
    csp$is_optimal[min_dEta_idx] <- TRUE
  }
  
  message(sprintf("Loaded CSP for %d codons, %d amino acids", 
                  nrow(csp), length(unique(csp$AA))))
  
  return(csp)
}

# ==============================================================================
# PREDICT CODON PROBABILITIES
# ==============================================================================

#' Calculate predicted codon probabilities for a range of phi values
#'
#' Uses the ROC multinomial model:
#'   P(codon_i | phi) = exp(-dM_i - dEta_i * phi) / Z
#' 
#' IMPORTANT: phi should be on LINEAR scale (not log10)
#' For plotting, we use log10(phi) on x-axis but compute model with linear phi
#'
#' @param phi_log10_values Numeric vector of log10(phi) values (for x-axis)
#' @param csp_df CSP parameters data frame
#' @return Data frame with phi_log10, AA, Codon, predicted_prob
predict_across_phi_range <- function(phi_log10_values, csp_df) {
  
  results <- list()
  
  for (phi_log10 in phi_log10_values) {
    # Convert log10(phi) to linear phi for model calculation
    phi_linear <- 10^phi_log10
    
    # For each AA, calculate multinomial probabilities
    aa_list <- split(csp_df, csp_df$AA)
    
    pred_list <- lapply(aa_list, function(aa_df) {
      # ROC model: log P(codon) = -dM - dEta * phi_linear - log(Z)
      log_unnorm <- -aa_df$dM - aa_df$dEta * phi_linear
      
      # Log-sum-exp for numerical stability
      max_log <- max(log_unnorm)
      log_Z <- max_log + log(sum(exp(log_unnorm - max_log)))
      
      aa_df$predicted_prob <- exp(log_unnorm - log_Z)
      aa_df$phi <- phi_log10  # Store log10 scale for plotting
      
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

#' Map observed codon data to AnaCoDa AA convention and recalculate frequencies
#' 
#' AnaCoDa uses "Z" for Serine AGN codons (AGC, AGT).
#' This function also RECALCULATES frequencies for Z codons within the Z group only,
#' since the original frequencies are computed relative to all 6 serine codons.
#'
#' @param codon_df Data frame with Gene, AA, Codon, Count columns
#' @return Data frame with AA_anacoda column and corrected Observed_freq
map_aa_to_anacoda <- function(codon_df) {
  
  # First, map AA to AnaCoDa convention
  codon_df$AA_anacoda <- codon_df$AA
  
  # AGC and AGT are coded as "Z" in AnaCoDa (separate from TCN serines)
  codon_df$AA_anacoda[codon_df$Codon %in% c("AGC", "AGT")] <- "Z"
  
  # CRITICAL: Recalculate frequencies within each AA_anacoda group

  # This fixes the Z (Ser2) problem where frequencies were computed
  # relative to all 6 serine codons instead of just AGC+AGT
  if ("Count" %in% names(codon_df)) {
    codon_df <- codon_df |>
      dplyr::group_by(Gene, AA_anacoda) |>
      dplyr::mutate(
        AA_total_anacoda = sum(Count, na.rm = TRUE),
        Observed_freq = ifelse(AA_total_anacoda > 0, 
                               Count / AA_total_anacoda, 
                               NA_real_)
      ) |>
      dplyr::ungroup()
    
    message("Recalculated codon frequencies within AnaCoDa AA groups (Z = AGC+AGT only)")
  }
  
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
# MAIN PLOTTING FUNCTION - AnaCoDa Style
# ==============================================================================

#' Plot codon frequency trajectories for a single amino acid
#'
#' @param aa_code Amino acid code (single letter)
#' @param pred_df Predictions for this AA
#' @param obs_df Observed data for this AA
#' @param codon_colors Named vector of colors for codons
#' @return ggplot object
plot_single_aa <- function(aa_code, pred_df, obs_df, codon_colors) {
  
  require(ggplot2)
  
  p <- ggplot() +
    # Predicted curves (smooth lines)
    geom_line(data = pred_df,
              aes(x = phi, y = predicted_prob, color = Codon_label),
              linewidth = 1.2) +
    # Observed points
    geom_point(data = obs_df,
               aes(x = Exp_mean, y = Observed_mean, color = Codon_label),
               size = 1.5) +
    # Error bars
    geom_errorbar(data = obs_df,
                  aes(x = Exp_mean,
                      ymin = pmax(0, Observed_mean - Observed_sd),
                      ymax = pmin(1, Observed_mean + Observed_sd),
                      color = Codon_label),
                  width = 0.05, linewidth = 0.7) +
    scale_color_manual(values = codon_colors) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = aa_code,
      x = NULL,
      y = NULL,
      color = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = 9),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  return(p)
}

#' Plot codon frequency trajectories across expression - AnaCoDa style
#'
#' Creates a multi-panel plot with each amino acid having its own color legend.
#' Predicted curves show the multinomial probability evolution across phi.
#' Asterisk marks the OPTIMAL (most preferred) codon, i.e. most negative dEta.
#'
#' @param csp_df CSP parameters (from load_csp_parameters)
#' @param obs_binned Observed frequencies by bin (from calculate_observed_by_bin)
#' @param phi_range Range of log10 phi values for prediction curves. 
#'                  Extended beyond data to show full sigmoid (-1 to 3 by default).
#' @param phi_n Number of points for smooth prediction curves
#' @param output_file Path to save PDF (optional)
#' @return List of ggplot objects (one per AA)
plot_codon_trajectories <- function(csp_df, obs_binned, 
                                     phi_range = NULL, 
                                     phi_n = 100,
                                     output_file = NULL) {
  
  require(ggplot2)
  require(gridExtra)
  
  # Determine phi range dynamically
  # Extend modestly beyond data to show sigmoid trends, but stay reasonable
  if (is.null(phi_range)) {
    data_range <- range(obs_binned$Exp_mean, na.rm = TRUE)
    data_width <- diff(data_range)
    
    # Extend by 30% of data width on each side, with reasonable bounds
    extension <- max(0.5, data_width * 0.5)
    phi_range <- c(
      data_range[1] - extension,
      data_range[2] + extension
    )
    
    message(sprintf("Phi range: [%.2f, %.2f] (data: [%.2f, %.2f])",
                    phi_range[1], phi_range[2], data_range[1], data_range[2]))
  }
  
  # Generate smooth prediction curves
  phi_seq <- seq(phi_range[1], phi_range[2], length.out = phi_n)
  pred_curves <- predict_across_phi_range(phi_seq, csp_df)
  pred_curves$AA_anacoda <- pred_curves$AA
  
  # Add codon labels - mark OPTIMAL codon (most negative dEta) with asterisk
  # This matches AnaCoDa's convention
  pred_curves <- merge(pred_curves, 
                       unique(csp_df[, c("AA", "Codon", "is_optimal")]),
                       by = c("AA", "Codon"),
                       all.x = TRUE)
  pred_curves$Codon_label <- ifelse(pred_curves$is_optimal,
                                     paste0(pred_curves$Codon, "*"),
                                     pred_curves$Codon)
  
  # Merge optimal info into observed
  obs_binned <- merge(obs_binned,
                      unique(csp_df[, c("AA", "Codon", "is_optimal")]),
                      by.x = c("AA_anacoda", "Codon"),
                      by.y = c("AA", "Codon"),
                      all.x = TRUE)
  obs_binned$is_optimal[is.na(obs_binned$is_optimal)] <- FALSE
  obs_binned$Codon_label <- ifelse(obs_binned$is_optimal,
                                    paste0(obs_binned$Codon, "*"),
                                    obs_binned$Codon)
  
  # Get unique amino acids
  aa_list <- sort(unique(pred_curves$AA_anacoda))
  
  # Create individual plots
  plot_list <- list()
  
  for (aa in aa_list) {
    # Subset data for this AA
    pred_aa <- pred_curves[pred_curves$AA_anacoda == aa, ]
    obs_aa <- obs_binned[obs_binned$AA_anacoda == aa, ]
    
    # Get codons for this AA and assign AnaCoDa colors
    codons_plain <- sort(unique(gsub("\\*$", "", pred_aa$Codon_label)))
    codons_labeled <- sort(unique(pred_aa$Codon_label))
    
    # Use AnaCoDa color scheme
    codon_colors <- sapply(codons_plain, function(c) {
      if (c %in% names(.codonColors)) .codonColors[[c]] else "gray50"
    })
    # Map colors to labeled codons (with potential asterisk)
    codon_colors_labeled <- setNames(
      sapply(codons_labeled, function(cl) {
        c_plain <- gsub("\\*$", "", cl)
        if (c_plain %in% names(.codonColors)) .codonColors[[c_plain]] else "gray50"
      }),
      codons_labeled
    )
    
    # Create plot
    plot_list[[aa]] <- plot_single_aa(aa, pred_aa, obs_aa, codon_colors_labeled)
  }
  
  # Arrange in grid
  n_plots <- length(plot_list)
  n_cols <- 5
  n_rows <- ceiling(n_plots / n_cols)
  
  # Add common axis labels
  combined_plot <- gridExtra::arrangeGrob(
    grobs = plot_list,
    ncol = n_cols,
    bottom = grid::textGrob(expression("Expression (" * log[10] * " CPM)"), 
                            gp = grid::gpar(fontsize = 12)),
    left = grid::textGrob("Proportion", rot = 90, 
                          gp = grid::gpar(fontsize = 12)),
    top = grid::textGrob("ROC Model: Predicted vs Observed Codon Usage\nLines = Model | Points = Observed ± SD | * = Optimal codon",
                         gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
  
  # Save if requested
  if (!is.null(output_file)) {
    ggsave(output_file, combined_plot, width = 26, height = 16)
    message(sprintf("Saved: %s", output_file))
  }
  
  return(list(combined = combined_plot, individual = plot_list))
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
