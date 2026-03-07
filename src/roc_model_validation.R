#' ROC Model Codon Trajectory Plot & Goodness-of-Fit - AnaCoDa Style
#' 
#' Functions to visualize predicted vs observed codon frequencies across
#' expression levels using the AnaCoDa ROC model, and to perform per-gene
#' goodness-of-fit tests comparing observed codon counts to ROC expectations.
#'
#' ROC model: P(codon_i | phi) = exp(-dM_i - dEta_i * phi) / Z
#' where phi is on linear scale (not log)
#'
#' @author Luis Javier Madrigal-Roca
#' @date 2024-12-04

# AnaCoDa COLOR SCHEME

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

# LOAD CSP PARAMETERS

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

# PREDICT CODON PROBABILITIES

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

# MAP AMINO ACIDS TO ANACODA CONVENTION

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

# CALCULATE OBSERVED FREQUENCIES BY EXPRESSION BIN

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

# MAIN PLOTTING FUNCTION - AnaCoDa Style

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

# CONVENIENCE WRAPPER

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

# LOAD PHI ESTIMATES

#' Load per-gene phi estimates from AnaCoDa gene_expression.txt
#'
#' @param phi_file Path to gene_expression.txt (AnaCoDa posterior output)
#' @return Data frame with GeneID, Phi (linear), Phi_log10
load_phi_estimates <- function(phi_file) {
  
  phi_df <- read.csv(phi_file, stringsAsFactors = FALSE)
  
  # AnaCoDa outputs: GeneID, Mean, Mean.log10 (and possibly SD, quantiles)
  if (!all(c("GeneID", "Mean") %in% names(phi_df))) {
    stop("phi_file must contain columns 'GeneID' and 'Mean'")
  }
  
  out <- data.frame(
    Gene = phi_df$GeneID,
    Phi = phi_df$Mean,          # Linear scale
    Phi_log10 = if ("Mean.log10" %in% names(phi_df)) phi_df$Mean.log10 
                else log10(phi_df$Mean),
    stringsAsFactors = FALSE
  )
  
  message(sprintf("Loaded phi estimates for %d genes", nrow(out)))
  return(out)
}

# PREDICT PER-GENE EXPECTED CODON COUNTS

#' Calculate ROC-predicted codon probabilities for a single phi value
#'
#' @param phi_linear Numeric scalar: gene's phi on linear scale
#' @param csp_df CSP data frame (from load_csp_parameters)
#' @return Data frame with AA, Codon, predicted_prob for each codon
predict_codon_probs <- function(phi_linear, csp_df) {
  
  aa_list <- split(csp_df, csp_df$AA)
  
  result <- do.call(rbind, lapply(aa_list, function(aa_df) {
    # ROC multinomial: log P(codon_i) = -dM_i - dEta_i * phi - log(Z)
    log_unnorm <- -aa_df$dM - aa_df$dEta * phi_linear
    
    # Log-sum-exp for numerical stability
    max_log <- max(log_unnorm)
    log_Z <- max_log + log(sum(exp(log_unnorm - max_log)))
    
    aa_df$predicted_prob <- exp(log_unnorm - log_Z)
    aa_df[, c("AA", "Codon", "predicted_prob")]
  }))
  
  rownames(result) <- NULL
  return(result)
}

#' Calculate expected codon counts for all genes from the ROC model
#'
#' For each gene, uses the gene-specific phi (expression level) to compute
#' the multinomial codon probabilities, then multiplies by the observed
#' amino acid totals to get expected counts per codon.
#'
#' Expected_count(codon_i, gene_g) = n_AA(gene_g) * P(codon_i | phi_g)
#'
#' @param codon_counts_long Data frame in long format with columns:
#'   Gene, Codon, AA, Count (raw codon counts per gene)
#' @param phi_df Data frame with Gene, Phi (linear scale) from load_phi_estimates
#' @param csp_df CSP parameters from load_csp_parameters
#' @param map_to_anacoda Logical; if TRUE, remap Ser AGN codons to AA = "Z"
#'   to match AnaCoDa convention (default TRUE)
#' @return Data frame with Gene, AA, Codon, Observed, Expected, residual columns
calculate_expected_counts <- function(codon_counts_long, phi_df, csp_df,
                                      map_to_anacoda = TRUE) {
  
  # Ensure we have required columns
  required_cols <- c("Gene", "Codon", "AA", "Count")
  missing <- setdiff(required_cols, names(codon_counts_long))
  if (length(missing) > 0) {
    stop(sprintf("codon_counts_long missing columns: %s", paste(missing, collapse = ", ")))
  }
  
  # Map to AnaCoDa AA convention (Z = AGC, AGT)
  if (map_to_anacoda) {
    codon_counts_long$AA_model <- codon_counts_long$AA
    codon_counts_long$AA_model[codon_counts_long$Codon %in% c("AGC", "AGT")] <- "Z"
  } else {
    codon_counts_long$AA_model <- codon_counts_long$AA
  }
  
  # Calculate AA totals per gene (using the model's AA grouping)
  aa_totals <- codon_counts_long |>
    dplyr::group_by(Gene, AA_model) |>
    dplyr::summarize(AA_total = sum(Count, na.rm = TRUE), .groups = "drop")
  
  # Merge phi values — restrict to genes with phi estimates
  genes_with_phi <- intersect(unique(codon_counts_long$Gene), phi_df$Gene)
  
  if (length(genes_with_phi) == 0) {
    stop("No overlap between gene names in codon counts and phi estimates. ",
         "Check that gene IDs match between datasets.")
  }
  
  message(sprintf("Computing expected counts for %d genes (of %d with phi estimates)",
                  length(genes_with_phi), nrow(phi_df)))
  
  # Pre-compute predicted probabilities for all unique phi values
  # (Many genes may share very similar phi, but each is unique from MCMC)
  # For efficiency, compute per gene
  phi_lookup <- setNames(phi_df$Phi, phi_df$Gene)
  
  # Get model codons (those in the CSP)
  model_codons <- unique(csp_df$Codon)
  
  # Process each gene
  result_list <- lapply(genes_with_phi, function(gene) {
    
    phi_g <- phi_lookup[gene]
    
    # Get predicted probabilities for this gene's phi
    pred <- predict_codon_probs(phi_g, csp_df)
    
    # Get observed counts for this gene
    obs_gene <- codon_counts_long[codon_counts_long$Gene == gene, ]
    
    # Merge observed counts with predictions
    merged <- merge(
      obs_gene[, c("Gene", "Codon", "AA_model", "Count")],
      pred[, c("AA", "Codon", "predicted_prob")],
      by.x = c("AA_model", "Codon"),
      by.y = c("AA", "Codon"),
      all.x = FALSE  # Only keep codons that are in the model
    )
    
    # Get AA totals for this gene
    aa_tot_gene <- aa_totals[aa_totals$Gene == gene, ]
    merged <- merge(merged, aa_tot_gene, by = c("Gene", "AA_model"))
    
    # Expected counts
    merged$Expected <- merged$AA_total * merged$predicted_prob
    merged$Observed <- merged$Count
    
    # Pearson residual: (O - E) / sqrt(E)
    merged$Pearson_residual <- ifelse(
      merged$Expected > 0,
      (merged$Observed - merged$Expected) / sqrt(merged$Expected),
      NA_real_
    )
    
    merged$Phi <- phi_g
    
    merged[, c("Gene", "AA_model", "Codon", "Observed", "Expected", 
               "AA_total", "predicted_prob", "Pearson_residual", "Phi")]
  })
  
  result <- do.call(rbind, result_list)
  rownames(result) <- NULL
  
  message(sprintf("Expected counts computed: %d gene-codon combinations", nrow(result)))
  
  return(result)
}

# GOODNESS-OF-FIT TESTS

#' Perform per-gene goodness-of-fit test (chi-squared or G-test)
#'
#' For each gene, tests whether the observed codon distribution within each
#' amino acid family differs significantly from the ROC model prediction.
#'
#' Tests are performed per amino acid family within each gene. AAs with
#' only 1 synonymous codon or zero total observations are skipped.
#'
#' @param expected_counts_df Output from calculate_expected_counts
#' @param test Type of test: "chisq" (Pearson chi-squared) or "g_test" 
#'   (log-likelihood ratio / G-test). Default: "chisq"
#' @param min_total Minimum total AA count in a gene to perform test (default: 5)
#' @param p_adjust_method Method for p.adjust (default: "BH" / Benjamini-Hochberg)
#' @return Data frame with per-gene-AA test results
gof_test_per_gene <- function(expected_counts_df, 
                               test = c("chisq", "g_test"),
                               min_total = 5,
                               p_adjust_method = "BH") {
  
  test <- match.arg(test)
  
  # Split by gene and AA
  gene_aa_groups <- split(expected_counts_df, 
                          list(expected_counts_df$Gene, expected_counts_df$AA_model),
                          drop = TRUE)
  
  results <- lapply(gene_aa_groups, function(df) {
    
    gene <- df$Gene[1]
    aa <- df$AA_model[1]
    n_codons <- nrow(df)
    aa_total <- df$AA_total[1]
    
    # Skip single-codon AAs (Met, Trp) and low-count groups
    if (n_codons < 2 || aa_total < min_total) {
      return(NULL)
    }
    
    obs <- df$Observed
    exp_counts <- df$Expected
    
    # Skip if any expected count is 0 (degenerate model prediction)
    if (any(exp_counts <= 0)) {
      return(NULL)
    }
    
    df_freedom <- n_codons - 1  # Multinomial constraint
    
    if (test == "chisq") {
      # Pearson chi-squared: sum((O - E)^2 / E)
      stat <- sum((obs - exp_counts)^2 / exp_counts)
      p_val <- pchisq(stat, df = df_freedom, lower.tail = FALSE)
    } else {
      # G-test (log-likelihood ratio): 2 * sum(O * log(O / E))
      # Handle zero observed counts safely
      nonzero <- obs > 0
      stat <- 2 * sum(obs[nonzero] * log(obs[nonzero] / exp_counts[nonzero]))
      p_val <- pchisq(stat, df = df_freedom, lower.tail = FALSE)
    }
    
    data.frame(
      Gene = gene,
      AA = aa,
      n_codons = n_codons,
      AA_total = aa_total,
      Statistic = stat,
      df = df_freedom,
      p_value = p_val,
      Phi = df$Phi[1],
      stringsAsFactors = FALSE
    )
  })
  
  results_df <- do.call(rbind, results[!sapply(results, is.null)])
  rownames(results_df) <- NULL
  
  # Multiple testing correction
  results_df$p_adj <- p.adjust(results_df$p_value, method = p_adjust_method)
  
  message(sprintf("GoF %s test: %d gene-AA tests performed, %d significant (FDR < 0.05)",
                  test,
                  nrow(results_df),
                  sum(results_df$p_adj < 0.05, na.rm = TRUE)))
  
  return(results_df)
}

#' Aggregate per-gene GoF results into a gene-level summary
#'
#' For each gene, aggregates across all amino acid families to produce
#' a single gene-level goodness-of-fit measure.
#'
#' @param gof_results Output from gof_test_per_gene
#' @param method Aggregation method: "sum" (sum chi-sq statistics and df),
#'   "fisher" (Fisher's method to combine p-values), or "both" (default "both")
#' @return Data frame with gene-level GoF summary
aggregate_gof_by_gene <- function(gof_results, method = c("both", "sum", "fisher")) {
  
  method <- match.arg(method)
  
  gene_summary <- gof_results |>
    dplyr::group_by(Gene) |>
    dplyr::summarize(
      n_aa_tested = dplyr::n(),
      total_codons_tested = sum(AA_total),
      Phi = Phi[1],
      
      # Sum of chi-sq statistics and df
      Stat_sum = sum(Statistic),
      df_sum = sum(df),
      
      # Gene-level p-value from summed statistics
      p_value_sum = pchisq(sum(Statistic), df = sum(df), lower.tail = FALSE),
      
      # Fisher's method: -2 * sum(log(p))
      Fisher_stat = -2 * sum(log(p_value)),
      Fisher_df = 2 * dplyr::n(),
      p_value_fisher = pchisq(-2 * sum(log(p_value)), 
                               df = 2 * dplyr::n(), 
                               lower.tail = FALSE),
      
      # Proportion of AA families significant at nominal level
      prop_aa_significant = mean(p_value < 0.05),
      
      # Mean Pearson residual (unsigned) — overall model deviation
      .groups = "drop"
    )
  
  # Multiple testing correction at gene level
  gene_summary$p_adj_sum <- p.adjust(gene_summary$p_value_sum, method = "BH")
  gene_summary$p_adj_fisher <- p.adjust(gene_summary$p_value_fisher, method = "BH")
  
  # Reduced chi-sq (stat / df) as a normalized measure of fit quality
  gene_summary$reduced_chisq <- gene_summary$Stat_sum / gene_summary$df_sum
  
  n_sig <- sum(gene_summary$p_adj_sum < 0.05, na.rm = TRUE)
  message(sprintf("Gene-level GoF: %d genes tested, %d (%.1f%%) significant (FDR < 0.05, sum method)",
                  nrow(gene_summary), n_sig, 100 * n_sig / nrow(gene_summary)))
  
  return(gene_summary)
}

# DIAGNOSTICS & VISUALIZATION

#' Summarize goodness-of-fit results with diagnostic plots
#'
#' @param gene_gof Gene-level GoF summary from aggregate_gof_by_gene
#' @param expected_counts Expected counts data frame from calculate_expected_counts
#' @param output_prefix Prefix for output files (optional)
#' @return List with summary statistics and plots
plot_gof_diagnostics <- function(gene_gof, expected_counts, output_prefix = NULL) {
  
  require(ggplot2)
  
  plots <- list()
  
  # 1. Reduced chi-sq distribution
  plots$reduced_chisq_hist <- ggplot(gene_gof, aes(x = reduced_chisq)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.8) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(
      title = "Reduced Chi-Squared Distribution",
      subtitle = "Red line = perfect fit (1.0); >1 = underfitting; <1 = overfitting",
      x = expression(chi[red]^2 ~ "(Statistic / df)"),
      y = "Number of genes"
    ) +
    theme_bw(base_size = 12)
  
  # 2. Reduced chi-sq vs expression (phi)
  plots$chisq_vs_phi <- ggplot(gene_gof, aes(x = Phi, y = reduced_chisq)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_smooth(method = "loess", color = "red", se = TRUE) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "blue") +
    scale_x_log10() +
    labs(
      title = "Model Fit vs Expression Level",
      subtitle = "Does the ROC model fit better for highly expressed genes?",
      x = expression(hat(phi) ~ "(linear scale, log axis)"),
      y = expression(chi[red]^2)
    ) +
    theme_bw(base_size = 12)
  
  # 3. P-value distribution (should be uniform under null)
  plots$pvalue_hist <- ggplot(gene_gof, aes(x = p_value_sum)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.8) +
    geom_hline(yintercept = nrow(gene_gof) / 50, linetype = "dashed", 
               color = "red", linewidth = 0.8) +
    labs(
      title = "P-value Distribution (Gene-level GoF)",
      subtitle = "Uniform = model fits well; left-skewed = systematic misfit",
      x = "p-value",
      y = "Count"
    ) +
    theme_bw(base_size = 12)
  
  # 4. Per-codon Pearson residuals across all genes (boxplot by codon)
  codon_residuals <- expected_counts |>
    dplyr::filter(!is.na(Pearson_residual)) |>
    dplyr::group_by(Codon) |>
    dplyr::mutate(median_resid = median(Pearson_residual, na.rm = TRUE)) |>
    dplyr::ungroup()
  
  plots$residuals_by_codon <- ggplot(codon_residuals, 
                                      aes(x = reorder(Codon, median_resid), 
                                          y = Pearson_residual)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, fill = "lightblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    coord_flip() +
    labs(
      title = "Pearson Residuals by Codon",
      subtitle = "Positive = model underestimates usage; Negative = model overestimates",
      x = "Codon",
      y = "Pearson residual (O - E) / sqrt(E)"
    ) +
    theme_bw(base_size = 10)
  
  # Save if requested
  if (!is.null(output_prefix)) {
    combined <- gridExtra::arrangeGrob(
      plots$reduced_chisq_hist, plots$chisq_vs_phi,
      plots$pvalue_hist, plots$residuals_by_codon,
      ncol = 2,
      top = grid::textGrob("ROC Model Goodness-of-Fit Diagnostics",
                           gp = grid::gpar(fontsize = 16, fontface = "bold"))
    )
    ggsave(paste0(output_prefix, "_gof_diagnostics.pdf"), combined, 
           width = 16, height = 12)
    message(sprintf("Saved: %s_gof_diagnostics.pdf", output_prefix))
  }
  
  # Print summary
  cat("\n=== ROC Model Goodness-of-Fit Summary ===\n")
  cat(sprintf("Genes tested: %d\n", nrow(gene_gof)))
  cat(sprintf("Median reduced chi-sq: %.3f\n", median(gene_gof$reduced_chisq, na.rm = TRUE)))
  cat(sprintf("Mean reduced chi-sq:   %.3f\n", mean(gene_gof$reduced_chisq, na.rm = TRUE)))
  cat(sprintf("Genes with good fit (reduced chi-sq < 2): %d (%.1f%%)\n",
              sum(gene_gof$reduced_chisq < 2, na.rm = TRUE),
              100 * mean(gene_gof$reduced_chisq < 2, na.rm = TRUE)))
  cat(sprintf("Significant misfit (FDR < 0.05): %d (%.1f%%)\n",
              sum(gene_gof$p_adj_sum < 0.05, na.rm = TRUE),
              100 * mean(gene_gof$p_adj_sum < 0.05, na.rm = TRUE)))
  
  # Spearman correlation between phi and reduced chi-sq
  cor_test <- cor.test(gene_gof$Phi, gene_gof$reduced_chisq, method = "spearman",
                       exact = F)
  cat(sprintf("Spearman rho (phi vs reduced chi-sq): %.3f (p = %.2e)\n",
              cor_test$estimate, cor_test$p.value))
  
  return(list(plots = plots, summary = gene_gof))
}

# COMPLETE GoF ANALYSIS WRAPPER

#' Run complete goodness-of-fit analysis for the ROC model
#'
#' Loads CSP parameters and phi estimates, computes expected codon counts
#' for each gene, performs per-gene chi-squared tests, aggregates results,
#' and generates diagnostic plots.
#'
#' @param mutation_file Path to Cluster_X_Mutation.csv
#' @param selection_file Path to Cluster_X_Selection.csv
#' @param phi_file Path to gene_expression.txt (AnaCoDa phi estimates)
#' @param codon_counts_long Data frame with Gene, Codon, AA, Count columns
#'   (long format, raw codon counts per gene — NOT frequencies)
#' @param test GoF test type: "chisq" or "g_test" (default: "chisq")
#' @param min_aa_total Minimum amino acid family size per gene to test (default: 5)
#' @param output_prefix Path prefix for saving diagnostics (optional)
#' @return List with expected_counts, gof_per_aa, gof_per_gene, and plots
run_gof_analysis <- function(mutation_file, selection_file, phi_file,
                              codon_counts_long,
                              test = "chisq",
                              min_aa_total = 5,
                              output_prefix = NULL) {
  
  # 1. Load model parameters
  csp <- load_csp_parameters(mutation_file, selection_file)
  phi <- load_phi_estimates(phi_file)
  
  # 2. Compute expected counts
  expected <- calculate_expected_counts(codon_counts_long, phi, csp,
                                         map_to_anacoda = TRUE)
  
  # 3. Per-gene-AA goodness-of-fit test
  gof_aa <- gof_test_per_gene(expected, test = test, min_total = min_aa_total)
  
  # 4. Aggregate to gene level
  gof_gene <- aggregate_gof_by_gene(gof_aa)
  
  # 5. Diagnostics
  diag <- plot_gof_diagnostics(gof_gene, expected, output_prefix = output_prefix)
  
  return(list(
    expected_counts = expected,
    gof_per_aa = gof_aa,
    gof_per_gene = gof_gene,
    diagnostics = diag,
    csp = csp,
    phi = phi
  ))
}

# =============================================================================
# INDEPENDENT MULTINOMIAL VALIDATION OF ROC PREFERRED CODONS
# =============================================================================
#
# The ROC model separates mutational bias (dM, fixed from introns) from
# selection (dEta) to identify preferred codons. A naive multinomial on
# expression alone confounds the two — it will pick T-ending (AT-biased)
# codons because Mimulus has strong AT mutational pressure.
#
# The correct validation is: "Does the ROC-predicted preferred codon's
# usage increase significantly with expression?"  This is tested via:
#   1. Per-AA binomial GLMs of preferred-codon proportion ~ expression
#   2. Full multinomial models (fixing the VGAM naming bug for 2-codon AAs)
#      to inspect delta_P for the ROC-preferred codon
# =============================================================================


#' Targeted validation: does ROC preferred codon usage increase with expression?
#'
#' For each amino acid family, computes per-gene frequency of the ROC-preferred
#' codon and fits a binomial GLM:
#'   cbind(preferred_count, other_count) ~ expression + log(CDS_length)
#' A positive expression coefficient confirms the ROC model's prediction.
#'
#' @param codon_counts_long Data frame with Gene, Codon, AA, Count columns
#' @param expression_df Data frame with Gene, Exp_log10 columns
#' @param length_df Optional data frame with Gene, CDS_length_nt columns
#' @param roc_preferred_df ROC preferred codons (columns: Family, Preferred_Codons)
#' @param min_aa_total Minimum AA total per gene to include (default: 5)
#' @return Data frame with per-AA test results
validate_roc_preferred_binomial <- function(
    codon_counts_long,
    expression_df,
    length_df = NULL,
    roc_preferred_df,
    min_aa_total = 5
) {

  # --- Parse ROC preferred codons ---
  roc <- roc_preferred_df
  if ("Preferred_Codons" %in% names(roc) && "Family" %in% names(roc)) {
    roc_map <- setNames(roc$Preferred_Codons, roc$Family)
  } else if ("Codon" %in% names(roc) && "aa" %in% names(roc)) {
    roc_map <- setNames(roc$Codon, roc$aa)
  } else {
    stop("Cannot parse ROC preferred codons format.")
  }

  # --- Merge data ---
  dat <- merge(
    codon_counts_long[, c("Gene", "Codon", "AA", "Count")],
    expression_df[, c("Gene", "Exp_log10")],
    by = "Gene"
  )
  has_length <- !is.null(length_df)
  if (has_length) {
    dat <- merge(dat, length_df[, c("Gene", "CDS_length_nt")], by = "Gene")
    dat$log_CDS_length <- log10(dat$CDS_length_nt)
  }

  # Map to AnaCoDa AA convention
  dat <- map_aa_to_anacoda(dat)

  # --- Per-gene preferred vs non-preferred counts ---
  aa_families <- sort(intersect(names(roc_map), unique(dat$AA_anacoda)))

  results_list <- list()

  for (aa in aa_families) {
    pref_codon <- roc_map[aa]
    aa_dat <- dat[dat$AA_anacoda == aa, ]

    # Compute per-gene: preferred count, total count
    gene_summary <- aa_dat |>
      dplyr::group_by(Gene) |>
      dplyr::summarize(
        pref_count = sum(Count[Codon == pref_codon], na.rm = TRUE),
        aa_total = sum(Count, na.rm = TRUE),
        Exp_log10 = Exp_log10[1],
        .groups = "drop"
      )

    if (has_length) {
      len_vals <- unique(aa_dat[, c("Gene", "log_CDS_length")])
      gene_summary <- merge(gene_summary, len_vals, by = "Gene")
    }

    gene_summary$other_count <- gene_summary$aa_total - gene_summary$pref_count
    gene_summary <- gene_summary[gene_summary$aa_total >= min_aa_total, ]

    if (nrow(gene_summary) < 50) next

    # Fit binomial GLM
    fit <- tryCatch({
      if (has_length) {
        glm(cbind(pref_count, other_count) ~ Exp_log10 + log_CDS_length,
            family = binomial(link = "logit"), data = gene_summary)
      } else {
        glm(cbind(pref_count, other_count) ~ Exp_log10,
            family = binomial(link = "logit"), data = gene_summary)
      }
    }, error = function(e) NULL)

    if (is.null(fit)) next

    summ <- summary(fit)
    exp_row <- which(rownames(summ$coefficients) == "Exp_log10")
    if (length(exp_row) == 0) next

    exp_beta <- summ$coefficients[exp_row, 1]
    exp_se   <- summ$coefficients[exp_row, 2]
    exp_z    <- summ$coefficients[exp_row, 3]
    exp_p    <- summ$coefficients[exp_row, 4]

    # Mean preferred frequency at low vs high expression
    q10 <- quantile(gene_summary$Exp_log10, 0.10, na.rm = TRUE)
    q90 <- quantile(gene_summary$Exp_log10, 0.90, na.rm = TRUE)

    low_exp  <- gene_summary[gene_summary$Exp_log10 <= q10, ]
    high_exp <- gene_summary[gene_summary$Exp_log10 >= q90, ]

    freq_low  <- sum(low_exp$pref_count) / sum(low_exp$aa_total)
    freq_high <- sum(high_exp$pref_count) / sum(high_exp$aa_total)

    results_list[[length(results_list) + 1]] <- data.frame(
      AA = aa,
      ROC_Preferred = pref_codon,
      n_genes = nrow(gene_summary),
      Exp_beta = exp_beta,
      Exp_SE = exp_se,
      Exp_z = exp_z,
      Exp_pvalue = exp_p,
      Freq_low_exp = freq_low,
      Freq_high_exp = freq_high,
      Delta_freq = freq_high - freq_low,
      Direction = ifelse(exp_beta > 0, "Increases", "Decreases"),
      stringsAsFactors = FALSE
    )
  }

  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  results_df$p_adj <- p.adjust(results_df$Exp_pvalue, method = "BH")
  results_df$Significant <- results_df$p_adj < 0.05
  results_df$Validated <- results_df$Significant & results_df$Exp_beta > 0

  n_validated <- sum(results_df$Validated, na.rm = TRUE)
  n_tested <- nrow(results_df)
  n_increases <- sum(results_df$Exp_beta > 0, na.rm = TRUE)

  message(sprintf("\n=== Binomial GLM Validation of ROC Preferred Codons ==="))
  message(sprintf("  %d / %d AA families: preferred codon frequency increases with expression",
                  n_increases, n_tested))
  message(sprintf("  %d / %d statistically significant (FDR < 0.05)",
                  n_validated, n_tested))

  return(results_df)
}


#' Fit full multinomial models per AA family
#'
#' For each amino acid family, fits a VGAM multinomial logistic regression
#' with expression and gene length as predictors. Reports delta_P for every
#' codon, including the ROC-preferred one.
#'
#' @param codon_counts_long Data frame with Gene, Codon, AA, Count columns
#' @param expression_df Data frame with Gene, Exp_log10 columns
#' @param length_df Optional data frame with Gene, CDS_length_nt columns
#' @param csp_df Optional CSP parameters for dEta correlation
#' @param min_aa_total Minimum AA occurrence per gene (default: 5)
#' @return List with per_codon results and model_details
fit_multinomial_per_aa <- function(
    codon_counts_long,
    expression_df,
    length_df = NULL,
    csp_df = NULL,
    min_aa_total = 5
) {

  require(VGAM)

  # --- 1. Merge and prepare data ---
  dat <- merge(
    codon_counts_long[, c("Gene", "Codon", "AA", "Count")],
    expression_df[, c("Gene", "Exp_log10")],
    by = "Gene"
  )

  has_length <- !is.null(length_df)
  if (has_length) {
    dat <- merge(dat, length_df[, c("Gene", "CDS_length_nt")], by = "Gene")
    dat$log_CDS_length <- log10(dat$CDS_length_nt)
  }

  dat <- map_aa_to_anacoda(dat)

  # --- 2. Identify multi-codon AA families ---
  aa_codons <- tapply(dat$Codon, dat$AA_anacoda, function(x) sort(unique(x)))
  multi_aas <- sort(names(aa_codons[sapply(aa_codons, length) >= 2]))

  exp_q10 <- quantile(dat$Exp_log10, 0.10, na.rm = TRUE)
  exp_q90 <- quantile(dat$Exp_log10, 0.90, na.rm = TRUE)
  mean_len <- if (has_length) mean(dat$log_CDS_length, na.rm = TRUE) else NULL

  message(sprintf("Fitting multinomial models for %d amino acid families", length(multi_aas)))
  message(sprintf("Prediction contrast: expression %.2f vs %.2f (10th vs 90th percentile)",
                  exp_q10, exp_q90))

  all_codon_results <- list()
  model_objects <- list()

  for (aa in multi_aas) {

    codons <- aa_codons[[aa]]
    k <- length(codons)
    aa_dat <- dat[dat$AA_anacoda == aa, ]

    wide <- aa_dat |>
      dplyr::select(Gene, Codon, Count, Exp_log10,
                    dplyr::any_of("log_CDS_length")) |>
      tidyr::pivot_wider(names_from = Codon, values_from = Count, values_fill = 0)

    wide$aa_total <- rowSums(wide[, codons, drop = FALSE])
    wide <- wide[wide$aa_total >= min_aa_total, ]

    if (nrow(wide) < 50) {
      message(sprintf("  %s: skipped (%d usable genes, need >= 50)", aa, nrow(wide)))
      next
    }

    Y <- as.matrix(wide[, codons, drop = FALSE])

    fit <- tryCatch({
      if (has_length) {
        VGAM::vglm(Y ~ Exp_log10 + log_CDS_length,
                    family = VGAM::multinomial(refLevel = k),
                    data = wide)
      } else {
        VGAM::vglm(Y ~ Exp_log10,
                    family = VGAM::multinomial(refLevel = k),
                    data = wide)
      }
    }, error = function(e) {
      message(sprintf("  %s: model failed - %s", aa, conditionMessage(e)))
      NULL
    })

    if (is.null(fit)) next

    # --- Extract expression coefficients ---
    # VGAM naming: for M>1 logits, "Exp_log10:1", "Exp_log10:2", ...
    #              for M=1 logit,  "Exp_log10" (no colon-suffix)
    cc <- coef(fit)
    M <- k - 1

    # Match with or without colon-suffix
    exp_idx <- grep("^Exp_log10(:|$)", names(cc))
    if (length(exp_idx) != M) {
      message(sprintf("  %s: unexpected coefficient structure (%d Exp coefs, expected %d)",
                      aa, length(exp_idx), M))
      next
    }
    exp_betas <- cc[exp_idx]

    all_exp_betas <- c(unname(exp_betas), 0)
    names(all_exp_betas) <- codons

    # --- Predicted probabilities at low and high expression ---
    newdata_low <- data.frame(Exp_log10 = exp_q10)
    newdata_high <- data.frame(Exp_log10 = exp_q90)
    if (has_length) {
      newdata_low$log_CDS_length <- mean_len
      newdata_high$log_CDS_length <- mean_len
    }

    p_low <- tryCatch(
      as.numeric(predict(fit, newdata = newdata_low, type = "response")[1, ]),
      error = function(e) rep(NA_real_, k)
    )
    p_high <- tryCatch(
      as.numeric(predict(fit, newdata = newdata_high, type = "response")[1, ]),
      error = function(e) rep(NA_real_, k)
    )
    names(p_low) <- names(p_high) <- codons

    delta_p <- p_high - p_low

    valid_delta <- delta_p[!is.na(delta_p) & !is.nan(delta_p)]
    if (length(valid_delta) == 0) {
      message(sprintf("  %s: skipped (could not compute probability predictions)", aa))
      next
    }
    multinom_preferred <- names(which.max(valid_delta))

    # --- Wald p-values ---
    p_values <- rep(NA_real_, M)
    summ <- tryCatch(summary(fit), error = function(e) NULL)
    if (!is.null(summ)) {
      coef_table <- coef(summ)
      exp_row_idx <- grep("^Exp_log10(:|$)", rownames(coef_table))
      for (j in seq_along(exp_row_idx)) {
        if (j <= M) p_values[j] <- coef_table[exp_row_idx[j], 4]
      }
    }
    p_values <- c(p_values, NA_real_)

    for (i in seq_along(codons)) {
      all_codon_results[[length(all_codon_results) + 1]] <- data.frame(
        AA = aa,
        Codon = codons[i],
        n_codons_family = k,
        n_genes = nrow(wide),
        Exp_beta = unname(all_exp_betas[i]),
        Exp_pvalue = p_values[i],
        Prob_low_exp = unname(p_low[i]),
        Prob_high_exp = unname(p_high[i]),
        Delta_prob = unname(delta_p[i]),
        is_multinom_preferred = (codons[i] == multinom_preferred),
        is_reference = (i == k),
        stringsAsFactors = FALSE
      )
    }

    model_objects[[aa]] <- list(
      fit = fit, codons = codons, preferred = multinom_preferred,
      prob_low = p_low, prob_high = p_high, delta_p = delta_p
    )

    message(sprintf("  %s [%d codons, %d genes]: multinom preferred = %s (Delta_P = %+.3f)",
                    aa, k, nrow(wide), multinom_preferred, max(valid_delta)))
  }

  results_df <- do.call(rbind, all_codon_results)
  rownames(results_df) <- NULL

  # Quantitative correlation: dEta vs delta_P
  deta_cor <- NULL
  if (!is.null(csp_df) && nrow(results_df) > 0) {
    csp_merge <- csp_df[, c("AA", "Codon", "dEta")]
    merged <- merge(results_df[, c("AA", "Codon", "Delta_prob")],
                    csp_merge, by = c("AA", "Codon"))
    merged <- merged[!is.na(merged$Delta_prob), ]

    if (nrow(merged) >= 5) {
      cor_result <- cor.test(merged$dEta, merged$Delta_prob,
                             method = "spearman", exact = FALSE)
      deta_cor <- list(
        rho = cor_result$estimate,
        p_value = cor_result$p.value,
        n = nrow(merged),
        data = merged
      )
      message(sprintf("\n=== dEta vs Delta_P Correlation ==="))
      message(sprintf("  Spearman rho = %.3f (p = %.2e, n = %d codons)",
                      cor_result$estimate, cor_result$p.value, nrow(merged)))
      message("  Expected: negative (more negative dEta -> larger Delta_P)")
    }
  }

  return(list(
    per_codon = results_df,
    model_details = model_objects,
    dEta_correlation = deta_cor
  ))
}


#' Visualize multinomial validation results
#'
#' Creates a multi-panel figure: (1) per-AA barplot of Delta_P highlighting
#' the ROC-preferred codon's behavior, (2) binomial validation forest plot,
#' and (3) dEta vs Delta_P scatter.
#'
#' @param binomial_results Output from validate_roc_preferred_binomial
#' @param multinom_results Output from fit_multinomial_per_aa
#' @param roc_preferred_df ROC preferred codons
#' @param output_file Path to save PDF (optional)
#' @return List of ggplot objects
plot_multinomial_validation <- function(binomial_results,
                                        multinom_results = NULL,
                                        roc_preferred_df = NULL,
                                        output_file = NULL) {

  require(ggplot2)
  require(gridExtra)

  plots <- list()

  # --- Panel 1: Forest plot of binomial validation ---
  binom_df <- binomial_results
  binom_df$AA_label <- paste0(binom_df$AA, " (", binom_df$ROC_Preferred, ")")
  binom_df$AA_label <- factor(binom_df$AA_label,
                               levels = binom_df$AA_label[order(binom_df$Exp_beta)])

  plots$binomial_forest <- ggplot(binom_df,
                                   aes(x = Exp_beta, y = AA_label,
                                       color = Validated)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(aes(xmin = Exp_beta - 1.96 * Exp_SE,
                        xmax = Exp_beta + 1.96 * Exp_SE),
                    size = 0.5) +
    scale_color_manual(values = c("TRUE" = "forestgreen", "FALSE" = "firebrick"),
                       labels = c("TRUE" = "Validated (FDR < 0.05)",
                                  "FALSE" = "Not validated"),
                       name = NULL) +
    labs(
      title = "ROC Preferred Codon Validation: Expression Effect",
      subtitle = paste0("Binomial GLM: preferred_count / AA_total ~ expression + log(CDS_length)\n",
                        "Positive = ROC-preferred codon usage increases with expression"),
      x = "Expression coefficient (logit scale)",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 13, face = "bold"))

  # --- Panel 2: Delta_P barplot (multinomial) with ROC preferred highlighted ---
  if (!is.null(multinom_results) && !is.null(roc_preferred_df)) {
    df <- multinom_results$per_codon

    roc <- roc_preferred_df
    if ("Preferred_Codons" %in% names(roc) && "Family" %in% names(roc)) {
      roc_map <- setNames(roc$Preferred_Codons, roc$Family)
    } else {
      roc_map <- setNames(roc$Codon, roc$aa)
    }

    df$is_ROC_preferred <- mapply(function(aa, codon) {
      if (aa %in% names(roc_map)) codon == roc_map[aa] else FALSE
    }, df$AA, df$Codon)

    df$Category <- "Other codon"
    df$Category[df$is_ROC_preferred & df$Delta_prob > 0] <- "ROC preferred (increases)"
    df$Category[df$is_ROC_preferred & df$Delta_prob <= 0] <- "ROC preferred (decreases)"
    df$Category[df$is_multinom_preferred] <- paste0(
      ifelse(df$is_ROC_preferred[df$is_multinom_preferred],
             "ROC preferred (increases)", "Highest Delta_P (not ROC)")
    )
    # Re-assign cleanly
    df$Category <- "Other codon"
    df$Category[df$is_multinom_preferred & !df$is_ROC_preferred] <- "Highest Delta_P"
    df$Category[df$is_ROC_preferred & df$Delta_prob > 0] <- "ROC preferred (increases)"
    df$Category[df$is_ROC_preferred & df$Delta_prob <= 0] <- "ROC preferred (decreases)"

    color_vals <- c("ROC preferred (increases)" = "forestgreen",
                    "ROC preferred (decreases)" = "firebrick",
                    "Highest Delta_P" = "steelblue",
                    "Other codon" = "gray70")

    plots$delta_prob <- ggplot(df, aes(x = reorder(Codon, Delta_prob),
                                       y = Delta_prob, fill = Category)) +
      geom_col(width = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black",
                 linewidth = 0.3) +
      facet_wrap(~ AA, scales = "free", ncol = 5) +
      coord_flip() +
      scale_fill_manual(values = color_vals) +
      labs(
        title = "Multinomial Model: Codon Probability Change with Expression",
        subtitle = paste0(
          "Delta_P = P(codon | high exp.) - P(codon | low exp.) | ",
          "Green = ROC preferred codon increases with expression"
        ),
        x = NULL,
        y = expression(Delta * "P (high exp. " - " low exp.)"),
        fill = NULL
      ) +
      theme_bw(base_size = 10) +
      theme(
        legend.position = "bottom",
        strip.text = element_text(face = "bold"),
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 9)
      )
  }

  # --- Panel 3: dEta vs Delta_P scatter ---
  if (!is.null(multinom_results$dEta_correlation) &&
      nrow(multinom_results$dEta_correlation$data) >= 5) {
    deta_dat <- multinom_results$dEta_correlation$data
    deta_rho <- multinom_results$dEta_correlation$rho
    deta_p   <- multinom_results$dEta_correlation$p_value

    plots$deta_scatter <- ggplot(deta_dat, aes(x = dEta, y = Delta_prob)) +
      geom_point(alpha = 0.6, size = 2, color = "steelblue") +
      geom_smooth(method = "lm", color = "red", se = TRUE, linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      annotate("text", x = Inf, y = Inf,
               label = sprintf("rho = %.3f\np = %.2e",
                               deta_rho, deta_p),
               hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
      labs(
        title = expression("ROC Selection (" * Delta * eta * ") vs Multinomial " * Delta * "P"),
        subtitle = "Expected: negative (selected-for codons increase with expression)",
        x = expression("ROC Selection Coefficient (" * Delta * eta * ")"),
        y = expression(Delta * "P (high " - " low expression)")
      ) +
      theme_bw(base_size = 11)
  }

  # --- Combine and save ---
  if (!is.null(output_file)) {
    grob_list <- lapply(plots, function(p) {
      if (inherits(p, "ggplot")) ggplotGrob(p) else p
    })

    n_panels <- length(grob_list)
    if (n_panels == 1) {
      combined <- grob_list[[1]]
      plot_height <- 10
    } else if (n_panels == 2) {
      combined <- gridExtra::arrangeGrob(grobs = grob_list, ncol = 1,
                                          heights = c(1, 1.5))
      plot_height <- 18
    } else {
      combined <- gridExtra::arrangeGrob(grobs = grob_list, ncol = 1,
                                          heights = c(1, 1.5, 0.8))
      plot_height <- 24
    }

    ggsave(output_file, combined, width = 18, height = plot_height)
    message(sprintf("Saved: %s", output_file))
  }

  return(plots)
}


#' Run complete independent multinomial validation
#'
#' Fits per-AA binomial GLMs (targeted test) and full multinomial models,
#' compares with ROC preferred codons, and generates diagnostic plots.
#'
#' @param codon_counts_long Data frame with Gene, Codon, AA, Count columns
#' @param expression_df Data frame with Gene, Exp_log10 columns
#' @param length_df Optional data frame with Gene, CDS_length_nt columns
#' @param roc_preferred_df ROC preferred codons (columns: Family, Preferred_Codons)
#' @param csp_df Optional CSP parameters for dEta correlation analysis
#' @param min_aa_total Minimum AA count per gene (default: 5)
#' @param output_prefix Path prefix for saving outputs (optional)
#' @return List with binomial_validation, multinomial_results, plots
run_multinomial_validation <- function(
    codon_counts_long,
    expression_df,
    length_df = NULL,
    roc_preferred_df = NULL,
    csp_df = NULL,
    min_aa_total = 5,
    output_prefix = NULL
) {

  # 1. Targeted binomial validation
  binomial_results <- validate_roc_preferred_binomial(
    codon_counts_long = codon_counts_long,
    expression_df = expression_df,
    length_df = length_df,
    roc_preferred_df = roc_preferred_df,
    min_aa_total = min_aa_total
  )

  # 2. Full multinomial models
  multinom_results <- fit_multinomial_per_aa(
    codon_counts_long = codon_counts_long,
    expression_df = expression_df,
    length_df = length_df,
    csp_df = csp_df,
    min_aa_total = min_aa_total
  )

  # 3. Build comparison table: merge binomial + multinomial preferred
  comparison <- binomial_results[, c("AA", "ROC_Preferred", "Exp_beta",
                                      "Direction", "Validated",
                                      "Delta_freq")]

  # Add multinomial preferred codon and ROC preferred's delta_P
  if (nrow(multinom_results$per_codon) > 0) {
    multinom_pref <- multinom_results$per_codon[
      multinom_results$per_codon$is_multinom_preferred,
      c("AA", "Codon", "Delta_prob")]
    names(multinom_pref) <- c("AA", "Multinom_Preferred", "Multinom_Delta_P")
    comparison <- merge(comparison, multinom_pref, by = "AA", all.x = TRUE)

    # ROC preferred codon's Delta_P from the multinomial
    roc_in_multinom <- merge(
      binomial_results[, c("AA", "ROC_Preferred")],
      multinom_results$per_codon[, c("AA", "Codon", "Delta_prob")],
      by.x = c("AA", "ROC_Preferred"), by.y = c("AA", "Codon"),
      all.x = TRUE
    )
    names(roc_in_multinom)[3] <- "ROC_Delta_P"
    comparison <- merge(comparison, roc_in_multinom[, c("AA", "ROC_Delta_P")],
                        by = "AA", all.x = TRUE)

    comparison$ROC_increases_in_multinom <- !is.na(comparison$ROC_Delta_P) &
                                             comparison$ROC_Delta_P > 0
  }

  # 4. Generate plots
  plot_file <- if (!is.null(output_prefix)) paste0(output_prefix, ".pdf") else NULL

  plots <- plot_multinomial_validation(
    binomial_results = binomial_results,
    multinom_results = multinom_results,
    roc_preferred_df = roc_preferred_df,
    output_file = plot_file
  )

  # 5. Save comparison table
  if (!is.null(output_prefix)) {
    write.csv(comparison, paste0(output_prefix, "_concordance.csv"), row.names = FALSE)
    message(sprintf("Saved: %s_concordance.csv", output_prefix))
  }

  # 6. Print summary
  cat("\n=== ROC Preferred Codon Validation Summary ===\n")
  n_tested <- nrow(binomial_results)
  n_increases <- sum(binomial_results$Exp_beta > 0, na.rm = TRUE)
  n_validated <- sum(binomial_results$Validated, na.rm = TRUE)

  cat(sprintf("AA families tested: %d\n", n_tested))
  cat(sprintf("ROC preferred increases with expression: %d / %d (%.1f%%)\n",
              n_increases, n_tested, 100 * n_increases / n_tested))
  cat(sprintf("Statistically significant (FDR < 0.05): %d / %d (%.1f%%)\n",
              n_validated, n_tested, 100 * n_validated / n_tested))
  cat("\nDetailed comparison:\n")
  print(comparison[order(comparison$AA), ])

  return(list(
    binomial_validation = binomial_results,
    multinomial_results = multinom_results,
    comparison = comparison,
    plots = plots
  ))
}
