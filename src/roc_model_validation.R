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
