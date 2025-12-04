#' ROC Model Validation Functions
#' 
#' Functions to validate AnaCoDa ROC model predictions using empirical expression data.
#' The ROC model predicts codon probabilities as:
#'   P(codon_i | phi) = exp(-dM_i - dEta_i * phi) / Z
#' where Z is the partition function (sum over all codons in the amino acid family).
#'
#' @author Luis Javier Madrigal-Roca
#' @date 2024-12-04

# ==============================================================================
# LOADING CSP PARAMETERS
# ==============================================================================

#' Load CSP (Codon-Specific Parameters) from AnaCoDa output
#'
#' @param mutation_file Path to Cluster_X_Mutation.csv file
#' @param selection_file Path to Cluster_X_Selection.csv file
#' @return Data frame with columns: AA, Codon, dM, dEta
#' @export
load_csp_parameters <- function(mutation_file, selection_file) {
  
  # Load mutation parameters (dM)
  dM_params <- read.csv(mutation_file, stringsAsFactors = FALSE)
  dM_params <- dM_params[, c("AA", "Codon", "Mean")]
  names(dM_params)[3] <- "dM"
  
  # Load selection parameters (dEta)
  dEta_params <- read.csv(selection_file, stringsAsFactors = FALSE)
  dEta_params <- dEta_params[, c("AA", "Codon", "Mean")]
  names(dEta_params)[3] <- "dEta"
  

  # Merge
  csp_params <- merge(dM_params, dEta_params, by = c("AA", "Codon"))
  
  # Handle the "Z" amino acid (Serine AGN codons in AnaCoDa convention)
  # AnaCoDa splits Serine into S (TCN) and Z (AGC, AGT)
  # For biological interpretation, we keep them separate as they have different

  # mutational neighborhoods
  
  message(sprintf("Loaded CSP parameters for %d codons across %d amino acid families",
                  nrow(csp_params), length(unique(csp_params$AA))))
  
  return(csp_params)
}

# ==============================================================================
# CODON PROBABILITY PREDICTION
# ==============================================================================

#' Predict codon probabilities for a single phi value
#'
#' Uses the ROC multinomial model:
#'   P(codon_i | phi) = exp(-dM_i - dEta_i * phi) / Z
#'
#' @param phi Expression level (typically log10 scale)
#' @param csp_df Data frame with columns: AA, Codon, dM, dEta
#' @return Data frame with columns: AA, Codon, predicted_prob
#' @export
predict_codon_probs_single <- function(phi, csp_df) {
  
  # Calculate log-unnormalized probabilities
  csp_df$log_unnorm <- -csp_df$dM - csp_df$dEta * phi
  
  # Split by amino acid and normalize within each family
  aa_list <- split(csp_df, csp_df$AA)
  
  result_list <- lapply(aa_list, function(aa_df) {
    # Log-sum-exp for numerical stability
    max_log <- max(aa_df$log_unnorm)
    log_Z <- max_log + log(sum(exp(aa_df$log_unnorm - max_log)))
    
    aa_df$predicted_prob <- exp(aa_df$log_unnorm - log_Z)
    aa_df[, c("AA", "Codon", "predicted_prob")]
  })
  
  do.call(rbind, result_list)
}

#' Predict codon probabilities for multiple genes
#'
#' @param gene_expr_df Data frame with columns: Gene, Exp_log10 (or phi column)
#' @param csp_df Data frame with CSP parameters
#' @param phi_col Name of the column containing phi values (default: "Exp_log10")
#' @param verbose Print progress (default: TRUE)
#' @return Data frame with Gene, AA, Codon, predicted_prob, Phi_empirical
#' @export
predict_codon_probs_batch <- function(gene_expr_df, csp_df, phi_col = "Exp_log10", 
                                       verbose = TRUE) {
  
  n_genes <- nrow(gene_expr_df)
  
  if (verbose) {
    message(sprintf("Predicting codon probabilities for %d genes...", n_genes))
  }
  
  # Pre-allocate list
  predictions_list <- vector("list", n_genes)
  
  # Progress tracking
  progress_step <- max(1, floor(n_genes / 10))
  
  for (i in seq_len(n_genes)) {
    gene_id <- gene_expr_df$Gene[i]
    phi <- gene_expr_df[[phi_col]][i]
    
    pred <- predict_codon_probs_single(phi, csp_df)
    pred$Gene <- gene_id
    pred$Phi_empirical <- phi
    
    predictions_list[[i]] <- pred
    
    if (verbose && i %% progress_step == 0) {
      message(sprintf("  Progress: %d/%d (%.0f%%)", i, n_genes, 100 * i / n_genes))
    }
  }
  
  result <- do.call(rbind, predictions_list)
  rownames(result) <- NULL
  
  if (verbose) {
    message("Done!")
  }
  
  return(result)
}

# ==============================================================================
# CODON COUNTING FROM SEQUENCES
# ==============================================================================

#' Count codons in a single DNA sequence
#'
#' @param seq DNAString or character string
#' @return Named vector of codon counts
#' @export
count_codons_in_sequence <- function(seq) {
  seq_str <- as.character(seq)
  n <- nchar(seq_str)
  
  # Need at least 3 nucleotides
  if (n < 3) return(NULL)
  
  # Extract codons (excluding incomplete final codon)
  n_codons <- floor(n / 3)
  
  starts <- seq(1, n_codons * 3, by = 3)
  ends <- seq(3, n_codons * 3, by = 3)
  codons <- substring(seq_str, starts, ends)
  
  # Remove codons with N
  codons <- codons[!grepl("N", codons)]
  
  if (length(codons) == 0) return(NULL)
  
  table(codons)
}

#' Calculate observed codon frequencies for all genes
#'
#' @param cds_seqs DNAStringSet of CDS sequences
#' @param genetic_code Named vector mapping codons to amino acids (default: standard)
#' @param verbose Print progress (default: TRUE)
#' @return Data frame with Gene, Codon, AA, Count, AA_total, Observed_freq
#' @export
calculate_observed_codon_frequencies <- function(cds_seqs, 
                                                   genetic_code = Biostrings::GENETIC_CODE,
                                                   verbose = TRUE) {
  
  gene_names <- names(cds_seqs)
  # Clean gene names (remove transcript suffix like .1, .2)
  gene_names_clean <- sub("\\.[0-9]+$", "", gene_names)
  
  n_genes <- length(cds_seqs)
  
  if (verbose) {
    message(sprintf("Counting codons for %d genes...", n_genes))
  }
  
  codon_freq_list <- vector("list", n_genes)
  progress_step <- max(1, floor(n_genes / 10))
  
  for (i in seq_len(n_genes)) {
    gene_id <- gene_names_clean[i]
    counts <- count_codons_in_sequence(cds_seqs[[i]])
    
    if (!is.null(counts)) {
      df <- data.frame(
        Gene = gene_id,
        Codon = names(counts),
        Count = as.numeric(counts),
        stringsAsFactors = FALSE
      )
      
      # Map codons to amino acids
      df$AA <- genetic_code[df$Codon]
      
      # Remove stop codons and unknown
      df <- df[!is.na(df$AA) & df$AA != "*", ]
      
      if (nrow(df) > 0) {
        # Calculate within-AA frequencies
        aa_totals <- tapply(df$Count, df$AA, sum)
        df$AA_total <- aa_totals[df$AA]
        df$Observed_freq <- df$Count / df$AA_total
        
        codon_freq_list[[i]] <- df
      }
    }
    
    if (verbose && i %% progress_step == 0) {
      message(sprintf("  Progress: %d/%d (%.0f%%)", i, n_genes, 100 * i / n_genes))
    }
  }
  
  result <- do.call(rbind, codon_freq_list)
  rownames(result) <- NULL
  
  if (verbose) {
    message(sprintf("Done! Calculated frequencies for %d gene-codon combinations", nrow(result)))
  }
  
  return(result)
}

# ==============================================================================
# MODEL VALIDATION
# ==============================================================================

#' Map standard amino acids to AnaCoDa convention
#' 
#' AnaCoDa uses "Z" for the Serine AGN codons (AGC, AGT) to separate them
#' from the TCN Serine codons (TCA, TCC, TCG, TCT)
#'
#' @param codon_df Data frame with AA and Codon columns
#' @return Data frame with AA column updated to AnaCoDa convention
#' @export
map_aa_to_anacoda <- function(codon_df) {
  # AGC and AGT are coded as "Z" in AnaCoDa
  codon_df$AA_anacoda <- codon_df$AA
  codon_df$AA_anacoda[codon_df$Codon %in% c("AGC", "AGT")] <- "Z"
  
  # Single-letter codes for other amino acids
  aa_map <- c(
    "A" = "A", "C" = "C", "D" = "D", "E" = "E", "F" = "F",
    "G" = "G", "H" = "H", "I" = "I", "K" = "K", "L" = "L",
    "M" = "M", "N" = "N", "P" = "P", "Q" = "Q", "R" = "R",
    "S" = "S", "T" = "T", "V" = "V", "W" = "W", "Y" = "Y"
  )
  
  # Replace full names with single letters if needed
  if (any(nchar(codon_df$AA_anacoda) > 1)) {
    codon_df$AA_anacoda <- aa_map[codon_df$AA_anacoda]
  }
  
  return(codon_df)
}

#' Validate ROC model predictions against observed frequencies
#'
#' @param observed_df Data frame with observed codon frequencies (from calculate_observed_codon_frequencies)
#' @param predictions_df Data frame with predicted probabilities (from predict_codon_probs_batch)
#' @param csp_df CSP parameters data frame
#' @return List with validation statistics and merged data
#' @export
validate_roc_model <- function(observed_df, predictions_df, csp_df) {
  
  # Map observed AA to AnaCoDa convention
  observed_df <- map_aa_to_anacoda(observed_df)
  
  # Merge predictions with observations
  validation_data <- merge(
    observed_df,
    predictions_df,
    by.x = c("Gene", "AA_anacoda", "Codon"),
    by.y = c("Gene", "AA", "Codon"),
    all = FALSE
  )
  
  message(sprintf("Validation dataset: %d observations from %d genes",
                  nrow(validation_data), length(unique(validation_data$Gene))))
  
  # Overall correlation
  cor_result <- cor.test(validation_data$predicted_prob, 
                         validation_data$Observed_freq,
                         method = "pearson")
  
  # Correlation by amino acid
  cor_by_aa <- do.call(rbind, lapply(split(validation_data, validation_data$AA_anacoda), function(df) {
    if (nrow(df) < 3) return(NULL)
    data.frame(
      AA = df$AA_anacoda[1],
      r = cor(df$predicted_prob, df$Observed_freq, use = "complete.obs"),
      n = nrow(df),
      stringsAsFactors = FALSE
    )
  }))
  cor_by_aa <- cor_by_aa[order(-cor_by_aa$r), ]
  
  # Mean absolute error
  validation_data$abs_error <- abs(validation_data$Observed_freq - validation_data$predicted_prob)
  mae <- mean(validation_data$abs_error, na.rm = TRUE)
  
  # Add selection direction based on dEta
  validation_data <- merge(validation_data, csp_df[, c("AA", "Codon", "dEta")],
                           by.x = c("AA_anacoda", "Codon"),
                           by.y = c("AA", "Codon"),
                           all.x = TRUE)
  
  validation_data$Selection_direction <- ifelse(
    validation_data$dEta < 0, "Preferred",
    ifelse(validation_data$dEta > 0, "Unpreferred", "Reference")
  )
  
  result <- list(
    data = validation_data,
    overall_correlation = cor_result,
    correlation_by_aa = cor_by_aa,
    mean_absolute_error = mae,
    n_genes = length(unique(validation_data$Gene)),
    n_observations = nrow(validation_data)
  )
  
  class(result) <- c("roc_validation", class(result))
  return(result)
}

#' Print summary of ROC validation results
#'
#' @param x roc_validation object
#' @param ... Additional arguments (ignored)
#' @export
print.roc_validation <- function(x, ...) {
  cat("\n=== ROC Model Validation Summary ===\n\n")
  cat(sprintf("Genes analyzed: %d\n", x$n_genes))
  cat(sprintf("Total observations: %d\n", x$n_observations))
  cat(sprintf("\nOverall Pearson r: %.4f (p = %.2e)\n", 
              x$overall_correlation$estimate, 
              x$overall_correlation$p.value))
  cat(sprintf("Mean absolute error: %.4f\n", x$mean_absolute_error))
  cat("\nCorrelation by amino acid (top 5):\n")
  print(head(x$correlation_by_aa, 5))
  invisible(x)
}

# ==============================================================================
# VISUALIZATION FUNCTIONS
# ==============================================================================

#' Create predicted vs observed scatter plot
#'
#' @param validation_result Result from validate_roc_model
#' @param output_file Path to save PDF (optional)
#' @return ggplot object
#' @export
plot_pred_vs_obs <- function(validation_result, output_file = NULL) {
  
  require(ggplot2)
  
  data <- validation_result$data
  r_val <- validation_result$overall_correlation$estimate
  n_genes <- validation_result$n_genes
  
  p <- ggplot(data, aes(x = predicted_prob, y = Observed_freq)) +
    geom_point(alpha = 0.05, size = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
    geom_smooth(method = "lm", color = "blue", se = TRUE, linewidth = 1) +
    labs(
      x = "Predicted Codon Frequency (ROC model)",
      y = "Observed Codon Frequency",
      title = "ROC Model Validation: Empirical Expression as Phi",
      subtitle = sprintf("r = %.3f, n = %d genes", r_val, n_genes)
    ) +
    theme_bw(base_size = 12) +
    coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1))
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 8, height = 8)
    message(sprintf("Saved: %s", output_file))
  }
  
  return(p)
}

#' Create faceted plot by amino acid
#'
#' @param validation_result Result from validate_roc_model
#' @param output_file Path to save PDF (optional)
#' @return ggplot object
#' @export
plot_by_amino_acid <- function(validation_result, output_file = NULL) {
  
  require(ggplot2)
  
  data <- validation_result$data
  
  p <- ggplot(data, aes(x = predicted_prob, y = Observed_freq)) +
    geom_point(alpha = 0.1, size = 0.3) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    facet_wrap(~ AA_anacoda, scales = "free") +
    labs(
      x = "Predicted Frequency",
      y = "Observed Frequency",
      title = "ROC Model Validation by Amino Acid"
    ) +
    theme_bw(base_size = 10) +
    theme(strip.text = element_text(size = 8))
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 14, height = 12)
    message(sprintf("Saved: %s", output_file))
  }
  
  return(p)
}

#' Plot codon usage vs expression by selection direction
#'
#' @param validation_result Result from validate_roc_model
#' @param n_quantiles Number of expression quantiles (default: 4)
#' @param output_file Path to save PDF (optional)
#' @return ggplot object
#' @export
plot_selection_vs_expression <- function(validation_result, n_quantiles = 4, 
                                          output_file = NULL) {
  
  require(ggplot2)
  
  data <- validation_result$data
  
  # Create expression quantiles
  data$Expression_Group <- cut(
    data$Phi_empirical,
    breaks = quantile(data$Phi_empirical, probs = seq(0, 1, length.out = n_quantiles + 1)),
    labels = paste0("Q", 1:n_quantiles),
    include.lowest = TRUE
  )
  
  # Summarize by expression group and selection direction
  usage_summary <- do.call(rbind, lapply(
    split(data[data$Selection_direction != "Reference", ], 
          list(data$Expression_Group[data$Selection_direction != "Reference"],
               data$Selection_direction[data$Selection_direction != "Reference"])),
    function(df) {
      if (nrow(df) == 0) return(NULL)
      data.frame(
        Expression_Group = df$Expression_Group[1],
        Selection_direction = df$Selection_direction[1],
        mean_freq = mean(df$Observed_freq, na.rm = TRUE),
        se_freq = sd(df$Observed_freq, na.rm = TRUE) / sqrt(nrow(df)),
        stringsAsFactors = FALSE
      )
    }
  ))
  
  p <- ggplot(usage_summary, 
              aes(x = Expression_Group, y = mean_freq, 
                  color = Selection_direction, group = Selection_direction)) +
    geom_point(size = 3) +
    geom_line(linewidth = 1) +
    geom_errorbar(aes(ymin = mean_freq - se_freq, ymax = mean_freq + se_freq), 
                  width = 0.2) +
    labs(
      x = "Expression Quantile",
      y = "Mean Observed Codon Frequency",
      color = "Selection Direction\n(based on dEta)",
      title = "Codon Usage vs Expression by Selection Direction",
      subtitle = "Preferred codons (dEta < 0) should increase with expression"
    ) +
    scale_color_manual(values = c("Preferred" = "#2166AC", "Unpreferred" = "#B2182B")) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 10, height = 7)
    message(sprintf("Saved: %s", output_file))
  }
  
  return(p)
}

#' Test biological prediction: preferred codons should increase with expression
#'
#' @param validation_result Result from validate_roc_model
#' @return List with test results
#' @export
test_selection_expression_relationship <- function(validation_result) {
  
  data <- validation_result$data
  
  # For preferred codons (dEta < 0), calculate mean usage per gene
  preferred_data <- data[data$dEta < 0, ]
  
  gene_preferred_usage <- do.call(rbind, lapply(
    split(preferred_data, preferred_data$Gene),
    function(df) {
      data.frame(
        Gene = df$Gene[1],
        Exp_log10 = df$Phi_empirical[1],
        mean_preferred_usage = mean(df$Observed_freq, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  ))
  
  # Correlation test
  cor_result <- cor.test(gene_preferred_usage$Exp_log10, 
                         gene_preferred_usage$mean_preferred_usage)
  
  # Interpretation
  if (cor_result$estimate > 0 && cor_result$p.value < 0.05) {
    interpretation <- "VALIDATED: Preferred codons (dEta < 0) are used MORE in highly expressed genes"
  } else if (cor_result$estimate < 0 && cor_result$p.value < 0.05) {
    interpretation <- "CONTRADICTED: Preferred codons (dEta < 0) are used LESS in highly expressed genes. Check dEta sign convention!"
  } else {
    interpretation <- "INCONCLUSIVE: No significant relationship between preferred codon usage and expression"
  }
  
  result <- list(
    correlation = cor_result,
    data = gene_preferred_usage,
    interpretation = interpretation
  )
  
  message("\n=== Selection-Expression Relationship Test ===")
  message(sprintf("Correlation: r = %.4f, p = %.2e", cor_result$estimate, cor_result$p.value))
  message(interpretation)
  
  return(result)
}
