plot_family_scurves <- function(model_result, meta_dt, 
                               expression_range = NULL,
                               output_file = NULL,
                               alpha_significance = 0.05,
                               preferred_codons_updated = NULL,
                               n_points = 1000) {
  #' Plots S-curves showing codon frequency vs expression for an amino acid family
  #' Handles multiple preferred codons and no-preference cases
  #'
  #' @param model_result Output from fit_pairwise_gams or fit_pairwise_glms
  #' @param meta_dt The integrated_data for calculating predictions
  #' @param expression_range Vector c(min, max) or NULL for auto
  #' @param output_file Path to save PDF, or NULL to return plot object
  #' @param alpha_significance Significance threshold (default 0.05)
  #' @param preferred_codons_updated Optional updated preferences from update_preferred_codons_from_models
  #' @param n_points Number of points to defined the prediction grid
  #' 
  #' @return ggplot object
  
  suppressPackageStartupMessages({
    require(ggplot2)
    require(dplyr)
  })
  
  family_name <- model_result$family
  baseline_codon <- model_result$baseline
  
  # Get predictions
  plot_data_long <- predict_codon_frequencies(
    model_result = model_result,
    meta_dt = meta_dt,
    expression_range = expression_range,
    n_points = n_points
  )
  
  # Get coefficient table with significance info
  coef_table <- model_result$coefficients
  
  # Determine significance for each codon
  # Use p_adj if available, otherwise p_value
  if ("p_adj" %in% names(coef_table)) {
    sig_info <- coef_table %>%
      dplyr::select(Codon, Selection_Slope, p_value, p_adj) %>%
      dplyr::mutate(
        Significant = p_adj < alpha_significance,
        Is_Baseline = (Codon == baseline_codon)
      )
  } else {
    sig_info <- coef_table %>%
      dplyr::select(Codon, Selection_Slope, p_value) %>%
      dplyr::mutate(
        p_adj = p_value,
        Significant = p_value < alpha_significance,
        Is_Baseline = (Codon == baseline_codon)
      )
  }
  
  # Identify preferred codon(s) from updated preferences or default to baseline
  if (!is.null(preferred_codons_updated)) {
    fam_prefs <- preferred_codons_updated %>%
      dplyr::filter(Family == family_name)
    
    if (nrow(fam_prefs) > 0) {
      pref_pattern <- unique(fam_prefs$Preference_Pattern)[1]
      preferred_codons <- fam_prefs$Preferred_Codons[1]
      
      # Parse multiple preferred codons
      if (grepl(";", preferred_codons)) {
        preferred_list <- strsplit(preferred_codons, ";")[[1]]
      } else if (preferred_codons == "None") {
        preferred_list <- character(0)  # No preference
      } else {
        preferred_list <- preferred_codons
      }
    } else {
      preferred_list <- baseline_codon
      pref_pattern <- "Single_Preferred_Baseline"
    }
  } else {
    # Default: use baseline
    preferred_list <- baseline_codon
    pref_pattern <- "Single_Preferred_Baseline"
  }
  
  # Merge significance info with plot data
  plot_data_long <- plot_data_long %>%
    dplyr::left_join(sig_info, by = "Codon")
  
  # Create visual attributes
  plot_data_long <- plot_data_long %>%
    dplyr::mutate(
      Is_Preferred = Codon %in% preferred_list,
      Line_Type = ifelse(Significant | Is_Baseline, "solid", "dashed"),
      Codon_Label = dplyr::case_when(
        length(preferred_list) == 0 ~ Codon,  # No preference
        Is_Preferred & length(preferred_list) > 1 ~ paste0(Codon, " ★"),  # Multiple
        Is_Preferred ~ paste0(Codon, " ★"),  # Single preferred
        TRUE ~ Codon
      )
    )
  
  # Create subtitle based on preference pattern
  if (length(preferred_list) == 0) {
    subtitle_text <- sprintf("No clear preference (all non-significant) | Baseline: %s", baseline_codon)
  } else if (length(preferred_list) > 1) {
    subtitle_text <- sprintf("Co-optimal codons ★: %s | Baseline: %s", 
                            paste(preferred_list, collapse = ", "), baseline_codon)
  } else {
    subtitle_text <- sprintf("Preferred: %s ★ | Baseline: %s | Controlling for gene length and GC12",
                            preferred_list[1], baseline_codon)
  }
  
  # Create plot with better visual encoding
  p <- ggplot(plot_data_long, 
              aes(x = High_exp_log2, y = Predicted_Frequency, 
                  color = Codon_Label, linetype = Line_Type)) +
    geom_line(linewidth = 1.2, alpha = 0.85) +
    scale_linetype_manual(
      values = c("solid" = "solid", "dashed" = "dashed"),
      labels = c("solid" = "Significant selection (FDR < 0.05) or baseline",
                "dashed" = "Non-significant (no clear selection)"),
      name = "Selection Evidence"
    ) +
    labs(
      title = paste("Codon Frequency vs Expression:", family_name),
      subtitle = subtitle_text,
      x = "Gene Expression (log2 CPM)",
      y = "Predicted Codon Frequency",
      color = "Codon"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      legend.box = "vertical",
      legend.spacing.y = unit(0.3, "cm")
    ) +
    guides(
      color = guide_legend(order = 1, title = "Codon"),
      linetype = guide_legend(order = 2, title = "Selection Evidence")
    )
  
  # Save if output file provided
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 10, height = 6)
    cat(sprintf("✓ Plot saved: %s\n", output_file))
  }
  
  return(p)
}


plot_all_families_panel <- function(all_model_results, meta_dt, 
                                   output_file = "./results/section_8.1_all_families_scurves.pdf",
                                   ncol = 4,
                                   preferred_codons_updated = NULL) {
  #' Creates a multi-panel plot with S-curves for all amino acid families
  #' Handles multiple preferred codons and no-preference cases
  #'
  #' @param all_model_results List of model results from lapply
  #' @param meta_dt The integrated_data
  #' @param output_file Path to save PDF
  #' @param ncol Number of columns in facet
  #' @param preferred_codons_updated Optional updated preferences
  #' @return ggplot object
  
  suppressPackageStartupMessages({
    require(ggplot2)
    require(dplyr)
    require(tidyr)
  })
  
  # Collect all predictions
  all_predictions <- lapply(names(all_model_results), function(fam_name) {
    result <- all_model_results[[fam_name]]
    if (is.null(result)) return(NULL)
    
    pred <- predict_codon_frequencies(
      model_result = result,
      meta_dt = meta_dt,
      n_points = n_points  # Fewer points for speed
    )
    pred$Family <- fam_name
    return(pred)
  })
  
  # Combine all predictions
  all_pred_df <- data.table::rbindlist(all_predictions)
  
  # Identify preferred codons and significance per family
  all_coef_tables <- lapply(all_model_results, function(x) {
    if (is.null(x)) return(NULL)
    x$coefficients %>% 
      dplyr::mutate(Family = x$family, Baseline = x$baseline)
  })
  
  all_coefs <- data.table::rbindlist(all_coef_tables)
  
  # Get preferred codons from updated preferences or default to baseline
  if (!is.null(preferred_codons_updated)) {
    # Use updated preferences
    preferred_codons <- preferred_codons_updated %>%
      dplyr::select(Family, Preferred_Codons, Preference_Pattern) %>%
      dplyr::distinct()
    
    # Expand families with multiple codons for joining
    preferred_expanded <- lapply(1:nrow(preferred_codons), function(i) {
      row <- preferred_codons[i, ]
      if (grepl(";", row$Preferred_Codons)) {
        codons <- strsplit(row$Preferred_Codons, ";")[[1]]
        data.frame(
          Family = rep(row$Family, length(codons)),
          Preferred_Codon = codons,
          Preference_Pattern = rep(row$Preference_Pattern, length(codons)),
          stringsAsFactors = FALSE
        )
      } else if (row$Preferred_Codons == "None") {
        data.frame(
          Family = row$Family,
          Preferred_Codon = NA,
          Preference_Pattern = row$Preference_Pattern,
          stringsAsFactors = FALSE
        )
      } else {
        data.frame(
          Family = row$Family,
          Preferred_Codon = row$Preferred_Codons,
          Preference_Pattern = row$Preference_Pattern,
          stringsAsFactors = FALSE
        )
      }
    })
    preferred_codons <- data.table::rbindlist(preferred_expanded)
  } else {
    # Default: use baseline as preferred
    preferred_codons <- all_coefs %>%
      dplyr::group_by(Family) %>%
      dplyr::slice(1) %>%
      dplyr::select(Family, Preferred_Codon = Baseline) %>%
      dplyr::distinct() %>%
      dplyr::mutate(Preference_Pattern = "Single_Preferred_Baseline")
  }
  
  # Get significance info
  sig_info <- all_coefs %>%
    dplyr::select(Family, Codon, Baseline) %>%
    dplyr::left_join(
      all_coefs %>% 
        dplyr::select(Family, Codon, Significant = Significant, p_adj),
      by = c("Family", "Codon")
    ) %>%
    dplyr::mutate(
      Significant = ifelse(is.na(Significant), FALSE, Significant),
      Is_Baseline = (Codon == Baseline)
    )
  
  # Merge with predictions
  all_pred_df <- all_pred_df %>%
    dplyr::left_join(preferred_codons, by = "Family") %>%
    dplyr::left_join(sig_info, by = c("Family", "Codon")) %>%
    dplyr::mutate(
      Is_Preferred = (Codon == Preferred_Codon),
      Codon_Type = dplyr::case_when(
        Is_Preferred ~ "Preferred",
        Significant | Is_Baseline ~ "Significant",
        TRUE ~ "Other"
      ),
      Line_Width = ifelse(Is_Preferred, 1.2, 0.7)
    )
  
  # Create faceted plot with enhanced visual encoding
  p <- ggplot(all_pred_df, 
              aes(x = High_exp_log2, y = Predicted_Frequency, 
                  color = Codon, group = Codon,
                  linetype = Codon_Type,
                  linewidth = Line_Width,
                  alpha = Codon_Type)) +
    geom_line() +
    scale_linetype_manual(
      values = c("Preferred" = "solid", "Significant" = "solid", "Other" = "dotted"),
      name = "",
      labels = c("Preferred" = "Preferred (highest slope)",
                "Significant" = "Significant selection",
                "Other" = "Non-significant")
    ) +
    scale_linewidth_identity() +
    scale_alpha_manual(
      values = c("Preferred" = 1.0, "Significant" = 0.7, "Other" = 0.4),
      guide = "none"
    ) +
    facet_wrap(~ Family, scales = "free_y", ncol = ncol) +
    labs(
      title = "Codon Frequencies vs Expression: All Amino Acid Families",
      subtitle = "Heavy solid = preferred codon | Light solid = significant | Dotted = not significant (FDR > 0.05)",
      x = "Gene Expression (log2 CPM)",
      y = "Predicted Frequency",
      color = "Codon"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold", size = 9),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "gray30")
    ) +
    guides(
      color = "none",
      linetype = guide_legend(order = 1, override.aes = list(linewidth = 1))
    )
  
  # Save
  ggsave(output_file, p, width = 14, height = 12)
  cat(sprintf("✓ Multi-panel plot saved: %s\n", output_file))
  
  return(p)
}


compare_gam_vs_glm <- function(gam_results, glm_results, 
                              output_file = "./results/section_8.1_GAM_vs_GLM_comparison.pdf") {
  #' Compares selection slopes from GAM vs GLM models
  #'
  #' @param gam_results List of GAM model results
  #' @param glm_results List of GLM (Box-Cox) model results
  #' @param output_file Path to save comparison plot
  #' @return Data frame with comparison statistics
  
  suppressPackageStartupMessages({
    require(ggplot2)
    require(dplyr)
  })
  
  # Extract coefficients from GAM models
  gam_coefs <- lapply(gam_results, function(x) {
    if (is.null(x)) return(NULL)
    x$coefficients %>% dplyr::mutate(Model = "GAM")
  })
  gam_df <- data.table::rbindlist(gam_coefs)
  
  # Extract coefficients from GLM models
  glm_coefs <- lapply(glm_results, function(x) {
    if (is.null(x)) return(NULL)
    x$coefficients %>% dplyr::mutate(Model = "GLM (Box-Cox)")
  })
  glm_df <- data.table::rbindlist(glm_coefs)
  
  # Merge by codon and family
  comparison_df <- gam_df %>%
    dplyr::select(Codon, Family, Baseline, GAM_Slope = Selection_Slope, 
                  GAM_SE = SE, GAM_p = p_value) %>%
    dplyr::inner_join(
      glm_df %>% dplyr::select(Codon, Family, GLM_Slope = Selection_Slope,
                               GLM_SE = SE, GLM_p = p_value),
      by = c("Codon", "Family")
    ) %>%
    dplyr::filter(Codon != Baseline)  # Remove baseline codons (slope = 0)
  
  # Calculate correlation
  cor_slopes <- cor(comparison_df$GAM_Slope, comparison_df$GLM_Slope, 
                   use = "complete.obs")
  
  # Create comparison plot
  p1 <- ggplot(comparison_df, aes(x = GAM_Slope, y = GLM_Slope)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "lm", se = TRUE, color = "blue", alpha = 0.2) +
    labs(
      title = "Selection Slopes: GAM vs GLM (Box-Cox)",
      subtitle = sprintf("Correlation: %.3f | Red line = perfect agreement | Blue = linear fit", 
                        cor_slopes),
      x = "Selection Slope (GAM with smoothers)",
      y = "Selection Slope (GLM with Box-Cox confounders)"
    ) +
    theme_bw()
  
  # P-value comparison
  comparison_df <- comparison_df %>%
    dplyr::mutate(
      GAM_Significant = GAM_p < 0.05,
      GLM_Significant = GLM_p < 0.05,
      Agreement = case_when(
        GAM_Significant & GLM_Significant ~ "Both Significant",
        !GAM_Significant & !GLM_Significant ~ "Both Non-Significant",
        GAM_Significant & !GLM_Significant ~ "GAM Only",
        !GAM_Significant & GLM_Significant ~ "GLM Only"
      )
    )
  
  # Agreement table
  agreement_table <- table(comparison_df$Agreement)
  
  p2 <- ggplot(comparison_df, aes(x = -log10(GAM_p), y = -log10(GLM_p), 
                                 color = Agreement)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Both Significant" = "#E41A1C",
                                  "Both Non-Significant" = "#999999",
                                  "GAM Only" = "#377EB8",
                                  "GLM Only" = "#4DAF4A")) +
    labs(
      title = "Significance Comparison: GAM vs GLM",
      subtitle = "Dashed lines = p = 0.05 threshold",
      x = "-log10(p-value) GAM",
      y = "-log10(p-value) GLM (Box-Cox)",
      color = "Agreement"
    ) +
    theme_bw() +
    theme(legend.position = "right")
  
  # Combine plots
  p_combined <- gridExtra::grid.arrange(p1, p2, ncol = 2)
  
  # Save
  ggsave(output_file, p_combined, width = 14, height = 6)
  cat(sprintf("✓ Comparison plot saved: %s\n", output_file))
  
  # Print summary statistics
  cat("\n=== GAM vs GLM Comparison Summary ===\n")
  cat(sprintf("Correlation of selection slopes: %.3f\n", cor_slopes))
  cat(sprintf("Mean absolute difference: %.4f\n", 
              mean(abs(comparison_df$GAM_Slope - comparison_df$GLM_Slope), na.rm = TRUE)))
  cat("\nSignificance agreement:\n")
  print(agreement_table)
  cat("\n")
  
  return(comparison_df)
}
