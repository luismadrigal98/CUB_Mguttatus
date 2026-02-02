# Enhanced biplot functions for CA/PCA analysis
# Simplified version: shows preferred vs non-preferred codons
# =============================================================================

create_preference_biplot <- function(
    ordination_result,
    gene_data,
    preferred_codons,
    dims = c(1, 2),
    arrow_scale = 1.0,
    title = "Codon Usage Biplot",
    subtitle = NULL,
    output_file = NULL
) {
  require(ggplot2)
  require(dplyr)
  
  # 1. Extract Coordinates
  if (inherits(ordination_result, "CA")) {
    gene_scores <- as.data.frame(ordination_result$row$coord[, dims])
    codon_loadings <- as.data.frame(ordination_result$col$coord[, dims])
  } else {
    gene_scores <- as.data.frame(ordination_result$ind$coord[, dims])
    codon_loadings <- as.data.frame(ordination_result$var$coord[, dims])
  }
  
  # Standardize internal names
  names(gene_scores) <- c("DimX", "DimY")
  names(codon_loadings) <- c("DimX", "DimY")
  
  # Clean gene names to ensure matching
  gene_scores$Gene_name <- sub("\\.1$", "", rownames(gene_scores))
  codon_loadings$Codon <- rownames(codon_loadings)
  
  # 2. Merge and Label the Middle 90%
  gene_scores <- gene_scores %>%
    dplyr::left_join(gene_data, by = "Gene_name") %>%
    dplyr::mutate(expression_group = ifelse(is.na(expression_group), 
                                            "Middle 90%", 
                                            as.character(expression_group)))
  
  # 3. Codon Classification
  preferred_list <- if ("Preferred_Codons" %in% names(preferred_codons)) {
    preferred_codons$Preferred_Codons
  } else {
    preferred_codons[[1]]
  }
  
  codon_loadings <- codon_loadings %>%
    dplyr::mutate(Preference = ifelse(Codon %in% preferred_list, "Preferred", "Non-preferred"))
  
  # 4. Scaling Logic for Arrows
  gene_range <- max(diff(range(gene_scores$DimX)), diff(range(gene_scores$DimY)))
  codon_range <- max(diff(range(codon_loadings$DimX)), diff(range(codon_loadings$DimY)))
  effective_scale <- arrow_scale * ((gene_range * 0.4) / codon_range)
  
  codon_loadings <- codon_loadings %>%
    dplyr::mutate(
      DimX_scaled = DimX * effective_scale, 
      DimY_scaled = DimY * effective_scale
    )
  
  # 5. Define Unified Color Palette
  # This map handles Points, Ellipses, and Arrows in one scale
  unified_colors <- c(
    "Top 5%" = "#E41A1C", 
    "Bottom 5%" = "#377EB8", 
    "Middle 90%" = "gray85",
    "Preferred" = "#E41A1C",
    "Non-preferred" = "#377EB8"
  )
  
  # 6. Build the Plot
  p <- ggplot() +
    # Reference Lines
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.2) +
    
    # Gene Points: Using color instead of fill to avoid scale conflicts
    # Alpha is low (0.05) because 25k points will overlap heavily
    geom_point(
      data = gene_scores, 
      aes(x = DimX, y = DimY, color = expression_group),
      alpha = 0.05, 
      size = 0.8
    ) +
    
    # Confidence Ellipses
    stat_ellipse(
      data = gene_scores, 
      aes(x = DimX, y = DimY, color = expression_group),
      level = 0.95, 
      linewidth = 0.8
    ) +
    
    # Codon Arrows
    geom_segment(
      data = codon_loadings, 
      aes(x = 0, y = 0, xend = DimX_scaled, yend = DimY_scaled, color = Preference),
      arrow = arrow(length = unit(0.2, "cm")), 
      linewidth = 0.8
    ) +
    
    # Codon Labels
    geom_text(
      data = codon_loadings, 
      aes(x = DimX_scaled * 1.1, y = DimY_scaled * 1.1, label = Codon, color = Preference), 
      size = 2.8, 
      fontface = "bold",
      show.legend = FALSE
    ) +
    
    # Unified Color Scale
    scale_color_manual(name = "Group / Preference", values = unified_colors) +
    
    # Formatting
    theme_custom() +
    labs(
      title = title,
      subtitle = subtitle,
      x = paste("Dimension", dims[1]),
      y = paste("Dimension", dims[2])
    ) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  # 7. Save and Return
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 10, height = 8, bg = "white")
    cat(sprintf("✓ Saved: %s\n", output_file))
  }
  
  return(p)
}


analyze_codon_loading_direction <- function(
    ordination_result,
    preferred_codons,
    dim = 1
) {
  #' Analyze whether preferred codons load in a consistent direction
  #' 
  #' @param ordination_result CA or PCA result object
  #' @param preferred_codons data.frame with preferred codon information
  #' @param dim which dimension to analyze (default = 1)
  #' @return list with test statistics and summary
  
  # Extract codon loadings
  if (inherits(ordination_result, "PCA")) {
    loadings <- ordination_result$var$coord[, dim]
  } else if (inherits(ordination_result, "CA")) {
    loadings <- ordination_result$col$coord[, dim]
  } else if ("co" %in% names(ordination_result)) {
    loadings <- ordination_result$co[, dim]
  } else {
    stop("Unsupported ordination result type")
  }
  
  codon_df <- data.frame(
    Codon = names(loadings),
    Loading = as.numeric(loadings)
  )
  
  # Get preferred codon list
  if ("Preferred_Codons" %in% names(preferred_codons)) {
    preferred_list <- preferred_codons$Preferred_Codons
  } else if ("codon" %in% names(preferred_codons)) {
    preferred_list <- preferred_codons$codon
  } else if ("Codon" %in% names(preferred_codons)) {
    preferred_list <- preferred_codons$Codon
  } else {
    preferred_list <- preferred_codons[[1]]
  }
  
  codon_df <- codon_df |>
    dplyr::mutate(Preference = ifelse(Codon %in% preferred_list, "Preferred", "Non-preferred"))
  
  # Statistical tests
  # 1. Wilcoxon test: do preferred and non-preferred codons differ in loading?
  wtest <- wilcox.test(Loading ~ Preference, data = codon_df)
  
  # 2. Sign test: do preferred codons consistently load in one direction?
  preferred_loadings <- codon_df$Loading[codon_df$Preference == "Preferred"]
  n_positive <- sum(preferred_loadings > 0)
  n_negative <- sum(preferred_loadings < 0)
  sign_test <- binom.test(n_positive, n_positive + n_negative)
  
  # Summary statistics
  summary_stats <- codon_df |>
    dplyr::group_by(Preference) |>
    dplyr::summarise(
      Mean = mean(Loading),
      Median = median(Loading),
      SD = sd(Loading),
      N_positive = sum(Loading > 0),
      N_negative = sum(Loading < 0),
      .groups = "drop"
    )
  
  result <- list(
    dimension = dim,
    wilcoxon_test = wtest,
    sign_test = sign_test,
    summary = summary_stats,
    codon_data = codon_df
  )
  
  # Print summary
  cat(sprintf("\n=== Codon Loading Direction Analysis (Dim %d) ===\n", dim))
  cat("\nSummary by preference group:\n")
  print(as.data.frame(summary_stats))
  cat(sprintf("\nWilcoxon test (preferred vs non-preferred): W = %.2f, p = %.4f\n",
              wtest$statistic, wtest$p.value))
  cat(sprintf("Sign test (preferred codons): %d positive, %d negative, p = %.4f\n",
              n_positive, n_negative, sign_test$p.value))
  
  if (sign_test$p.value < 0.05) {
    direction <- ifelse(n_positive > n_negative, "positive", "negative")
    cat(sprintf("→ Preferred codons significantly load in the %s direction\n", direction))
  }
  
  return(result)
}
