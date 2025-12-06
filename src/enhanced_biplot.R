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
  #' Create a clean biplot showing preferred vs non-preferred codons

  #' 
  #' @param ordination_result CA or PCA result object from FactoMineR

  #' @param gene_data data.frame with Gene_name and expression_group columns

  #' @param preferred_codons data.frame with Preferred_Codons column (from ROC)
  #' @param dims numeric vector of length 2, which dimensions to plot
  #' @param arrow_scale numeric, scaling factor for codon arrows
  #' @param title character, plot title
  #' @param subtitle character, optional subtitle
  #' @param output_file character, path to save PDF (NULL = no save)
  #' @return ggplot object
  
  require(ggplot2)
  require(dplyr)
  require(ggnewscale)
  
  # Detect analysis type and extract coordinates
  if (inherits(ordination_result, "PCA")) {
    # PCA from FactoMineR
    gene_scores <- as.data.frame(ordination_result$ind$coord[, dims])
    codon_loadings <- as.data.frame(ordination_result$var$coord[, dims])
    var_explained <- ordination_result$eig[dims, 2]
    analysis_type <- "PCA"
  } else if (inherits(ordination_result, "CA")) {
    # CA from FactoMineR
    gene_scores <- as.data.frame(ordination_result$row$coord[, dims])
    codon_loadings <- as.data.frame(ordination_result$col$coord[, dims])
    var_explained <- ordination_result$eig[dims, 2]
    analysis_type <- "CA"
  } else if (inherits(ordination_result, "coa") || "li" %in% names(ordination_result)) {
    # ade4 style or converted CA object
    gene_scores <- as.data.frame(ordination_result$li[, dims])
    codon_loadings <- as.data.frame(ordination_result$co[, dims])
    if ("eig" %in% names(ordination_result) && is.matrix(ordination_result$eig)) {
      var_explained <- ordination_result$eig[dims, 2]
    } else if ("eig" %in% names(ordination_result)) {
      total <- sum(ordination_result$eig)
      var_explained <- ordination_result$eig[dims] / total * 100
    } else {
      var_explained <- c(NA, NA)
    }
    analysis_type <- "CA"
  } else {
    stop("Unsupported ordination result type")
  }
  
  # Standardize column names
  names(gene_scores) <- c("Dim1", "Dim2")
  names(codon_loadings) <- c("Dim1", "Dim2")
  
  # Add gene names
  gene_scores$Gene_name <- sub("\\.1$", "", rownames(gene_scores))
  codon_loadings$Codon <- rownames(codon_loadings)
  
  # Merge with gene data (expression groups)
  gene_scores <- gene_scores |>
    dplyr::left_join(gene_data, by = "Gene_name") |>
    dplyr::filter(!is.na(expression_group))
  
  # Classify codons as preferred vs non-preferred
  # Handle different column name formats
  if ("Preferred_Codons" %in% names(preferred_codons)) {
    preferred_list <- preferred_codons$Preferred_Codons
  } else if ("codon" %in% names(preferred_codons)) {
    preferred_list <- preferred_codons$codon
  } else if ("Codon" %in% names(preferred_codons)) {
    preferred_list <- preferred_codons$Codon
  } else {
    preferred_list <- preferred_codons[[1]]
  }
  
  codon_loadings <- codon_loadings |>
    dplyr::mutate(
      Preference = ifelse(Codon %in% preferred_list, "Preferred", "Non-preferred")
    )
  
  # Auto-scale arrows to be visible relative to gene score spread
  # Calculate the range of gene scores to determine appropriate arrow scaling
  gene_range <- max(
    diff(range(gene_scores$Dim1, na.rm = TRUE)),
    diff(range(gene_scores$Dim2, na.rm = TRUE))
  )
  codon_range <- max(
    diff(range(codon_loadings$Dim1, na.rm = TRUE)),
    diff(range(codon_loadings$Dim2, na.rm = TRUE))
  )
  
  # Scale so arrows span ~40% of the gene score range
  auto_scale <- (gene_range * 0.4) / codon_range
  effective_scale <- arrow_scale * auto_scale
  
  # Scale arrows
  codon_loadings <- codon_loadings |>
    dplyr::mutate(
      Dim1_scaled = Dim1 * effective_scale,
      Dim2_scaled = Dim2 * effective_scale
    )
  
  # Define colors - Preferred = red (matches high expression/selection)
  preference_colors <- c("Preferred" = "#E41A1C", "Non-preferred" = "#377EB8")
  
  # Build axis labels
  if (!is.na(var_explained[1])) {
    x_label <- sprintf("Dim %d (%.1f%%)", dims[1], var_explained[1])
    y_label <- sprintf("Dim %d (%.1f%%)", dims[2], var_explained[2])
  } else {
    x_label <- sprintf("Dim %d", dims[1])
    y_label <- sprintf("Dim %d", dims[2])
  }
  
  # Define ellipse colors (stronger for visibility)
  ellipse_colors <- c(
    "Top 5%" = "#E41A1C",
    "Bottom 5%" = "#377EB8",
    "High Selection (Top 5%)" = "#E41A1C",
    "Low Selection (Bottom 5%)" = "#377EB8"
  )
  
  # Create the plot
  p <- ggplot() +
    # Gene points with expression group colors (light, transparent)
    geom_point(
      data = gene_scores,
      aes(x = Dim1, y = Dim2, fill = expression_group),
      alpha = 0.25,
      size = 1.5,
      shape = 21,
      color = NA
    ) +
    scale_fill_manual(
      name = "Gene Group",
      values = c("Top 5%" = "#FCBBA1", "Bottom 5%" = "#9ECAE1",
                 "High Selection (Top 5%)" = "#FCBBA1", 
                 "Low Selection (Bottom 5%)" = "#9ECAE1")
    ) +
    # Add new fill scale for points
    ggnewscale::new_scale_fill() +
    # Add confidence ellipses for each group (strong colors)
    stat_ellipse(
      data = gene_scores,
      aes(x = Dim1, y = Dim2, color = expression_group),
      level = 0.95,
      linewidth = 1.2,
      linetype = "solid"
    ) +
    # Codon arrows colored by preference (thicker, more visible)
    geom_segment(
      data = codon_loadings,
      aes(x = 0, y = 0, xend = Dim1_scaled, yend = Dim2_scaled, 
          color = Preference),
      arrow = arrow(length = unit(0.2, "cm")),
      linewidth = 1.0,
      alpha = 1.0
    ) +
    # Codon labels (larger)
    geom_text(
      data = codon_loadings,
      aes(x = Dim1_scaled * 1.08, y = Dim2_scaled * 1.08, 
          label = Codon, color = Preference),
      size = 3.0,
      fontface = "bold",
      show.legend = FALSE
    ) +
    # Combined color scale for ellipses and arrows
    scale_color_manual(
      name = "",
      values = c(preference_colors, ellipse_colors),
      breaks = c("Preferred", "Non-preferred", 
                 names(ellipse_colors)[names(ellipse_colors) %in% unique(gene_scores$expression_group)])
    ) +
    # Axis labels
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label
    ) +
    # Reference lines
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
    # Apply custom theme
    theme_custom() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "right"
    )
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 10, height = 8, dpi = 300)
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
