##' Create enhanced CA/PCA biplot with codon classification coloring
##' 
##' @description Creates biplots showing gene clouds and codon vectors
##' with coloring based on codon classification (selection, preference, AT/GC)
##' 
##' @param ordination_result CA or PCA result object
##' @param gene_data Data frame with gene scores and expression groups
##' @param codon_test_results Output from test_codon_proportions()
##' @param w_table CAI weight table with preferred codons
##' @param dims Dimensions to plot (default: c(1,2))
##' @param color_by How to color codons: "selection", "preference", "ending", "combined"
##' @param show_only_significant Show only significant codons? (default: FALSE)
##' @param arrow_scale Scaling factor for codon vectors
##' @param title Plot title
##' @param output_file Path to save plot
##' 
##' @return ggplot object
##' 
##' @author Luis J. Madrigal-Roca
##' @date November 12, 2025

create_enhanced_biplot <- function(ordination_result,
                                   gene_data,
                                   codon_test_results,
                                   w_table,
                                   dims = c(1, 2),
                                   color_by = "combined",
                                   show_only_significant = FALSE,
                                   arrow_scale = 1.0,
                                   title = "Enhanced Biplot",
                                   output_file = NULL) {
  
  require(ggplot2)
  require(dplyr)
  require(data.table)
  
  cat(sprintf("\n=== Creating Enhanced Biplot ===\n"))
  cat(sprintf("Ordination type: %s\n", class(ordination_result)[1]))
  cat(sprintf("Dimensions: %d vs %d\n", dims[1], dims[2]))
  cat(sprintf("Color scheme: %s\n", color_by))
  
  # Extract ordination type
  is_ca <- inherits(ordination_result, "coa")
  is_pca <- inherits(ordination_result, c("prcomp", "princomp"))
  
  # Get gene scores
  if (is_ca) {
    gene_scores <- as.data.frame(ordination_result$li)
    codon_loadings <- as.data.frame(ordination_result$co)
    variance_explained <- ordination_result$eig / sum(ordination_result$eig) * 100
  } else if (is_pca) {
    gene_scores <- as.data.frame(ordination_result$x)
    codon_loadings <- as.data.frame(ordination_result$rotation)
    variance_explained <- (ordination_result$sdev^2) / sum(ordination_result$sdev^2) * 100
  } else {
    stop("Ordination result must be from CA (ade4::dudi.coa) or PCA (prcomp)")
  }
  
  # Add gene names
  gene_scores$Gene_name <- rownames(gene_scores)
  codon_loadings$Codon <- rownames(codon_loadings)
  
  # Merge with gene expression groups
  gene_plot_data <- gene_scores %>%
    left_join(gene_data %>% select(Gene_name, expression_group), 
              by = "Gene_name")
  
  # Merge codon loadings with test results and w values
  codon_plot_data <- codon_loadings %>%
    left_join(codon_test_results %>% 
                select(Codon, Classification, Significant, 
                       Difference, Ending, Amino_Acid),
              by = "Codon") %>%
    left_join(w_table %>% select(codon, relative_adaptiveness),
              by = c("Codon" = "codon")) %>%
    mutate(
      Preferred = relative_adaptiveness == 1.0,
      # Use corrected classification if available
      Category = if ("Combined_Classification" %in% names(codon_test_results)) {
        Classification
      } else {
        # Fallback to simple classification
        case_when(
          !Significant ~ "Neutral",
          Difference > 0 & Preferred ~ "Selection + Preferred (w=1)",
          Difference > 0 & !Preferred ~ "Under Selection (not pref)",
          Difference < 0 & Preferred ~ "Rel. Preferred (AA avoided)",
          Difference < 0 ~ "Avoided in High Expr",
          TRUE ~ "Neutral"
        )
      }
    )
  
  # Filter to significant only if requested
  if (show_only_significant) {
    codon_plot_data <- codon_plot_data %>% filter(Significant)
    cat(sprintf("Showing only significant codons: %d\n", nrow(codon_plot_data)))
  }
  
  # Select dimensions
  dim_names <- colnames(gene_scores)[dims]
  x_var <- variance_explained[dims[1]]
  y_var <- variance_explained[dims[2]]
  
  # Create column name mapping
  gene_plot_data$x <- gene_plot_data[[dim_names[1]]]
  gene_plot_data$y <- gene_plot_data[[dim_names[2]]]
  codon_plot_data$x <- codon_plot_data[[dim_names[1]]] * arrow_scale
  codon_plot_data$y <- codon_plot_data[[dim_names[2]]] * arrow_scale
  
  # Define color schemes
  if (color_by == "selection") {
    # Color by selection status only
    colors <- c(
      "Selection + Preferred (w=1)" = "#d73027",
      "Under Selection (not pref)" = "#fc8d59",
      "Avoided in High Expr" = "#4575b4",
      "Neutral" = "gray70"
    )
    codon_plot_data$Color_Variable <- codon_plot_data$Category
    legend_title <- "Selection Status"
    
  } else if (color_by == "preference") {
    # Color by CAI preference only
    colors <- c("TRUE" = "#d73027", "FALSE" = "gray70")
    codon_plot_data$Color_Variable <- codon_plot_data$Preferred
    legend_title <- "Preferred (w=1)"
    
  } else if (color_by == "ending") {
    # Color by AT vs GC ending
    colors <- c("AT" = "#fdae61", "GC" = "#abd9e9")
    codon_plot_data$Color_Variable <- codon_plot_data$Ending
    legend_title <- "Codon Ending"
    
  } else {  # combined
    # Combined classification
    codon_plot_data <- codon_plot_data %>%
      mutate(
        Combined = case_when(
          !Significant ~ "Neutral",
          Preferred & Difference > 0 & Ending == "GC" ~ "Sel + Pref + GC",
          Preferred & Difference > 0 & Ending == "AT" ~ "Sel + Pref + AT",
          !Preferred & Difference > 0 & Ending == "GC" ~ "Sel (non-pref) + GC",
          !Preferred & Difference > 0 & Ending == "AT" ~ "Sel (non-pref) + AT",
          Difference < 0 ~ "Avoided",
          TRUE ~ "Neutral"
        )
      )
    colors <- c(
      "Sel + Pref + GC" = "#d73027",
      "Sel + Pref + AT" = "#fc8d59",
      "Sel (non-pref) + GC" = "#fee090",
      "Sel (non-pref) + AT" = "#ffffbf",
      "Avoided" = "#4575b4",
      "Neutral" = "gray80"
    )
    codon_plot_data$Color_Variable <- codon_plot_data$Combined
    legend_title <- "Codon Classification"
  }
  
  # Gene cloud colors
  gene_colors <- c(
    "Top 5%" = "#e41a1c",
    "Bottom 5%" = "#377eb8",
    "Middle 90%" = "gray60"
  )
  
  # Create plot
  p <- ggplot() +
    # Gene clouds
    geom_point(data = gene_plot_data,
               aes(x = x, y = y, color = expression_group),
               alpha = 0.3, size = 1.5) +
    scale_color_manual(values = gene_colors, name = "Expression Group") +
    # Codon vectors
    geom_segment(data = codon_plot_data,
                 aes(x = 0, y = 0, xend = x, yend = y,
                     color = Color_Variable),
                 arrow = arrow(length = unit(0.02, "npc"), type = "closed"),
                 alpha = 0.7, size = 0.5,
                 inherit.aes = FALSE) +
    geom_text(data = codon_plot_data,
              aes(x = x * 1.1, y = y * 1.1, label = Codon),
              size = 2.5, alpha = 0.8) +
    # Styling
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
    labs(
      title = title,
      subtitle = sprintf("%d genes, %d codons | %s coloring",
                        nrow(gene_plot_data),
                        nrow(codon_plot_data),
                        color_by),
      x = sprintf("%s (%.1f%%)", dim_names[1], x_var),
      y = sprintf("%s (%.1f%%)", dim_names[2], y_var)
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "right",
      panel.grid.minor = element_blank()
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 3, alpha = 1))
    )
  
  # Add second color scale for codons (using ggnewscale if available)
  if (requireNamespace("ggnewscale", quietly = TRUE)) {
    p <- p + 
      ggnewscale::new_scale_color() +
      geom_segment(data = codon_plot_data,
                   aes(x = 0, y = 0, xend = x, yend = y,
                       color = Color_Variable),
                   arrow = arrow(length = unit(0.02, "npc"), type = "closed"),
                   alpha = 0.7, size = 0.5) +
      scale_color_manual(values = colors, name = legend_title)
  } else {
    cat("\nNote: Install 'ggnewscale' package for better legend handling\n")
  }
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 12, height = 10)
    cat(sprintf("\n✓ Biplot saved: %s\n", output_file))
  }
  
  return(p)
}


##' Create panel of biplots with different codon coloring schemes
##'
##' @param ordination_result CA or PCA result
##' @param gene_data Gene data with expression groups
##' @param codon_test_results Codon test results
##' @param w_table CAI weights
##' @param dims Dimensions to plot
##' @param output_file Path to save
##'
##' @return Combined plot object

create_biplot_panel <- function(ordination_result, gene_data, 
                                codon_test_results, w_table,
                                dims = c(1, 2), output_file = NULL) {
  
  require(cowplot)
  
  cat("\n=== Creating Biplot Panel ===\n")
  
  # Create plots with different color schemes
  p1 <- create_enhanced_biplot(
    ordination_result, gene_data, codon_test_results, w_table,
    dims = dims, color_by = "selection", show_only_significant = FALSE,
    title = "A) Selection Status"
  )
  
  p2 <- create_enhanced_biplot(
    ordination_result, gene_data, codon_test_results, w_table,
    dims = dims, color_by = "preference", show_only_significant = FALSE,
    title = "B) CAI Preference (w=1)"
  )
  
  p3 <- create_enhanced_biplot(
    ordination_result, gene_data, codon_test_results, w_table,
    dims = dims, color_by = "ending", show_only_significant = FALSE,
    title = "C) AT vs GC Ending"
  )
  
  p4 <- create_enhanced_biplot(
    ordination_result, gene_data, codon_test_results, w_table,
    dims = dims, color_by = "selection", show_only_significant = TRUE,
    title = "D) Significant Codons Only"
  )
  
  # Combine
  combined <- plot_grid(p1, p2, p3, p4, ncol = 2, nrow = 2)
  
  if (!is.null(output_file)) {
    ggsave(output_file, combined, width = 20, height = 18)
    cat(sprintf("\n✓ Panel saved: %s\n", output_file))
  }
  
  return(combined)
}


##' Analyze codon loading patterns in relation to gene groups
##'
##' @param ordination_result CA or PCA result
##' @param codon_test_results Codon test results
##' @param dims Dimensions to analyze
##'
##' @return Data frame with loading analysis

analyze_codon_loading_patterns <- function(ordination_result, 
                                           codon_test_results,
                                           dims = c(1, 2)) {
  
  cat("\n=== Analyzing Codon Loading Patterns ===\n")
  
  # Extract loadings
  if (inherits(ordination_result, "coa")) {
    loadings <- as.data.frame(ordination_result$co)
  } else {
    loadings <- as.data.frame(ordination_result$rotation)
  }
  
  loadings$Codon <- rownames(loadings)
  
  # Merge with test results
  analysis <- loadings %>%
    left_join(codon_test_results, by = "Codon") %>%
    mutate(
      Loading_Dim1 = .data[[colnames(loadings)[dims[1]]]],
      Loading_Dim2 = .data[[colnames(loadings)[dims[2]]]],
      Loading_Magnitude = sqrt(Loading_Dim1^2 + Loading_Dim2^2)
    )
  
  # Statistical tests
  cat("\n--- Testing Loading Patterns ---\n")
  
  # Test 1: Do significant codons have higher loadings?
  sig_loadings <- analysis$Loading_Magnitude[analysis$Significant]
  nonsig_loadings <- analysis$Loading_Magnitude[!analysis$Significant]
  wilcox_test <- wilcox.test(sig_loadings, nonsig_loadings)
  
  cat(sprintf("Significant vs non-significant loading magnitude:\n"))
  cat(sprintf("  Median (sig): %.3f\n", median(sig_loadings, na.rm = TRUE)))
  cat(sprintf("  Median (non-sig): %.3f\n", median(nonsig_loadings, na.rm = TRUE)))
  cat(sprintf("  Wilcoxon p-value: %.4f\n\n", wilcox_test$p.value))
  
  # Test 2: Do preferred codons load in positive direction?
  preferred_codons <- analysis %>% 
    filter(relative_adaptiveness == 1.0, Significant)
  
  if (nrow(preferred_codons) > 0) {
    cat(sprintf("Preferred + significant codons (n=%d):\n", nrow(preferred_codons)))
    cat(sprintf("  Mean Dim1 loading: %.3f\n", 
                mean(preferred_codons$Loading_Dim1, na.rm = TRUE)))
    cat(sprintf("  Mean Dim2 loading: %.3f\n", 
                mean(preferred_codons$Loading_Dim2, na.rm = TRUE)))
    
    # Test if mean loading on Dim1 is > 0
    t_test <- t.test(preferred_codons$Loading_Dim1, mu = 0)
    cat(sprintf("  t-test (Dim1 > 0): p = %.4f\n\n", t_test$p.value))
  }
  
  # Test 3: AT vs GC ending patterns
  at_loadings <- analysis %>% filter(Ending == "AT", Significant)
  gc_loadings <- analysis %>% filter(Ending == "GC", Significant)
  
  if (nrow(at_loadings) > 0 & nrow(gc_loadings) > 0) {
    cat("AT vs GC ending (significant codons):\n")
    cat(sprintf("  AT mean Dim1: %.3f (n=%d)\n", 
                mean(at_loadings$Loading_Dim1), nrow(at_loadings)))
    cat(sprintf("  GC mean Dim1: %.3f (n=%d)\n", 
                mean(gc_loadings$Loading_Dim1), nrow(gc_loadings)))
    
    wilcox_ending <- wilcox.test(at_loadings$Loading_Dim1, 
                                  gc_loadings$Loading_Dim1)
    cat(sprintf("  Wilcoxon p-value: %.4f\n", wilcox_ending$p.value))
  }
  
  return(analysis)
}
