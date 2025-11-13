create_simple_3d_gif <- function(pca_result,
                                 gene_data,
                                 codon_test_results = NULL,
                                 preferred_codons = NULL,
                                 dims = c(1, 2, 3),
                                 color_by = "expression",
                                 show_loadings = TRUE,
                                 loading_scale = 1.0,
                                 title = "3D PCA",
                                 output_file = "pca_3d.gif",
                                 n_frames = 60,
                                 width = 800,
                                 height = 600,
                                 point_size = 2,
                                 resolution = 150) {
  #' Create simple 3D PCA rotating GIF using static plots
  #' 
  #' Uses ggplot2 to create 2D projections at different angles
  #' No dependencies on rgl or other heavy 3D libraries
  #' 
  #' @param pca_result PCA result from FactoMineR::PCA()
  #' @param gene_data Data frame with gene metadata
  #' @param codon_test_results Optional: results from test_codon_proportions()
  #' @param preferred_codons Optional: preferred codons data frame
  #' @param dims Vector of 3 dimensions to plot
  #' @param color_by How to color points
  #' @param show_loadings Show codon loading vectors?
  #' @param loading_scale Scaling factor for loading arrows
  #' @param title Plot title
  #' @param output_file Output GIF file path
  #' @param n_frames Number of frames in rotation
  #' @param width Plot width in pixels
  #' @param height Plot height in pixels
  #' @param point_size Point size for genes
  #' @param resolution DPI for output
  #' 
  #' @return NULL (saves GIF file)
  
  suppressPackageStartupMessages({
    require(ggplot2)
    require(dplyr)
    require(magick)
  })
  
  cat(sprintf("Creating simple 3D GIF with %d frames...\n", n_frames))
  
  # Extract gene scores
  gene_scores <- as.data.frame(pca_result$ind$coord[, dims])
  colnames(gene_scores) <- c("PC1", "PC2", "PC3")
  gene_scores$Gene_name <- rownames(gene_scores)
  
  # Merge with metadata
  plot_data <- gene_scores %>%
    left_join(gene_data %>% dplyr::select(Gene_name, expression_group),
              by = "Gene_name")
  
  # Get variance explained
  var_exp <- pca_result$eig[dims, "percentage of variance"]
  
  # Prepare loadings if requested
  if (show_loadings) {
    loadings <- as.data.frame(pca_result$var$coord[, dims])
    colnames(loadings) <- c("PC1", "PC2", "PC3")
    loadings$Codon <- rownames(loadings)
    
    # Scale loadings
    loadings <- loadings %>%
      mutate(
        PC1 = PC1 * loading_scale,
        PC2 = PC2 * loading_scale,
        PC3 = PC3 * loading_scale
      )
    
    # Add classification if available
    if (!is.null(codon_test_results)) {
      loadings <- loadings %>%
        left_join(codon_test_results %>% 
                    dplyr::select(Codon, Classification, Significant, Difference),
                  by = "Codon")
    }
  }
  
  # Color scheme
  color_values <- c("Top 5%" = "#D62728", "Bottom 95%" = "#1F77B4")
  
  # Rotation function: rotate around Z-axis
  rotate_points <- function(x, y, z, angle_deg) {
    angle_rad <- angle_deg * pi / 180
    x_rot <- x * cos(angle_rad) - y * sin(angle_rad)
    y_rot <- x * sin(angle_rad) + y * cos(angle_rad)
    z_rot <- z
    return(data.frame(x = x_rot, y = y_rot, z = z_rot))
  }
  
  # Create frames
  angles <- seq(0, 360, length.out = n_frames + 1)[1:n_frames]
  temp_dir <- tempdir()
  frame_files <- character(n_frames)
  
  cat("Rendering frames...\n")
  
  for (i in seq_along(angles)) {
    angle <- angles[i]
    
    # Rotate gene data
    rotated_genes <- rotate_points(plot_data$PC1, plot_data$PC2, plot_data$PC3, angle)
    plot_data$x <- rotated_genes$x
    plot_data$y <- rotated_genes$z  # Use Z for vertical axis
    
    # Create plot
    p <- ggplot(plot_data, aes(x = x, y = y, color = expression_group)) +
      geom_point(size = point_size, alpha = 0.6) +
      scale_color_manual(values = color_values) +
      labs(
        title = sprintf("%s (Angle: %.0f°)", title, angle),
        subtitle = sprintf("PC%d vs PC%d (rotated)", dims[1], dims[3]),
        x = sprintf("Rotated Axis (%.1f%% + %.1f%% variance)", var_exp[1], var_exp[2]),
        y = sprintf("PC%d (%.1f%% variance)", dims[3], var_exp[3]),
        color = "Expression Group"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "gray70")
      ) +
      coord_fixed(ratio = 1)
    
    # Add loadings if requested
    if (show_loadings) {
      # Rotate loading arrows
      rotated_loadings <- rotate_points(loadings$PC1, loadings$PC2, loadings$PC3, angle)
      loadings$x <- rotated_loadings$x
      loadings$y <- rotated_loadings$z
      
      # Filter to most important loadings (top 20 by magnitude)
      loadings$magnitude <- sqrt(loadings$x^2 + loadings$y^2)
      top_loadings <- loadings %>%
        arrange(desc(magnitude)) %>%
        head(20)
      
      # Color by significance
      if (!is.null(codon_test_results)) {
        top_loadings$arrow_color <- ifelse(
          !is.na(top_loadings$Significant) & top_loadings$Significant,
          ifelse(top_loadings$Difference > 0, "Enriched", "Depleted"),
          "Non-significant"
        )
      } else {
        top_loadings$arrow_color <- "Codon"
      }
      
      # Add arrows
      p <- p +
        geom_segment(
          data = top_loadings,
          aes(x = 0, y = 0, xend = x, yend = y, color = NULL, linetype = arrow_color),
          arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
          linewidth = 0.7,
          alpha = 0.7,
          inherit.aes = FALSE
        ) +
        geom_text(
          data = top_loadings,
          aes(x = x * 1.1, y = y * 1.1, label = Codon, color = NULL),
          size = 3,
          inherit.aes = FALSE,
          alpha = 0.8
        ) +
        scale_linetype_manual(
          name = "Codon Status",
          values = c("Enriched" = "solid", "Depleted" = "dashed", "Non-significant" = "dotted", "Codon" = "solid")
        )
    }
    
    # Save frame
    frame_file <- file.path(temp_dir, sprintf("frame_%04d.png", i))
    ggsave(frame_file, plot = p, width = width/resolution, height = height/resolution, 
           dpi = resolution, bg = "white")
    frame_files[i] <- frame_file
    
    if (i %% 10 == 0) cat(sprintf("  Frame %d/%d\n", i, n_frames))
  }
  
  # Combine frames into GIF
  cat("Creating GIF...\n")
  frames <- image_read(frame_files)
  
  # Animate with 25 fps (factor of 100)
  animation <- image_animate(frames, fps = 25)
  image_write(animation, output_file)
  
  # Clean up
  unlink(frame_files)
  
  cat(sprintf("✓ Simple 3D GIF saved: %s\n", output_file))
  cat(sprintf("  Duration: %.1f seconds\n", n_frames / 25))
  cat(sprintf("  Size: %d x %d pixels\n", width, height))
  
  return(invisible(NULL))
}
