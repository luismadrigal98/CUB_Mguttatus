create_3d_pca_plot <- function(pca_result,
                              gene_data,
                              codon_test_results = NULL,
                              preferred_codons = NULL,
                              dims = c(1, 2, 3),
                              color_by = "expression",
                              show_loadings = TRUE,
                              loading_scale = 1.0,
                              title = "3D PCA Plot") {
  #' Create interactive 3D PCA plot with gene scores and codon loadings
  #' 
  #' @param pca_result PCA result from FactoMineR::PCA()
  #' @param gene_data Data frame with gene metadata (Gene_name, expression_group, etc.)
  #' @param codon_test_results Optional: results from test_codon_proportions()
  #' @param preferred_codons Optional: preferred codons data frame
  #' @param dims Vector of 3 dimensions to plot
  #' @param color_by How to color points: "expression", "selection"
  #' @param show_loadings Show codon loading vectors?
  #' @param loading_scale Scaling factor for loading arrows
  #' @param title Plot title
  #' 
  #' @return plotly object
  
  suppressPackageStartupMessages({
    require(plotly)
    require(dplyr)
  })
  
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
  
  # Create axis labels
  xlab <- sprintf("PC%d (%.1f%%)", dims[1], var_exp[1])
  ylab <- sprintf("PC%d (%.1f%%)", dims[2], var_exp[2])
  zlab <- sprintf("PC%d (%.1f%%)", dims[3], var_exp[3])
  
  # Color scheme
  if (color_by == "expression") {
    colors <- c("Top 5%" = "#D62728", "Bottom 95%" = "#1F77B4")
  } else {
    colors <- c("Selected" = "#2CA02C", "Neutral" = "#FF7F0E")
  }
  
  # Create 3D scatter plot for genes
  p <- plot_ly() %>%
    add_trace(
      data = plot_data,
      x = ~PC1, y = ~PC2, z = ~PC3,
      type = "scatter3d",
      mode = "markers",
      color = ~expression_group,
      colors = colors,
      marker = list(size = 3, opacity = 0.6),
      text = ~paste0("Gene: ", Gene_name, "<br>Group: ", expression_group),
      hoverinfo = "text",
      name = "Genes"
    )
  
  # Add codon loadings if requested
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
    
    # Add preferred codon info if available
    if (!is.null(preferred_codons)) {
      codon_col <- if("Codon" %in% colnames(preferred_codons)) "Codon" else "codon"
      loadings <- loadings %>%
        left_join(preferred_codons %>% 
                    dplyr::select(!!sym(codon_col), relative_adaptiveness),
                  by = c("Codon" = codon_col))
      
      loadings <- loadings %>%
        mutate(is_preferred = !is.na(relative_adaptiveness) & relative_adaptiveness == 1.0)
    }
    
    # Add loading vectors as cones/arrows
    for (i in 1:nrow(loadings)) {
      codon_color <- if (!is.null(codon_test_results) && !is.na(loadings$Significant[i]) && loadings$Significant[i]) {
        if (loadings$Difference[i] > 0) "#2CA02C" else "#FF7F0E"  # Green for enriched, orange for depleted
      } else {
        "#999999"  # Gray for non-significant
      }
      
      p <- p %>%
        add_trace(
          x = c(0, loadings$PC1[i]),
          y = c(0, loadings$PC2[i]),
          z = c(0, loadings$PC3[i]),
          type = "scatter3d",
          mode = "lines+text",
          line = list(color = codon_color, width = 3),
          text = c("", loadings$Codon[i]),
          textposition = "top center",
          textfont = list(size = 10, color = codon_color),
          hoverinfo = "text",
          hovertext = paste0(
            "Codon: ", loadings$Codon[i], "<br>",
            "PC", dims[1], ": ", round(loadings$PC1[i], 3), "<br>",
            "PC", dims[2], ": ", round(loadings$PC2[i], 3), "<br>",
            "PC", dims[3], ": ", round(loadings$PC3[i], 3),
            if (!is.null(codon_test_results) && !is.na(loadings$Classification[i])) 
              paste0("<br>Class: ", loadings$Classification[i]) else ""
          ),
          showlegend = FALSE,
          name = loadings$Codon[i]
        )
    }
  }
  
  # Layout
  p <- p %>%
    layout(
      title = title,
      scene = list(
        xaxis = list(title = xlab),
        yaxis = list(title = ylab),
        zaxis = list(title = zlab),
        camera = list(
          eye = list(x = 1.5, y = 1.5, z = 1.5)
        )
      ),
      showlegend = TRUE,
      legend = list(x = 0.02, y = 0.98)
    )
  
  return(p)
}


create_3d_pca_animation <- function(pca_result,
                                    gene_data,
                                    codon_test_results = NULL,
                                    preferred_codons = NULL,
                                    dims = c(1, 2, 3),
                                    color_by = "expression",
                                    show_loadings = TRUE,
                                    loading_scale = 1.0,
                                    title = "3D PCA Animation",
                                    output_file = "pca_3d_animation.html",
                                    n_frames = 360,
                                    frame_duration = 50) {
  #' Create rotating 3D PCA animation and save as HTML
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
  #' @param output_file Output HTML file path
  #' @param n_frames Number of frames in rotation (default: 360 for 1 degree steps)
  #' @param frame_duration Duration of each frame in ms (default: 50ms)
  #' 
  #' @return plotly object (also saves HTML file)
  
  suppressPackageStartupMessages({
    require(plotly)
    require(dplyr)
    require(htmlwidgets)
  })
  
  cat(sprintf("Creating 3D animation with %d frames...\n", n_frames))
  
  # Create base plot
  p <- create_3d_pca_plot(
    pca_result = pca_result,
    gene_data = gene_data,
    codon_test_results = codon_test_results,
    preferred_codons = preferred_codons,
    dims = dims,
    color_by = color_by,
    show_loadings = show_loadings,
    loading_scale = loading_scale,
    title = title
  )
  
  # Add animation frames (rotation around Z-axis)
  angles <- seq(0, 360, length.out = n_frames)
  
  # Generate camera positions for rotation
  frames_list <- lapply(1:length(angles), function(i) {
    angle <- angles[i]
    rad <- angle * pi / 180
    
    list(
      name = as.character(i),
      layout = list(
        scene = list(
          camera = list(
            eye = list(
              x = 1.5 * cos(rad),
              y = 1.5 * sin(rad),
              z = 1.5
            ),
            center = list(x = 0, y = 0, z = 0),
            up = list(x = 0, y = 0, z = 1)
          )
        )
      )
    )
  })
  
  # Add frames to plot
  p$x$frames <- frames_list
  
  # Add animation configuration
  p$x$layout$updatemenus <- list(
    list(
      type = "buttons",
      direction = "left",
      showactive = FALSE,
      x = 0.1,
      y = 0,
      xanchor = "right",
      yanchor = "top",
      pad = list(t = 87, r = 10),
      buttons = list(
        list(
          label = "Play",
          method = "animate",
          args = list(
            NULL,
            list(
              frame = list(duration = frame_duration, redraw = FALSE),
              fromcurrent = TRUE,
              transition = list(duration = 0)
            )
          )
        ),
        list(
          label = "Pause",
          method = "animate",
          args = list(
            list(NULL),
            list(
              frame = list(duration = 0, redraw = FALSE),
              mode = "immediate",
              transition = list(duration = 0)
            )
          )
        )
      )
    )
  )
  
  # Add slider
  p$x$layout$sliders <- list(
    list(
      active = 0,
      steps = lapply(1:length(angles), function(i) {
        list(
          label = sprintf("%.0f°", angles[i]),
          method = "animate",
          args = list(
            list(as.character(i)),
            list(
              frame = list(duration = frame_duration, redraw = FALSE),
              mode = "immediate",
              transition = list(duration = 0)
            )
          )
        )
      }),
      x = 0.1,
      y = 0,
      len = 0.9,
      xanchor = "left",
      yanchor = "top",
      pad = list(t = 50, b = 10),
      currentvalue = list(
        visible = TRUE,
        prefix = "Angle: ",
        xanchor = "right"
      )
    )
  )
  
  # Save as HTML
  htmlwidgets::saveWidget(
    widget = p,
    file = output_file,
    selfcontained = TRUE
  )
  
  cat(sprintf("✓ 3D animation saved: %s\n", output_file))
  cat("  Open in browser to view interactive rotating plot\n\n")
  
  return(p)
}


create_3d_pca_video <- function(pca_result,
                                gene_data,
                                codon_test_results = NULL,
                                preferred_codons = NULL,
                                dims = c(1, 2, 3),
                                color_by = "expression",
                                show_loadings = TRUE,
                                loading_scale = 1.0,
                                title = "3D PCA",
                                output_file = "pca_3d_video.gif",
                                n_frames = 120,
                                fps = 30,
                                width = 800,
                                height = 600) {
  #' Create rotating 3D PCA video as GIF
  #' 
  #' Requires: rgl, magick packages
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
  #' @param fps Frames per second
  #' @param width Video width in pixels
  #' @param height Video height in pixels
  #' 
  #' @return NULL (saves GIF file)
  
  suppressPackageStartupMessages({
    require(rgl)
    require(magick)
    require(dplyr)
  })
  
  cat(sprintf("Creating 3D video with %d frames at %d fps...\n", n_frames, fps))
  
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
  
  # Create axis labels
  xlab <- sprintf("PC%d (%.1f%%)", dims[1], var_exp[1])
  ylab <- sprintf("PC%d (%.1f%%)", dims[2], var_exp[2])
  zlab <- sprintf("PC%d (%.1f%%)", dims[3], var_exp[3])
  
  # Setup RGL device
  open3d(windowRect = c(0, 0, width, height))
  bg3d("white")
  
  # Color mapping
  color_map <- ifelse(plot_data$expression_group == "Top 5%", "#D62728", "#1F77B4")
  
  # Plot gene points
  plot3d(
    plot_data$PC1, plot_data$PC2, plot_data$PC3,
    col = color_map,
    size = 5,
    xlab = xlab, ylab = ylab, zlab = zlab,
    main = title,
    type = "s",
    radius = 0.02
  )
  
  # Add codon loadings if requested
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
    
    # Draw loading arrows
    for (i in 1:nrow(loadings)) {
      arrow_color <- if (!is.null(codon_test_results) && !is.na(loadings$Significant[i]) && loadings$Significant[i]) {
        if (loadings$Difference[i] > 0) "#2CA02C" else "#FF7F0E"
      } else {
        "#999999"
      }
      
      # Draw arrow
      arrow3d(
        c(0, 0, 0),
        c(loadings$PC1[i], loadings$PC2[i], loadings$PC3[i]),
        type = "rotation",
        col = arrow_color,
        barblen = 0.05
      )
      
      # Add text label
      text3d(
        loadings$PC1[i] * 1.1,
        loadings$PC2[i] * 1.1,
        loadings$PC3[i] * 1.1,
        loadings$Codon[i],
        col = arrow_color,
        cex = 0.8
      )
    }
  }
  
  # Create rotation animation
  cat("Rendering frames...\n")
  angles <- seq(0, 360, length.out = n_frames)
  
  # Temporary directory for frames
  temp_dir <- tempdir()
  frame_files <- character(n_frames)
  
  for (i in 1:n_frames) {
    view3d(theta = angles[i], phi = 20, fov = 60, zoom = 0.8)
    frame_file <- file.path(temp_dir, sprintf("frame_%04d.png", i))
    rgl.snapshot(frame_file)
    frame_files[i] <- frame_file
    
    if (i %% 10 == 0) cat(sprintf("  Frame %d/%d\n", i, n_frames))
  }
  
  # Close RGL device
  close3d()
  
  # Combine frames into GIF
  cat("Creating GIF...\n")
  frames <- image_read(frame_files)
  
  # magick requires fps to be a factor of 100 (e.g., 10, 20, 25, 50, 100)
  # Convert fps to delay (1/100ths of a second)
  delay <- round(100 / fps)
  
  animation <- image_animate(frames, fps = 100/delay)
  image_write(animation, output_file)
  
  # Clean up temporary files
  unlink(frame_files)
  
  cat(sprintf("✓ 3D video saved: %s\n", output_file))
  cat(sprintf("  Duration: %.1f seconds\n", n_frames / fps))
  
  return(invisible(NULL))
}
