plot_multivariate_analysis <- function(coord_data, 
                                       dims = c("Dim.1", "Dim.2", "Dim.3"),
                                       group_var = "Expression_Group",
                                       analysis_type = "PCA",
                                       plot_type = "bivariate_ellipse",
                                       confidence_level = 0.95,
                                       output_file = NULL,
                                       point_alpha = 0.6,
                                       point_size = 2,
                                       ellipse_alpha = 0.3,
                                       colors = NULL)
{
  #' Create publication-quality multivariate analysis plots
  #' 
  #' @description Creates various visualizations for PCA, CA, or other 
  #' multivariate analyses including bivariate ellipse plots and 3D scatter plots.
  #' 
  #' @param coord_data Data frame with dimension coordinates and grouping variable
  #' @param dims Character vector of dimension column names (length 2 or 3)
  #' @param group_var Column name for grouping variable (e.g., "Expression_Group")
  #' @param analysis_type String: "PCA", "CA", or custom label for titles
  #' @param plot_type One of: "bivariate_ellipse", "bivariate_pairs", 
  #'                  "3D_static", "3D_interactive"
  #' @param confidence_level Confidence level for ellipses (default 0.95)
  #' @param output_file Path to save plot (NULL for no save)
  #' @param point_alpha Transparency for points (0-1)
  #' @param point_size Size of points
  #' @param ellipse_alpha Transparency for ellipses (0-1)
  #' @param colors Named vector of colors for groups (NULL for default)
  #' 
  #' @return ggplot object or plotly object (for 3D interactive)
  #' ___________________________________________________________________________
  
  require(ggplot2)
  require(dplyr)
  
  # Input validation
  if (!all(dims %in% names(coord_data))) {
    stop("Dimension columns not found in coord_data")
  }
  
  if (!group_var %in% names(coord_data)) {
    stop(paste("Grouping variable", group_var, "not found in coord_data"))
  }
  
  # Set default colors if not provided
  if (is.null(colors)) {
    n_groups <- length(unique(coord_data[[group_var]]))
    colors <- scales::hue_pal()(n_groups)
    names(colors) <- unique(coord_data[[group_var]])
  }
  
  # Create plots based on type
  if (plot_type == "bivariate_ellipse") {
    # Bivariate plot with confidence ellipses
    if (length(dims) < 2) {
      stop("Need at least 2 dimensions for bivariate plot")
    }
    
    p <- ggplot(coord_data, aes_string(x = dims[1], y = dims[2], 
                                        color = group_var, fill = group_var)) +
      geom_point(alpha = point_alpha, size = point_size) +
      stat_ellipse(geom = "polygon", level = confidence_level, 
                   alpha = ellipse_alpha, show.legend = FALSE) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      labs(title = paste(analysis_type, "- Bivariate Plot with", 
                        confidence_level * 100, "% Confidence Ellipses"),
           x = paste(analysis_type, "Dimension", gsub("[^0-9]", "", dims[1])),
           y = paste(analysis_type, "Dimension", gsub("[^0-9]", "", dims[2])),
           color = group_var,
           fill = group_var) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right",
        panel.grid.minor = element_blank()
      )
    
  } else if (plot_type == "bivariate_pairs") {
    # Pairwise bivariate plots with ellipses
    if (length(dims) < 3) {
      stop("Need at least 3 dimensions for pairwise bivariate plots")
    }
    
    require(GGally)
    
    # Select relevant columns
    plot_data <- coord_data[, c(dims, group_var)]
    
    # Create pairs plot
    p <- GGally::ggpairs(
      plot_data,
      columns = 1:length(dims),
      mapping = aes_string(color = group_var, fill = group_var),
      upper = list(continuous = "points"),
      lower = list(continuous = wrap("points", alpha = point_alpha)),
      diag = list(continuous = wrap("densityDiag", alpha = 0.5)),
      title = paste(analysis_type, "- Pairwise Bivariate Plots")
    ) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      theme_minimal()
    
  } else if (plot_type == "3D_static") {
    # 3D plot using scatterplot3d
    if (length(dims) < 3) {
      stop("Need 3 dimensions for 3D plot")
    }
    
    require(scatterplot3d)
    
    # Convert group to numeric for colors
    group_numeric <- as.numeric(as.factor(coord_data[[group_var]]))
    plot_colors <- colors[as.factor(coord_data[[group_var]])]
    
    # Create 3D scatter plot
    if (!is.null(output_file)) {
      pdf(output_file, width = 10, height = 8)
    }
    
    s3d <- scatterplot3d(
      x = coord_data[[dims[1]]],
      y = coord_data[[dims[2]]],
      z = coord_data[[dims[3]]],
      color = plot_colors,
      pch = 19,
      cex.symbols = 1.5,
      angle = 55,
      main = paste(analysis_type, "- 3D Scatter Plot"),
      xlab = paste(analysis_type, "Dim", gsub("[^0-9]", "", dims[1])),
      ylab = paste(analysis_type, "Dim", gsub("[^0-9]", "", dims[2])),
      zlab = paste(analysis_type, "Dim", gsub("[^0-9]", "", dims[3])),
      grid = TRUE,
      box = TRUE
    )
    
    # Add legend
    legend("topright",
           legend = names(colors),
           col = colors,
           pch = 19,
           cex = 0.8,
           title = group_var,
           bty = "n")
    
    if (!is.null(output_file)) {
      dev.off()
      message(paste("3D static plot saved to:", output_file))
    }
    
    return(invisible(s3d))
    
  } else if (plot_type == "3D_interactive") {
    # Interactive 3D plot using plotly
    if (length(dims) < 3) {
      stop("Need 3 dimensions for 3D plot")
    }
    
    require(plotly)
    
    # Create interactive 3D scatter plot
    p <- plot_ly(
      data = coord_data,
      x = ~get(dims[1]),
      y = ~get(dims[2]),
      z = ~get(dims[3]),
      color = ~get(group_var),
      colors = colors,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 3, opacity = point_alpha)
    ) %>%
      layout(
        title = paste(analysis_type, "- Interactive 3D Scatter Plot"),
        scene = list(
          xaxis = list(title = paste(analysis_type, "Dim", 
                                     gsub("[^0-9]", "", dims[1]))),
          yaxis = list(title = paste(analysis_type, "Dim", 
                                     gsub("[^0-9]", "", dims[2]))),
          zaxis = list(title = paste(analysis_type, "Dim", 
                                     gsub("[^0-9]", "", dims[3])))
        )
      )
    
    if (!is.null(output_file)) {
      htmlwidgets::saveWidget(p, output_file)
      message(paste("Interactive 3D plot saved to:", output_file))
    }
    
    return(p)
    
  } else {
    stop(paste("Unknown plot_type:", plot_type))
  }
  
  # Save non-3D plots if output file specified
  if (!is.null(output_file) && plot_type %in% c("bivariate_ellipse", "bivariate_pairs")) {
    ggsave(output_file, plot = p, width = 10, height = 8, dpi = 300)
    message(paste("Plot saved to:", output_file))
  }
  
  return(p)
}


plot_variance_explained <- function(analysis_result, 
                                   analysis_type = "PCA",
                                   n_dims = 10,
                                   output_file = NULL)
{
  #' Plot variance explained by each dimension
  #' 
  #' @description Creates a bar plot showing the variance explained by
  #' each dimension in PCA or CA analysis
  #' 
  #' @param analysis_result Result from PCA() or CA() function
  #' @param analysis_type "PCA" or "CA"
  #' @param n_dims Number of dimensions to plot
  #' @param output_file Path to save plot (NULL for no save)
  #' 
  #' @return ggplot object
  #' ___________________________________________________________________________
  
  require(ggplot2)
  
  # Extract eigenvalues/variance explained
  if (analysis_type == "PCA") {
    variance <- analysis_result$eig[1:n_dims, 2]  # Percentage of variance
    dim_names <- paste0("PC", 1:n_dims)
  } else if (analysis_type == "CA") {
    variance <- analysis_result$eig[1:n_dims, 2]  # Percentage of variance
    dim_names <- paste0("CA", 1:n_dims)
  } else {
    stop("analysis_type must be 'PCA' or 'CA'")
  }
  
  # Create data frame
  var_df <- data.frame(
    Dimension = factor(dim_names, levels = dim_names),
    Variance_Explained = variance
  )
  
  # Create plot
  p <- ggplot(var_df, aes(x = Dimension, y = Variance_Explained)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.1f%%", Variance_Explained)), 
              vjust = -0.5, size = 3) +
    labs(title = paste(analysis_type, "- Variance Explained by Each Dimension"),
         x = "Dimension",
         y = "Percentage of Variance Explained (%)") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  if (!is.null(output_file)) {
    ggsave(output_file, plot = p, width = 10, height = 6, dpi = 300)
    message(paste("Variance plot saved to:", output_file))
  }
  
  return(p)
}


create_biplot <- function(analysis_result,
                         coord_data,
                         group_var = "Expression_Group",
                         analysis_type = "PCA",
                         dims = c(1, 2),
                         n_loadings = 10,
                         colors = NULL,
                         output_file = NULL)
{
  #' Create a biplot showing both observations and variable loadings
  #' 
  #' @description Creates a biplot for PCA or CA showing the relationship
  #' between observations (genes) and variables (codons)
  #' 
  #' @param analysis_result Result from PCA() or CA() function
  #' @param coord_data Data frame with coordinates and grouping variable
  #' @param group_var Column name for grouping variable
  #' @param analysis_type "PCA" or "CA"
  #' @param dims Vector of 2 dimension numbers to plot
  #' @param n_loadings Number of top variable loadings to show
  #' @param colors Named vector of colors for groups
  #' @param output_file Path to save plot (NULL for no save)
  #' 
  #' @return ggplot object
  #' ___________________________________________________________________________
  
  require(ggplot2)
  require(ggrepel)
  
  # Extract loadings/column coordinates
  if (analysis_type == "PCA") {
    loadings <- analysis_result$var$coord[, dims]
    dim_labels <- paste0("PC", dims)
  } else if (analysis_type == "CA") {
    loadings <- analysis_result$col$coord[, dims]
    dim_labels <- paste0("CA", dims)
  } else {
    stop("analysis_type must be 'PCA' or 'CA'")
  }
  
  # Calculate loading magnitudes and select top N
  loading_magnitude <- sqrt(rowSums(loadings^2))
  top_indices <- order(loading_magnitude, decreasing = TRUE)[1:n_loadings]
  loadings_top <- loadings[top_indices, ]
  
  # Scale loadings for visualization
  scale_factor <- max(abs(coord_data[[paste0("Dim.", dims[1])]])) / 
                  max(abs(loadings_top[, 1])) * 0.8
  loadings_scaled <- loadings_top * scale_factor
  
  loadings_df <- data.frame(
    Variable = rownames(loadings_top),
    Dim1 = loadings_scaled[, 1],
    Dim2 = loadings_scaled[, 2]
  )
  
  # Set default colors
  if (is.null(colors)) {
    n_groups <- length(unique(coord_data[[group_var]]))
    colors <- scales::hue_pal()(n_groups)
    names(colors) <- unique(coord_data[[group_var]])
  }
  
  # Create biplot
  dim_cols <- paste0("Dim.", dims)
  
  p <- ggplot() +
    # Plot observations
    geom_point(data = coord_data,
               aes_string(x = dim_cols[1], y = dim_cols[2], 
                         color = group_var),
               alpha = 0.4, size = 1.5) +
    # Plot loading arrows
    geom_segment(data = loadings_df,
                aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
                arrow = arrow(length = unit(0.3, "cm")),
                color = "red", alpha = 0.7, linewidth = 0.5) +
    # Label loadings
    geom_text_repel(data = loadings_df,
                   aes(x = Dim1, y = Dim2, label = Variable),
                   color = "darkred", size = 3, 
                   box.padding = 0.5, max.overlaps = 20) +
    scale_color_manual(values = colors) +
    labs(title = paste(analysis_type, "Biplot -", 
                      paste(dim_labels, collapse = " vs ")),
         x = paste(analysis_type, "Dimension", dims[1]),
         y = paste(analysis_type, "Dimension", dims[2]),
         color = group_var) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  if (!is.null(output_file)) {
    ggsave(output_file, plot = p, width = 12, height = 10, dpi = 300)
    message(paste("Biplot saved to:", output_file))
  }
  
  return(p)
}
