analyze_codon_loadings <- function(analysis_result, 
                                    analysis_type = "PCA",
                                    dims = c(1, 2),
                                    genetic_code,
                                    output_file = NULL)
{
  #' Analyze and visualize codon loading patterns
  #' 
  #' @description Categorizes codons by their 3rd position nucleotide
  #' and shows how AT-ending vs GC-ending codons load on dimensions.
  #' Useful for understanding mutational bias patterns.
  #' 
  #' @param analysis_result Result from PCA() or CA() function
  #' @param analysis_type "PCA" or "CA"
  #' @param dims Vector of dimension numbers to analyze
  #' @param genetic_code Named vector mapping codons to amino acids
  #' @param output_file Path to save plot (NULL for display only)
  #' 
  #' @return List with loading_data and plot
  #' ___________________________________________________________________________
  
  require(ggplot2)
  require(dplyr)
  
  # Extract loadings
  if (analysis_type == "PCA") {
    loadings <- analysis_result$var$coord[, dims, drop = FALSE]
    dim_labels <- paste0("PC", dims)
  } else if (analysis_type == "CA") {
    loadings <- analysis_result$col$coord[, dims, drop = FALSE]
    dim_labels <- paste0("CA", dims)
  } else {
    stop("analysis_type must be 'PCA' or 'CA'")
  }
  
  # Create data frame with codon information
  loading_df <- data.frame(
    Codon = rownames(loadings),
    loadings,
    stringsAsFactors = FALSE
  )
  
  colnames(loading_df)[2:(length(dims)+1)] <- paste0("Dim", dims)
  
  # Add nucleotide information
  loading_df <- loading_df %>%
    mutate(
      Third_Position = substr(Codon, 3, 3),
      GC_Type = case_when(
        Third_Position %in% c("G", "C") ~ "GC-ending",
        Third_Position %in% c("A", "T") ~ "AT-ending",
        TRUE ~ "Other"
      ),
      AA = genetic_code[Codon]
    ) %>%
    filter(AA != "STOP")  # Remove stop codons
  
  # Calculate mean loadings by GC type
  gc_summary <- loading_df %>%
    group_by(GC_Type) %>%
    summarise(
      across(starts_with("Dim"), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"),
      n = n()
    )
  
  cat("\n=== Codon Loading Analysis ===\n")
  cat(sprintf("%s Analysis - Dimensions: %s\n", 
              analysis_type, paste(dim_labels, collapse = ", ")))
  cat("\nMean loadings by codon type:\n")
  print(gc_summary)
  
  # Create visualization for first two dimensions
  if (length(dims) >= 2) {
    dim1 <- paste0("Dim", dims[1])
    dim2 <- paste0("Dim", dims[2])
    
    p <- ggplot(loading_df, aes(x = .data[[dim1]], y = .data[[dim2]], 
                                color = .data[["GC_Type"]])) +
      geom_point(alpha = 0.6, size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
      # Add arrows showing mean direction for each group
      geom_segment(data = gc_summary,
                   aes(x = 0, y = 0, 
                       xend = get(paste0(dim1, "_mean")), 
                       yend = get(paste0(dim2, "_mean")),
                       color = GC_Type),
                   arrow = arrow(length = unit(0.5, "cm"), type = "closed"),
                   linewidth = 1.5, alpha = 0.8) +
      scale_color_manual(values = c("AT-ending" = "#E41A1C",
                                    "GC-ending" = "#377EB8")) +
      labs(title = paste(analysis_type, "Codon Loadings: AT-ending vs GC-ending"),
           subtitle = "Arrows show mean direction for each group",
           x = paste(analysis_type, "Dimension", dims[1]),
           y = paste(analysis_type, "Dimension", dims[2]),
           color = "Codon Type (3rd position)") +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right"
      )
    
    if (!is.null(output_file)) {
      ggsave(output_file, plot = p, width = 10, height = 8, dpi = 300)
      message(paste("Codon loading plot saved to:", output_file))
    }
    
    # Test for significant difference in loadings
    for (d in dims) {
      dim_name <- paste0("Dim", d)
      at_loadings <- loading_df[[dim_name]][loading_df$GC_Type == "AT-ending"]
      gc_loadings <- loading_df[[dim_name]][loading_df$GC_Type == "GC-ending"]
      
      wtest <- wilcox.test(at_loadings, gc_loadings)
      cat(sprintf("\n%s Dimension %d: AT vs GC loadings\n", analysis_type, d))
      cat(sprintf("  AT-ending mean: %.4f\n", mean(at_loadings)))
      cat(sprintf("  GC-ending mean: %.4f\n", mean(gc_loadings)))
      cat(sprintf("  Wilcoxon p-value: %.4e %s\n", 
                  wtest$p.value,
                  ifelse(wtest$p.value < 0.05, "***", "")))
    }
    
    return(list(loading_data = loading_df, plot = p, summary = gc_summary))
  }
  
  return(list(loading_data = loading_df, summary = gc_summary))
}
