neutrality_plot <- function(gc_content, output_file = "neutrality_plot.pdf")
{
  #' Create neutrality plot (GC12 vs GC3)
  #' 
  #' @description The neutrality plot helps distinguish between mutation pressure
  #' and selection. A strong correlation suggests mutation bias, while deviation
  #' suggests selection pressure. The slope indicates the relative influence of
  #' mutation vs selection.
  #' 
  #' @param gc_content Data frame with GC content metrics (from calculate_gc_content)
  #' @param output_file Output file path for the plot
  #' 
  #' @return ggplot object
  #' ___________________________________________________________________________
  
  library(ggplot2)
  library(data.table)
  
  # Remove any NA or infinite values
  plot_data <- gc_content[is.finite(gc_content$GC12) & is.finite(gc_content$GC3), ]
  
  # Calculate correlation
  if(nrow(plot_data) > 0)
  {
    cor_test <- cor.test(plot_data$GC12, plot_data$GC3, method = "pearson")
    cor_val <- cor_test$estimate
    p_val <- cor_test$p.value
    
    # Fit linear regression
    lm_fit <- lm(GC3 ~ GC12, data = plot_data)
    slope <- coef(lm_fit)[2]
    intercept <- coef(lm_fit)[1]
    
    # Create plot
    p <- ggplot(plot_data, aes(x = GC12, y = GC3)) +
      geom_point(alpha = 0.3, size = 1) +
      geom_smooth(method = "lm", color = "red", se = TRUE) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
      theme_minimal() +
      labs(title = "Neutrality Plot (GC12 vs GC3)",
           subtitle = sprintf("r = %.3f, p < %.2e, slope = %.3f", 
                            cor_val, p_val, slope),
           x = "GC12 (1st and 2nd codon positions)",
           y = "GC3 (3rd codon position)") +
      annotate("text", x = min(plot_data$GC12) + 0.05, 
               y = max(plot_data$GC3) - 0.05,
               label = paste("Slope:", round(slope, 3), 
                           "\nInterpretation:",
                           ifelse(slope > 0.5, "Mutation dominates", 
                                  "Selection dominates")),
               hjust = 0, size = 3.5)
    
    ggsave(output_file, p, width = 8, height = 6)
    message(paste("Neutrality plot saved to:", output_file))
    
    return(p)
  }
  else
  {
    stop("No valid data for neutrality plot")
  }
}


enc_plot <- function(enc_values, gc_content, output_file = "enc_plot.pdf")
{
  #' Create ENC plot (ENC vs GC3s)
  #' 
  #' @description The ENC plot helps identify genes under selection for codon
  #' usage. The expected ENC curve represents mutation-drift equilibrium. Genes
  #' below the curve are under selection for codon bias.
  #' 
  #' @param enc_values Data frame with ENC values (from calculate_enc)
  #' @param gc_content Data frame with GC content metrics (from calculate_gc_content)
  #' @param output_file Output file path for the plot
  #' 
  #' @return ggplot object
  #' ___________________________________________________________________________
  
  library(ggplot2)
  library(data.table)
  
  # Merge data
  plot_data <- merge(enc_values, gc_content[, c("Gene_name", "GC3s")], 
                     by = "Gene_name")
  
  # Remove invalid values
  plot_data <- plot_data[is.finite(plot_data$ENC) & 
                         is.finite(plot_data$GC3s) &
                         plot_data$ENC > 0 & plot_data$ENC <= 61, ]
  
  # Calculate expected ENC under mutation-drift equilibrium (Wright 1990)
  # ENC_expected = 2 + GC3s + 29/(GC3s^2 + (1-GC3s)^2)
  gc3s_range <- seq(0, 1, by = 0.01)
  enc_expected <- 2 + gc3s_range + 29 / (gc3s_range^2 + (1 - gc3s_range)^2)
  
  expected_curve <- data.frame(
    GC3s = gc3s_range,
    ENC_expected = enc_expected
  )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = GC3s, y = ENC)) +
    geom_point(alpha = 0.3, size = 1, color = "darkgray") +
    geom_line(data = expected_curve, aes(x = GC3s, y = ENC_expected), 
              color = "red", linewidth = 1) +
    theme_minimal() +
    labs(title = "ENC Plot (Effective Number of Codons vs GC3s)",
         subtitle = "Red curve: expected ENC under mutation-drift equilibrium",
         x = "GC3s (GC content at synonymous 3rd positions)",
         y = "ENC (Effective Number of Codons)") +
    ylim(20, 61) +
    xlim(0, 1) +
    annotate("text", x = 0.1, y = 25,
             label = "Genes below curve:\nunder selection for codon bias",
             hjust = 0, size = 3.5, color = "red")
  
  ggsave(output_file, p, width = 8, height = 6)
  message(paste("ENC plot saved to:", output_file))
  
  return(p)
}


pr2_bias_plot <- function(codon_counts, output_file = "pr2_plot.pdf")
{
  #' Create PR2 bias plot (Parity Rule 2 analysis)
  #' 
  #' @description PR2 plot analyzes the bias between purines (A, G) and 
  #' pyrimidines (T, C) at the 3rd codon position. Center (0.5, 0.5) indicates
  #' no bias. Deviation suggests selection or mutation pressure.
  #' 
  #' @param codon_counts Data table with codon counts per gene
  #' @param output_file Output file path for the plot
  #' 
  #' @return ggplot object
  #' ___________________________________________________________________________
  
  library(ggplot2)
  library(data.table)
  
  codon_cols <- setdiff(names(codon_counts), "Gene_name")
  
  results <- data.frame(
    Gene_name = codon_counts$Gene_name,
    AU3 = numeric(nrow(codon_counts)),  # A/(A+T) at 3rd position
    GC3 = numeric(nrow(codon_counts))   # G/(G+C) at 3rd position
  )
  
  for(gene_idx in 1:nrow(codon_counts))
  {
    a3 <- 0
    t3 <- 0
    g3 <- 0
    c3 <- 0
    
    for(codon in codon_cols)
    {
      count <- as.numeric(codon_counts[gene_idx, codon, with = FALSE])
      
      if(count > 0)
      {
        base3 <- substr(codon, 3, 3)
        
        if(base3 == "A") a3 <- a3 + count
        else if(base3 == "T") t3 <- t3 + count
        else if(base3 == "G") g3 <- g3 + count
        else if(base3 == "C") c3 <- c3 + count
      }
    }
    
    # Calculate ratios
    at_total <- a3 + t3
    gc_total <- g3 + c3
    
    results$AU3[gene_idx] <- ifelse(at_total > 0, a3 / at_total, 0.5)
    results$GC3[gene_idx] <- ifelse(gc_total > 0, g3 / gc_total, 0.5)
  }
  
  # Remove invalid values
  plot_data <- results[is.finite(results$AU3) & is.finite(results$GC3), ]
  
  # Create plot
  p <- ggplot(plot_data, aes(x = AU3, y = GC3)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(title = "PR2 Bias Plot (Parity Rule 2)",
         subtitle = "Analysis of purine/pyrimidine bias at 3rd codon position",
         x = "A3/(A3+T3) - A vs T at 3rd position",
         y = "G3/(G3+C3) - G vs C at 3rd position") +
    xlim(0, 1) +
    ylim(0, 1) +
    annotate("text", x = 0.5, y = 0.5, label = "No bias",
             size = 3, color = "red", vjust = -0.5)
  
  ggsave(output_file, p, width = 8, height = 6)
  message(paste("PR2 plot saved to:", output_file))
  
  return(p)
}
