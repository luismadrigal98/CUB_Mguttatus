visualize_codon_usage <- function(codon_counts, genetic_code, 
                                  output_file = "codon_usage_heatmap.pdf",
                                  type = "heatmap")
{
  #' Visualize codon usage patterns
  #' 
  #' @description Creates visualizations of codon usage bias patterns.
  #' Supports heatmap and logo-style representations.
  #' 
  #' @param codon_counts Data table with codon counts per gene
  #' @param genetic_code Named vector mapping codons to amino acids
  #' @param output_file Output file path for the plot
  #' @param type Type of visualization: "heatmap" or "barplot"
  #' 
  #' @return Creates a plot file and returns invisible
  #' ___________________________________________________________________________
  
  library(data.table)
  library(ggplot2)
  
  # Calculate genome-wide codon frequencies
  codon_cols <- setdiff(names(codon_counts), "Gene_name")
  
  # Sum across all genes
  total_counts <- colSums(codon_counts[, codon_cols, with = FALSE])
  
  # Create data frame for plotting
  plot_data <- data.frame(
    Codon = names(total_counts),
    Count = as.numeric(total_counts),
    AA = genetic_code[names(total_counts)],
    stringsAsFactors = FALSE
  )
  
  # Remove STOP codons
  plot_data <- plot_data[plot_data$AA != "STOP", ]
  
  # Calculate frequencies within each amino acid
  plot_data <- as.data.table(plot_data)
  plot_data[, Frequency := Count / sum(Count), by = AA]
  plot_data[, RSCU := Frequency / (1 / .N), by = AA]
  
  if(type == "heatmap")
  {
    # Order amino acids and codons
    plot_data <- plot_data[order(AA, -RSCU)]
    plot_data$Codon <- factor(plot_data$Codon, levels = unique(plot_data$Codon))
    
    p <- ggplot(plot_data, aes(x = Codon, y = AA, fill = RSCU)) +
      geom_tile(color = "white", linewidth = 0.5) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                          midpoint = 1, name = "RSCU") +
      theme_custom() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
            axis.text.y = element_text(size = 10)) +
      labs(title = "Genome-wide Codon Usage Bias",
           x = "Codon", y = "Amino Acid") +
      facet_grid(. ~ AA, scales = "free_x", space = "free_x")
    
    ggsave(output_file, p, width = 16, height = 8)
  }
  else if(type == "barplot")
  {
    # Create barplot for each amino acid
    p <- ggplot(plot_data, aes(x = Codon, y = RSCU, fill = AA)) +
      geom_bar(stat = "identity") +
      facet_wrap(~ AA, scales = "free_x", ncol = 5) +
      theme_custom() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
            legend.position = "none") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      labs(title = "Relative Synonymous Codon Usage (RSCU) by Amino Acid",
           x = "Codon", y = "RSCU")
    
    ggsave(output_file, p, width = 14, height = 10)
  }
  
  message(paste("Plot saved to:", output_file))
  return(invisible(plot_data))
}

create_codon_logo <- function(codon_counts, genetic_code, amino_acid,
                              output_file = NULL)
{
  #' Create sequence logo-style visualization for a specific amino acid
  #' 
  #' @description Creates a logo-style visualization showing nucleotide
  #' preferences at each codon position for a given amino acid.
  #' 
  #' @param codon_counts Data table with codon counts per gene
  #' @param genetic_code Named vector mapping codons to amino acids
  #' @param amino_acid Three-letter amino acid code (e.g., "Leu", "Ser")
  #' @param output_file Optional output file path
  #' 
  #' @return ggplot object
  #' ___________________________________________________________________________
  
  library(data.table)
  library(ggplot2)
  
  # Get codons for this amino acid
  aa_codons <- names(genetic_code)[genetic_code == amino_acid]
  
  if(length(aa_codons) == 0)
  {
    stop(paste("No codons found for amino acid:", amino_acid))
  }
  
  # Get counts
  codon_cols <- intersect(aa_codons, names(codon_counts))
  total_counts <- colSums(codon_counts[, codon_cols, with = FALSE])
  
  # Parse codon positions
  logo_data <- data.frame()
  
  for(codon in names(total_counts))
  {
    bases <- strsplit(codon, "")[[1]]
    count <- total_counts[codon]
    
    for(pos in 1:3)
    {
      logo_data <- rbind(logo_data, data.frame(
        Position = pos,
        Base = bases[pos],
        Count = count
      ))
    }
  }
  
  # Calculate frequencies per position
  logo_data <- as.data.table(logo_data)
  logo_data[, Frequency := Count / sum(Count), by = Position]
  
  # Calculate information content (bits)
  logo_data[, IC := -Frequency * log2(Frequency)]
  logo_data[IC < 0 | is.na(IC) | is.infinite(IC), IC := 0]
  
  # Calculate height (frequency * IC)
  max_ic_per_pos <- logo_data[, .(MaxIC = sum(IC)), by = Position]
  logo_data <- merge(logo_data, max_ic_per_pos, by = "Position")
  logo_data[, Height := Frequency * 2]  # Scale for visualization
  
  # Order bases by frequency within each position
  logo_data <- logo_data[order(Position, -Frequency)]
  
  p <- ggplot(logo_data, aes(x = Position, y = Height, fill = Base)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    scale_fill_manual(values = c("A" = "#109648", "T" = "#F00000", 
                                  "G" = "#F59E00", "C" = "#255FBE")) +
    theme_custom() +
    labs(title = paste("Codon Usage Pattern for", amino_acid),
         x = "Codon Position", y = "Frequency") +
    scale_x_continuous(breaks = 1:3) +
    theme(panel.grid.minor = element_blank())
  
  if(!is.null(output_file))
  {
    ggsave(output_file, p, width = 6, height = 4)
    message(paste("Logo saved to:", output_file))
  }
  
  return(p)
}
