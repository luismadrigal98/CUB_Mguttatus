cub_summary <- function(codon_counts, genetic_code, output_dir = "./results")
{
  #' Comprehensive CUB analysis summary
  #' 
  #' @description Performs complete codon usage bias analysis including:
  #' - RSCU calculation
  #' - ENC calculation
  #' - GC content metrics
  #' - Creates all standard plots
  #' - Generates summary statistics
  #' 
  #' @param codon_counts Data table with codon counts per gene
  #' @param genetic_code Named vector mapping codons to amino acids
  #' @param output_dir Directory for output files
  #' 
  #' @return List with all analysis results
  #' ___________________________________________________________________________
  
  library(data.table)
  
  # Create output directory if needed
  if(!dir.exists(output_dir))
  {
    dir.create(output_dir, recursive = TRUE)
  }
  
  message("Starting comprehensive CUB analysis...")
  
  # 1. Calculate RSCU
  message("Calculating RSCU...")
  rscu_results <- calculate_rscu(codon_counts, genetic_code)
  
  # 2. Calculate ENC
  message("Calculating ENC...")
  enc_results <- calculate_enc(codon_counts, genetic_code)
  
  # 3. Calculate GC content
  message("Calculating GC content metrics...")
  gc_results <- calculate_gc_content(codon_counts)
  
  # 4. Merge all results
  message("Merging results...")
  summary_table <- merge(codon_counts, enc_results, by = "Gene_name")
  summary_table <- merge(summary_table, gc_results, by = "Gene_name")
  
  # 5. Generate visualizations
  message("Creating visualizations...")
  
  # Codon usage heatmap
  visualize_codon_usage(codon_counts, genetic_code, 
                       file.path(output_dir, "codon_usage_heatmap.pdf"),
                       type = "heatmap")
  
  # Codon usage barplot
  visualize_codon_usage(codon_counts, genetic_code, 
                       file.path(output_dir, "codon_usage_barplot.pdf"),
                       type = "barplot")
  
  # Neutrality plot
  neutrality_plot(gc_results, file.path(output_dir, "neutrality_plot.pdf"))
  
  # ENC plot
  enc_plot(enc_results, gc_results, file.path(output_dir, "enc_plot.pdf"))
  
  # PR2 bias plot
  pr2_bias_plot(codon_counts, file.path(output_dir, "pr2_plot.pdf"))
  
  # 6. Generate summary statistics
  message("Calculating summary statistics...")
  
  stats <- list(
    n_genes = nrow(codon_counts),
    enc_mean = mean(enc_results$ENC, na.rm = TRUE),
    enc_median = median(enc_results$ENC, na.rm = TRUE),
    enc_sd = sd(enc_results$ENC, na.rm = TRUE),
    gc_mean = mean(gc_results$GC, na.rm = TRUE),
    gc3s_mean = mean(gc_results$GC3s, na.rm = TRUE),
    gc12_vs_gc3_cor = cor(gc_results$GC12, gc_results$GC3, 
                          use = "complete.obs")
  )
  
  # 7. Save results
  message("Saving results...")
  fwrite(summary_table, file.path(output_dir, "cub_analysis_complete.csv"))
  fwrite(enc_results, file.path(output_dir, "enc_values.csv"))
  fwrite(gc_results, file.path(output_dir, "gc_content.csv"))
  
  # Save summary statistics
  stats_df <- data.frame(
    Metric = names(stats),
    Value = unlist(stats)
  )
  fwrite(stats_df, file.path(output_dir, "summary_statistics.csv"))
  
  # Print summary
  message("\n=== CUB Analysis Summary ===")
  message(sprintf("Number of genes analyzed: %d", stats$n_genes))
  message(sprintf("Mean ENC: %.2f (SD: %.2f)", stats$enc_mean, stats$enc_sd))
  message(sprintf("Median ENC: %.2f", stats$enc_median))
  message(sprintf("Mean GC content: %.2f%%", stats$gc_mean * 100))
  message(sprintf("Mean GC3s: %.2f%%", stats$gc3s_mean * 100))
  message(sprintf("GC12 vs GC3 correlation: %.3f", stats$gc12_vs_gc3_cor))
  message(sprintf("\nAll results saved to: %s", output_dir))
  
  return(list(
    summary_table = summary_table,
    enc_results = enc_results,
    gc_results = gc_results,
    rscu_results = rscu_results,
    statistics = stats
  ))
}


create_aa_specific_logos <- function(codon_counts, genetic_code, 
                                     output_dir = "./results/codon_logos")
{
  #' Create codon logos for all amino acids
  #' 
  #' @description Generates sequence logo-style plots for each amino acid
  #' showing nucleotide preferences at each codon position
  #' 
  #' @param codon_counts Data table with codon counts per gene
  #' @param genetic_code Named vector mapping codons to amino acids
  #' @param output_dir Directory for output files
  #' 
  #' @return Invisible
  #' ___________________________________________________________________________
  
  if(!dir.exists(output_dir))
  {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get unique amino acids (excluding STOP)
  amino_acids <- unique(genetic_code[genetic_code != "STOP"])
  
  message(sprintf("Creating codon logos for %d amino acids...", 
                  length(amino_acids)))
  
  for(aa in amino_acids)
  {
    tryCatch({
      output_file <- file.path(output_dir, paste0("codon_logo_", aa, ".pdf"))
      create_codon_logo(codon_counts, genetic_code, aa, output_file)
    }, error = function(e) {
      message(sprintf("Warning: Could not create logo for %s: %s", aa, e$message))
    })
  }
  
  message(sprintf("Codon logos saved to: %s", output_dir))
  
  return(invisible())
}
