cub_summary <- function(codon_counts, genetic_code, output_dir = "./results",
                        aa_group = NULL)
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
  #' @param aa_group Grouping of amino acids based on chemistry. Expected
  #' format is a data.frame with two columns, AA and class.
  #' 
  #' @return List with all analysis results
  #' ___________________________________________________________________________
  
  require(data.table)
  
  # Create output directory if needed
  if(!dir.exists(output_dir))
  {
    dir.create(output_dir, recursive = TRUE)
  }
  
  message("Starting comprehensive CUB analysis...")
  
  # 1. Calculate RSCU ----
  message("Calculating RSCU...")
  rscu_results <- calculate_rscu(codon_counts, genetic_code)
  
  # 2. Calculate ENC ----
  message("Calculating ENC...")
  enc_results <- calculate_enc(codon_counts, genetic_code)
  
  # 3. Calculate the RF ----
  message("Calculating the RF...")
  rf_results <- calculate_rf(codon_counts, genetic_code)
  
  # 4. Calculate the PSPM ----
  message("Calculating the PSPM...")
  pspm_overall <- calculate_overall_PSPM(rf_results, genetic_code)
   
  # 5. Calculate GC content ----
  message("Calculating GC content metrics...")
  gc_results <- calculate_gc_content(codon_counts)
  
  # 6. Merge all results ----
  message("Merging results...")
  summary_table <- merge(codon_counts, enc_results, by = "Gene_name")
  summary_table <- merge(summary_table, gc_results, by = "Gene_name")
  
  # 7. Generate visualizations ----
  message("Creating visualizations...")
  
  # Codon usage heatmap
  visualize_codon_usage(codon_counts, genetic_code, 
                       file.path(output_dir, "codon_usage_heatmap.pdf"),
                       type = "heatmap")
  
  # Codon usage barplot
  visualize_codon_usage(codon_counts, genetic_code, 
                       file.path(output_dir, "codon_usage_barplot.pdf"),
                       type = "barplot",
                       aa_grouping = aa_group)
  
  # Neutrality plot
  neutrality_plot(gc_results, file.path(output_dir, "neutrality_plot.pdf"))
  
  # ENC plot
  enc_plot(enc_results, gc_results, file.path(output_dir, "enc_plot.pdf"))
  
  # PR2 bias plot
  pr2_bias_plot(codon_counts, file.path(output_dir, "pr2_plot.pdf"))
  
  # 8. Perform goodness of fit tests (G test) ----
  
  message("Performing G tests...")
  
  # Test 1: By gene (default)
  message("  - Testing individual genes...")
  g_test_by_gene <- CUB_g_test(codon_counts, genetic_code, 
                                mode = "by_gene", correct.p_values = TRUE)
  
  # Test 2: Genome-wide
  message("  - Testing genome-wide codon usage...")
  g_test_genome <- CUB_g_test(codon_counts, genetic_code, 
                               mode = "by_genome")
  
  # Test 3: By amino acid
  message("  - Testing codon usage per amino acid...")
  g_test_by_aa <- CUB_g_test(codon_counts, genetic_code, 
                              mode = "by_aminoacid", correct.p_values = TRUE)
  
  # Test 4: Heterogeneity across genes
  message("  - Testing codon usage heterogeneity across genes...")
  g_test_heterogeneity <- CUB_g_test(codon_counts, genetic_code, 
                                      mode = "heterogeneity_per_aa", 
                                      correct.p_values = TRUE)
  
  # Save all G-test results
  fwrite(g_test_by_gene, file.path(output_dir, "g_test_by_gene.csv"))
  fwrite(g_test_genome, file.path(output_dir, "g_test_genome_wide.csv"))
  fwrite(g_test_by_aa, file.path(output_dir, "g_test_by_aminoacid.csv"))
  fwrite(g_test_heterogeneity, file.path(output_dir, "g_test_heterogeneity.csv"))
  
  # 9. Generate summary statistics ----
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
  
  # 10. Save results ----
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
    rf_results = rf_results,
    pspm_results = pspm_overall,
    gc_results = gc_results,
    rscu_results = rscu_results,
    g_test_by_gene = g_test_by_gene,
    g_test_genome = g_test_genome,
    g_test_by_aa = g_test_by_aa,
    g_test_heterogeneity = g_test_heterogeneity,
    statistics = stats
  ))
}