create_aa_logo <- function(pspm_list, output_dir = "./results/codon_logos",
                           nrows = 6, width = 8, height = 8)
{
  #' Create codon logos for all amino acids
  #' 
  #' @description Generates sequence logo-style plots for each amino acid
  #' showing nucleotide preferences at each codon position
  #' 
  #' @param pspm_list List with PSPM per amino acid
  #' @param output_dir Directory for output files
  #' 
  #' @return Invisible
  #' ___________________________________________________________________________
  
  if(!dir.exists(output_dir))
  {
    dir.create(output_dir, recursive = TRUE)
  }
  
  plot_logos <- gridExtra::grid.arrange(grobs = lapply(X = names(pspm_list), 
                                                       FUN = function(aa)
    {
      pspm_aa <- pspm_list[[aa]]
                                                         
      ggseqlogo(pspm_aa, method = "prob") + 
      ggtitle(aa) +
      theme(plot.title = element_text(hjust = 0.5, size = 10))
    }), 
    nrow = nrows)
  
  ggsave(filename = file.path(output_dir, "aa_logos.pdf"), plot = plot_logos, 
         height = height, width = width)
  
  message(sprintf("Codon logos saved to: %s", output_dir))
  
  return(invisible())
}