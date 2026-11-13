theme_custom <- function() 
{
  #' Custom theme for ggplot plots
  #' ___________________________________________________________________________
  
  theme_bw(base_size = 14) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.box = "horizontal",
      axis.text = element_text(color = "black")
    )
}
nucleotide_palette <- function(alphabet = c("DNA", "RNA"))
{
  #' One colour per third-position base, shared by every figure that colours
  #' codons by their ending (codon share by expression; preferred codons across
  #' plant species). Referee 2 asked for a single scheme: the two figures used
  #' to code A, C, G and T/U with different colours.
  #'
  #' Hues follow the usual sequence-logo convention (A green, C blue, G orange,
  #' T/U red), taken from the dark half of ColorBrewer "Paired" so lines stay
  #' legible; tile fills use the same colours with transparency.
  #' ___________________________________________________________________________

  alphabet <- match.arg(alphabet)
  pal <- c(A = "#33A02C", C = "#1F78B4", G = "#FF7F00", T = "#E31A1C")
  if (alphabet == "RNA") names(pal)[4] <- "U"
  pal
}
