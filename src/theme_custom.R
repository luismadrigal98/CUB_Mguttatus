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

save_figure_panel <- function(plot, name, width, height, dir = "./results/figure_panels")
{
  #' Save a copy of a manuscript panel under its figure label (e.g. "Fig7C_..."),
  #' so all panels for the Inkscape assembly sit in one folder. The pipeline's
  #' descriptive filenames stay as they are; only these copies carry figure
  #' numbers, which change with the manuscript and are cheap to rename here.
  #' ___________________________________________________________________________

  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  path <- file.path(dir, paste0(name, ".pdf"))
  ggplot2::ggsave(path, plot, width = width, height = height, device = dev)
  invisible(path)
}
