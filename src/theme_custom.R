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