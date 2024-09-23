# Loading env
source("Source/init_knitr.R")
source("Source/set_ggplot_theme.R")


# Figure S1
source("Jobs/FigS1.R")

fontsize = 18
dotsize = 0.2
linesize = 0.6

p = ggdraw() +
  draw_text("Supplementary PDF") +
  bgcolor("white")


ggsave(file = paste("Figure/Paper/FigS1.jpeg", sep = ""), plot = p, width = 10, height = 5, dpi = 300)