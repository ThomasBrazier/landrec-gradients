# Loading env
source("Source/init_knitr.R")
source("Source/set_ggplot_theme.R")


fontsize = 16
dotsize = 0.2
linesize = 1.5


gene_range = readRDS("Data/Recombination/Gradient/gradients_deciles.rds")

gene_range$species = factor(gene_range$dataset,
levels = c("Arabidopsis_thaliana_1001genomes", "Glycine_max_Yang2021", "Populus_tremula_Liu2022",
"Homo_sapiens_Sudmant2015", "Camellia_sinensis_Zhang2021", "Citrullus_lanatus_Guo2019", "Malus_sieversii_Sun2020",
"Oryza_sativa_Wang2018", "Phaseolus_vulgaris_Wu2020", "Sorghum_bicolor_Lozano2021",
"Spinacia_oleracea_Cai2021", "Triticum_aestivum_Zhou2020"),
labels = list_species_ordered)


rhogradient = aggregate(mean.rho ~ index + nb_exons + species, gene_range, median)
rhogradient$index = as.factor(rhogradient$index)
rhogradient$nb_exons = as.factor(rhogradient$nb_exons)
ggplot(data = rhogradient, aes(x = index, y = mean.rho, group = nb_exons, colour = nb_exons)) +
  geom_point() +
  geom_line() +
  xlab("Decile of distance along CDS (10 bins)") +
  facet_wrap(~ species, scales = "free") +
  scale_color_viridis_d() +
  ylab("Median recombination rate (ρ/kb)") +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        axis.text.x=element_text(),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize - 2, face = "italic"),
        legend.title=element_text(size=fontsize - 2),
        legend.position='right')

ggsave(file = paste("Figure/Paper/FigS11.jpeg", sep = ""), width = 17, height = 11, dpi = 300)
