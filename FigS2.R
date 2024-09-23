# Loading env
source("Source/init_knitr.R")
source("Source/set_ggplot_theme.R")

fontsize = 14

# Figure S2
ldmaps = read.table(gzfile("Data/Recombination/LD/ldhat/LD_maps.csv.gz"), header = T, sep = "\t")


sp = unique(ldmaps$set)

for (i in 1:length(sp)) {
  sp[i] = dataset_metadata$species[which(dataset_metadata$dataset == sp[i])]
}

ldmaps$set = as.factor(ldmaps$set)

ldmaps$species = factor(ldmaps$set,
                        levels = levels(ldmaps$set),
                        labels = sp)


ldmaps$species = factor(ldmaps$species,
                        levels = list_species_ordered)


fontsize = 18
dotsize = 0.2
linesize = 0.6

p1 = ggplot(data = ldmaps, aes(x = species, y = Mean_rho, fill = species)) +
  geom_boxplot(width = 0.4) + xlab("Species") + ylab("ρ/kb (log scale)") +
  geom_violin(alpha = 0.4) +
  coord_trans(y = "log10") +
  scale_y_continuous(breaks = c(1, 10, 50, 150)) +
  scale_fill_brewer(palette="Paired") +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text = element_text(size=fontsize, colour="black"),
        axis.text.x = element_text(angle = 90, face = "italic"),
        strip.text.x = element_text(size = fontsize),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize, face = "italic"),
        legend.title=element_text(size=fontsize),
        legend.position='none')
p1


# p2 = ggplot(data = ldmaps, aes(x = Mean_rho)) +
#   geom_density() + ylab("Density") + xlab("ρ/kb") +
#   # coord_trans(x = "log") +
#   facet_wrap(~ species, ncol = 3, scale = "free") +
#   theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
#         axis.title.x = element_text(color="black", size=fontsize),
#         axis.title.y = element_text(color="black", size=fontsize),
#         axis.text = element_text(size=fontsize, colour="black"),
#         axis.text.x = element_text(angle = 90, face = "italic"),
#         strip.text.x = element_text(size = fontsize),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.height = unit(2,"line"),
#         legend.key.width = unit(3,"line"),
#         legend.text=element_text(size=fontsize, face = "italic"),
#         legend.title=element_text(size=fontsize),
#         legend.position='none')
# p2


# p3 = ggplot(data = ldmaps, aes(x = Mean_rho)) +
#   geom_density() + ylab("Density") + xlab("ρ/kb") +
#   xlim(0, 25) + 
#   # coord_trans(x = "log") +
#   facet_wrap(~ species, ncol = 3, scale = "free") +
#   theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
#         axis.title.x = element_text(color="black", size=fontsize),
#         axis.title.y = element_text(color="black", size=fontsize),
#         axis.text = element_text(size=fontsize, colour="black"),
#         axis.text.x = element_text(angle = 90, face = "italic"),
#         strip.text.x = element_text(size = fontsize),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.height = unit(2,"line"),
#         legend.key.width = unit(3,"line"),
#         legend.text=element_text(size=fontsize, face = "italic"),
#         legend.title=element_text(size=fontsize),
#         legend.position='none')
# p3

# p = ggarrange(p1, p2, p3, nrow = 3, labels = "AUTO")
# p

ggsave(file = paste("Figure/Paper/FigS2.jpeg", sep = ""), plot = p1, width = 10, height = 8, dpi = 300)
