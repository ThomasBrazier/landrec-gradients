# Loading env
source("Source/init_knitr.R")
source("Source/set_ggplot_theme.R")


df_subset = data_all[which(data_all$feature == "CDS" & data_all$nb_exons <= max.exons),]

df_subset$hotoverlap = ifelse(df_subset$hotspot_overlap_intensity2 > 0, TRUE, FALSE)
# df_subset$hotoverlap = ifelse(df_subset$hotspot_overlap_raw > 0, TRUE, FALSE)

hotspot_overlap = aggregate(hotoverlap ~ species + rank + nb_exons, df_subset, sum)

hotspot_overlap$nb_exons = factor(hotspot_overlap$nb_exons)

hotspot_n = aggregate(hotoverlap ~ species + rank + nb_exons, df_subset, length)
hotspot_overlap$n_sample = hotspot_n$hotoverlap


fontsize = 18
dotsize = 0.2
linesize = 1.5

ggplot(hotspot_overlap, aes(x = rank, y = hotoverlap/n_sample, group = nb_exons, colour = nb_exons)) +
  geom_point() +
  geom_line() +
  scale_colour_viridis_d(option = "B") +
  xlab("CDS part rank") + ylab("Proportion of CDS overlapped by a hotspot") +
  labs(colour = "# exons") +
  facet_wrap(~ species, ncol = 4, scales = "free") +
  xlim(1, 14) +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        axis.text.x=element_text(),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize, face = "italic"),
        legend.title=element_text(size=fontsize),
        legend.position='right')

ggsave(file = paste("Figure/Paper/FigS9.jpeg", sep = ""), width = 18, height = 13, dpi = 300)
