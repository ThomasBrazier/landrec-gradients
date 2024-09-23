# Loading env
source("Source/init_knitr.R")
source("Source/set_ggplot_theme.R")


fontsize = 16
dotsize = 0.2
linesize = 1.5


ld = data.frame()
for (ds in list_dataset) {
  bin_file = paste0("Data/Recombination/LD/r2/", ds, ".ld_decay_bins")
  if (file.exists(bin_file)) {
    ld_bins = read_tsv(bin_file)
    
    bins = data.frame(start = seq(0, 200000, 1000),
                      end = (seq(0, 200000, 1000) + 999))
    bins$pos = round((bins$start + bins$end)/2, digits = 0) / 10^3
    bins$r2 = NA
    for (i in 1:nrow(bins)) {
      bins$r2[i] = mean(ld_bins$avg_R2[which(ld_bins$distance >= bins$start[i] & ld_bins$distance <= bins$end[i])])
    }
    bins$dataset = ds
    ld = rbind(ld, bins)
  }
}

ltype = (rep(c("solid", "dashed", "twodash"), 4))


p1 = ggplot(ld, aes(x = pos, y = r2, group = dataset, color = dataset, linetype = dataset)) +
  geom_line() +
  scale_linetype_manual(values = ltype) +
  xlab("Distance (kb)") + ylab(expression(italic(r)^2)) +
  xlim(0, 200) +
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
        legend.position='none')
p1



# Genic LD ----
ld_genic = data.frame()
for (ds in list_dataset) {
  bin_file = paste0("Data/Recombination/LD/r2/", ds, "_genic.ld_decay_bins")
  if (file.exists(bin_file)) {
    ld_bins = read_tsv(bin_file)
    
    bins = data.frame(start = seq(0, 200000, 1000),
                      end = (seq(0, 200000, 1000) + 999))
    bins$pos = round((bins$start + bins$end)/2, digits = 0) / 10^3
    bins$r2 = NA
    for (i in 1:nrow(bins)) {
      bins$r2[i] = mean(ld_bins$avg_R2[which(ld_bins$distance >= bins$start[i] & ld_bins$distance <= bins$end[i])])
    }
    bins$dataset = ds
    ld_genic = rbind(ld_genic, bins)
  }
}

for (i in 1:nrow(ld_genic)) {
  ld_genic$species[i] = dataset_metadata$species[which(dataset_metadata$dataset == ld_genic$dataset[i])]
}

ld_genic$species = as.factor(ld_genic$species)
ld_genic$species = factor(ld_genic$species,
                          levels = list_species_ordered,
                          labels = list_species_short)

p2 = ggplot(ld_genic, aes(x = pos, y = r2, group = species, color = species, linetype = species)) +
  geom_line() +
  scale_linetype_manual(values = ltype) +
  xlab("Distance (kb)") + ylab(expression(italic(r)^2)) +
  xlim(0, 200) +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        axis.text.x=element_text(),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize - 4, face = "italic"),
        legend.title=element_text(size=fontsize - 2),
        legend.position='right')
p2

p = ggarrange(p1, p2, ncol = 2, widths = c(4,5), labels = "AUTO")
p

ggsave(file = paste("Figure/Paper/FigS8.jpeg", sep = ""), plot = p, width = 21, height = 7, dpi = 300)
