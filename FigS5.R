# Loading env
source("Source/init_knitr.R")
source("Source/set_ggplot_theme.R")


dist2hotspot = readRDS("Data/Recombination/rho_hotspot_center.size100000000.intensity0.maxintensity100000000_All.rds")
# Soft filtering: Intensity > 2 and < 400
dist2hotspot_filtered_2 = readRDS("Data/Recombination/rho_hotspot_center.size10000.intensity0.maxintensity100000000_All.rds")
# Hard filtering: Intensity > 4 and < 200
dist2hotspot_filtered_4 = readRDS("Data/Recombination/rho_hotspot_center.size10000.intensity4.maxintensity200_All.rds")



for (i in 1:nrow(dist2hotspot)) {
  genus = strsplit(dist2hotspot$dataset[i], split = "_")[[1]][1]
  sp = strsplit(dist2hotspot$dataset[i], split = "_")[[1]][2]
  dist2hotspot$species[i] = paste(genus, sp, sep = " ")
}
dist2hotspot = dist2hotspot[which(dist2hotspot$dataset %in% list_dataset),]


for (i in 1:nrow(dist2hotspot_filtered_2)) {
  genus = strsplit(dist2hotspot_filtered_2$dataset[i], split = "_")[[1]][1]
  sp = strsplit(dist2hotspot_filtered_2$dataset[i], split = "_")[[1]][2]
  dist2hotspot_filtered_2$species[i] = paste(genus, sp, sep = " ")
}
dist2hotspot_filtered_2 = dist2hotspot_filtered_2[which(dist2hotspot_filtered_2$dataset %in% list_dataset),]


for (i in 1:nrow(dist2hotspot_filtered_4)) {
  genus = strsplit(dist2hotspot_filtered_4$dataset[i], split = "_")[[1]][1]
  sp = strsplit(dist2hotspot_filtered_4$dataset[i], split = "_")[[1]][2]
  dist2hotspot_filtered_4$species[i] = paste(genus, sp, sep = " ")
}
dist2hotspot_filtered_4 = dist2hotspot_filtered_4[which(dist2hotspot_filtered_4$dataset %in% list_dataset),]

# Scale by the min
for (i in 1:nrow(dist2hotspot_filtered_4)) {
  dist2hotspot_filtered_4$minrho[i] = min(dist2hotspot_filtered_4$weighted.meanRho[dist2hotspot_filtered_4$dataset == dist2hotspot_filtered_4$dataset[i]])
  dist2hotspot_filtered_4$maxrho[i] = max(dist2hotspot_filtered_4$weighted.meanRho[dist2hotspot_filtered_4$dataset == dist2hotspot_filtered_4$dataset[i]])
  dist2hotspot_filtered_4$sd.rho[i] = sd(dist2hotspot_filtered_4$weighted.meanRho[dist2hotspot_filtered_4$dataset == dist2hotspot_filtered_4$dataset[i]])
  
  dist2hotspot_filtered_4$minrho.control[i] = min(dist2hotspot_filtered_4$weighted.meanRho.control[dist2hotspot_filtered_4$dataset == dist2hotspot_filtered_4$dataset[i]])
  # dist2hotspot_filtered_4$mean.Rho.genome[i] = dataset_infos$mean.Rho.genome[dataset_infos$dataset == dist2hotspot_filtered_4$dataset[i]]

  dist2hotspot_filtered_4$mean.Rho.bg[i] = mean(dist2hotspot_filtered_4$weighted.meanRho.control[dist2hotspot_filtered_4$dataset == dist2hotspot_filtered_4$dataset[i]])
  dist2hotspot_filtered_4$bg.Rho[i] = min(dist2hotspot_filtered_4$weighted.meanRho[dist2hotspot_filtered_4$dataset == dist2hotspot_filtered_4$dataset[i]])
  
  dist2hotspot_filtered_4$max.Rho.hotspot[i] = max(dist2hotspot_filtered_4$weighted.meanRho[dist2hotspot_filtered_4$dataset == dist2hotspot_filtered_4$dataset[i]])
  
  dist2hotspot_filtered_4$SNPCount.genome[i] = mean(dist2hotspot_filtered_4$SNPCount.control[dist2hotspot_filtered_4$dataset == dist2hotspot_filtered_4$dataset[i]])
  dist2hotspot_filtered_4$geneCount.genome[i] = mean(dist2hotspot_filtered_4$geneCount.control[dist2hotspot_filtered_4$dataset == dist2hotspot_filtered_4$dataset[i]])
  dist2hotspot_filtered_4$geneCount.start.genome[i] = mean(dist2hotspot_filtered_4$geneCount.start.control[dist2hotspot_filtered_4$dataset == dist2hotspot_filtered_4$dataset[i]])
  dist2hotspot_filtered_4$geneCount.end.genome[i] = mean(dist2hotspot_filtered_4$geneCount.end.control[dist2hotspot_filtered_4$dataset == dist2hotspot_filtered_4$dataset[i]])
}



rho.species = factor(dist2hotspot_filtered_4$dataset)
rsp = levels(rho.species)
for (i in 1:length(rsp)) {
  map = read.table(gzfile(paste0("Data/Recombination/LD/ldhat/", rsp[i],".csv.gz")), header = TRUE)
  # rsp[i] = mean(map$Mean_rho, na.rm = TRUE)
  rsp[i] = weighted.mean(map$Mean_rho, (map$end - map$start), na.rm = TRUE)
}

dist2hotspot_filtered_4$rho.species = as.numeric(as.character(factor(rho.species,
                        levels = levels(rho.species),
                        labels = rsp)))
dist2hotspot_filtered_2$rho.species = as.numeric(as.character(factor(rho.species,
                        levels = levels(rho.species),
                        labels = rsp)))
dist2hotspot$rho.species = as.numeric(as.character(factor(rho.species,
                        levels = levels(rho.species),
                        labels = rsp)))


rho.species = factor(dist2hotspot_filtered_4$dataset)
rsp = levels(rho.species)
for (i in 1:length(rsp)) {
  map = read.table(gzfile(paste0("Data/Recombination/LD/ldhat/", rsp[i],".csv.gz")), header = TRUE)
  rsp[i] = sd(map$Mean_rho, na.rm = TRUE)
}

dist2hotspot_filtered_4$sd.rho.species = as.numeric(as.character(factor(rho.species,
                        levels = levels(rho.species),
                        labels = rsp)))
dist2hotspot_filtered_2$sd.rho.species = as.numeric(as.character(factor(rho.species,
                        levels = levels(rho.species),
                        labels = rsp)))


dist2hotspot_filtered_4$weighted.meanRho.scaled = dist2hotspot_filtered_4$weighted.meanRho / dist2hotspot_filtered_4$rho.species
dist2hotspot_filtered_4$weighted.meanRho.scaled.control = dist2hotspot_filtered_4$meanRho.control / dist2hotspot_filtered_4$rho.species


dist2hotspot_filtered_2$weighted.meanRho.scaled = dist2hotspot_filtered_2$weighted.meanRho / dist2hotspot_filtered_2$rho.species
dist2hotspot_filtered_2$weighted.meanRho.scaled.control = dist2hotspot_filtered_2$meanRho.control / dist2hotspot_filtered_2$rho.species


dist2hotspot$weighted.meanRho.scaled = dist2hotspot$weighted.meanRho / dist2hotspot$rho.species
dist2hotspot$weighted.meanRho.scaled.control = dist2hotspot$weighted.meanRho.control / dist2hotspot$rho.species

dist2hotspot_filtered_2$meanRho.normalized = (dist2hotspot_filtered_2$meanRho - dist2hotspot_filtered_2$rho.species) / dist2hotspot_filtered_2$sd.rho.species

dist2hotspot_filtered_4$meanRho.normalized = (dist2hotspot_filtered_4$meanRho - dist2hotspot_filtered_4$rho.species) / dist2hotspot_filtered_4$sd.rho.species



ld_maps_10 = readRDS("Data/Recombination/ldmap_windows_10kb.rds")

ld_maps_10$species_short = factor(ld_maps_10$species,
                            levels = list_species_ordered,
                            labels = list_species_short)

dist2hotspot$legend = "Unfiltered"
dist2hotspot_filtered_2$legend = "Soft"
dist2hotspot_filtered_4$legend = "Hard"

df = bind_rows(dist2hotspot,
           dist2hotspot_filtered_2,
           dist2hotspot_filtered_4)

df$species = factor(df$species,
                            levels = list_species_ordered,
                            labels = list_species_ordered)


df$species_short = factor(df$species,
                            levels = list_species_ordered,
                            labels = list_species_short)

                            
fontsize = 12
dotsize = 0.2
linesize = 1


dist2hotspot$species_label = factor(dist2hotspot$species,
                            levels = list_species_ordered,
                            labels = list_species_short)

dist2hotspot_filtered_2$species_label = factor(dist2hotspot_filtered_2$species,
                            levels = list_species_ordered,
                            labels = list_species_short)

dist2hotspot_filtered_4$species_label = factor(dist2hotspot_filtered_4$species,
                            levels = list_species_ordered,
                            labels = list_species_short)


# All in one figure
p3 = ggplot(data = dist2hotspot, aes(x = idx/10^3, y = geneCount.start/n.hotspots)) +
  geom_line(colour = brewer.pal(3, "Set1")[1]) +
  geom_line(data = dist2hotspot, aes(y = geneCount.start.control/n.hotspots), colour = brewer.pal(3, "Set1")[1], alpha = 0.4) +
  geom_line(data = dist2hotspot_filtered_2, aes(x = idx/10^3, y = geneCount.start/n.hotspots), colour = brewer.pal(3, "Set1")[2]) +
  geom_line(data = dist2hotspot_filtered_2, aes(y = geneCount.start.control/n.hotspots), colour = brewer.pal(3, "Set1")[2], alpha = 0.4) +
  geom_line(data = dist2hotspot_filtered_4, aes(x = idx/10^3, y = geneCount.start/n.hotspots), colour = brewer.pal(3, "Set1")[3]) +
  geom_line(data = dist2hotspot_filtered_4, aes(y = geneCount.start.control/n.hotspots), colour = brewer.pal(3, "Set1")[3], alpha = 0.4) +
  xlab("Distance to hotspot centre (kb)") + ylab("Overlapping TSS density") +
  facet_wrap(~ species_label, ncol = 4, scale = "free_y") +
  ggtitle("5' end") +
  scale_x_continuous(breaks = c(-5, 0, 5), labels = c("-5", "0", "5")) +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text = element_text(size=fontsize, colour="black"),
        axis.text.x = element_text(),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize, face = "italic"),
        legend.title=element_text(size=fontsize),
        legend.position='none')
p3


df_hotspot = bind_rows(dist2hotspot,
                       dist2hotspot_filtered_2,
                       dist2hotspot_filtered_4)

p4 = ggplot(data = df_hotspot, aes(x = idx/10^3, y = geneCount.end/n.hotspots, group = legend, colour = legend)) +
  geom_line() +
  geom_line(aes(y = geneCount.end.control/n.hotspots), alpha = 0.4) +
  xlab("Distance to hotspot centre (kb)") + ylab("Overlapping TTS density") +
  scale_color_manual(name='Filtering strategy',
                     breaks=c('Unfiltered', 'Soft', 'Hard'),
                     values=c('Unfiltered'=brewer.pal(3, "Set1")[1],
                              'Soft'=brewer.pal(3, "Set1")[2],
                              'Hard'=brewer.pal(3, "Set1")[3])) +
  facet_wrap(~ species_label, ncol = 4, scale = "free_y") +
  ggtitle("3' end") +
  scale_x_continuous(breaks = c(-5, 0, 5), labels = c("-5", "0", "5")) +
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
p4

p = ggpubr::ggarrange(p3, p4, ncol = 2, widths = c(3,4), labels = "AUTO")
p

ggsave(file = paste("Figure/Paper/FigS5.jpeg", sep = ""), plot = p, width = 23, height = 9, dpi = 300)
