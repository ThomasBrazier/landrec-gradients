# Loading env
source("Source/init_knitr.R")
source("Source/set_ggplot_theme.R")


# Figure 1

# I tested also the median Rho, which gave the same results
# Cleaner signal, yet strange GW median RHo values in C. sinensis and Fig 1C



# Mean rho count at +- 5kb of hotspot center
# No filtering
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

table(dist2hotspot_filtered_4$rho.species, dist2hotspot$dataset)
table(dist2hotspot_filtered_2$rho.species, dist2hotspot$dataset)
table(dist2hotspot$rho.species, dist2hotspot$dataset)


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
dist2hotspot$sd.rho.species = as.numeric(as.character(factor(rho.species,
                        levels = levels(rho.species),
                        labels = rsp)))

table(dist2hotspot_filtered_4$sd.rho.species, dist2hotspot$dataset)
table(dist2hotspot_filtered_2$sd.rho.species, dist2hotspot$dataset)
table(dist2hotspot$sd.rho.species, dist2hotspot$dataset)


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

fontsize = 14
dotsize = 0.2
linesize = 0.6


# CDF curves
p1 = ggplot(data = ld_maps_10, aes(x = cumulative.genome.proportion,
                              y = cumulative.recombination.proportion)) +
  geom_line(aes(group = seqnames)) +
  facet_wrap(~ species_short, nrow = 3) +
  xlab("Cumulative genome position") +
  ylab("Cumulative recombination") +
  geom_vline(xintercept = 0.2, linetype = "dashed") +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  theme(plot.title = element_text(color="black", size=fontsize,hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        axis.text.x=element_text(),
        strip.text.x = element_text(size = fontsize-2),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize, face = "italic"),
        legend.title=element_text(size=fontsize),
        legend.position='none')

p1


rho.species = factor(df$dataset)
rsp = levels(rho.species)
for (i in 1:length(rsp)) {
  map = read.table(gzfile(paste0("Data/Recombination/LD/ldhat/", rsp[i],".csv.gz")), header = TRUE)
  rsp[i] = median(map$Mean_rho, na.rm = TRUE)
}

df$median.rho.species = factor(rho.species,
                        levels = levels(rho.species),
                        labels = rsp)

df$median.rho.species = as.numeric(as.character(df$median.rho.species))

df$speciesRho[which(df$species == "Glycine max")] = mean(df$speciesRho[which(df$species == "Glycine max")], na.rm = TRUE) 

# Rho around hotspot centre - different filtering strategies
p2 = ggplot(data = df, aes(x = idx/10^3, y = meanRho, group = legend, colour = legend)) +
  geom_line() +
  geom_line(data = df, aes(y = meanRho.control, group = legend, colour = legend), alpha = 0.4) +
  geom_hline(aes(yintercept = rho.species, group = NA, colour = NA), colour = "black", linetype = "dashed") +
  xlab("Distance to hotspot centre (kb)") + ylab("Median ρ/kb") +
  scale_x_continuous(limits = c(-5, 5), breaks = c(-5, 0, 5)) +
  facet_wrap(~ species_short, ncol = 4, scale = "free") +
  theme(legend.position = "none") +
  scale_color_manual(name='Filtering strategy  ',
                     breaks=c('Unfiltered', 'Soft', 'Hard'),
                     values=c('Unfiltered'=brewer.pal(3, "Set1")[1],
                              'Soft'=brewer.pal(3, "Set1")[2],
                              'Hard'=brewer.pal(3, "Set1")[3])) +
  theme(plot.title = element_text(color="black", size=fontsize,hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        axis.text.x=element_text(),
        strip.text.x = element_text(size = fontsize-2),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize, face = "italic"),
        legend.title=element_text(size=fontsize),
        legend.position='bottom')
p2


dist2hotspot_filtered_2$species = factor(dist2hotspot_filtered_2$species,
                                    levels = list_species_ordered,
                                    labels = list_species_short)

# Order colors by peak size
colpal = hcl.colors(12, palette = "Dark3")
colorder = aggregate(meanRho.normalized ~ species, dist2hotspot_filtered_2, max)
ord = colorder$species[order(colorder$meanRho.normalized, decreasing = TRUE)]
colorder = arrange(colorder, desc(meanRho.normalized))
colorder$color = colpal

# Human in black
colorder$color[which(colorder$species == "H. sapiens")] = "#000000"

dist2hotspot_filtered_2_p3 = dist2hotspot_filtered_2 
dist2hotspot_filtered_2_p3$species = factor(dist2hotspot_filtered_2_p3$species,
                                    levels = list_species_short[ord],
                                    labels = list_species_short[ord])

# 
# colorder$species = factor(colorder$species,
#                           levels = list_species_short[ord],
#                           labels = list_species_short[ord])

ltype = (rep(c("solid", "dashed", "twodash"), 4))

# Rho/kb around hotspot centre
# Fixed scale to compare intensities
p3 = ggplot(data = dist2hotspot_filtered_2_p3, aes(x = idx/10^3, y = meanRho.normalized, group = species, color = species, linetype = species)) +
  geom_line() +
  scale_color_manual(values = colorder$color) +
  scale_linetype_manual(values = ltype) +
  # geom_line(data = dist2hotspot_filtered_4, aes(y = meanRho.scaled.control), alpha = 0.4, linetype = "dashed") +
  # geom_hline(yintercept = 4, linetype = "dashed") +
  # geom_hline(aes(yintercept = mean.Rho.genome, colour = species)) +
  xlab("Distance to hotspot centre (kb)") + ylab("Standardized mean recombination rate") +
  scale_x_continuous(limits = c(-5, 5), breaks = c(-5, -2.5, 0, 2.5, 5), labels = c("-5", "-2.5", "0", "2.5", "5")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        axis.text.x=element_text(),
        strip.text.x = element_text(size = fontsize),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(1.5,"line"),
        legend.key.width = unit(1.5,"line"),
        legend.text=element_text(size=fontsize-2, face = "italic"),
        legend.title=element_text(size=fontsize-2, colour = "white"),
        legend.position='right')
p3

# lab = levels(ldhot_trimmed_hard$dataset)
# for (i in 1:length(lab)) {
#   lab[i] = dataset_metadata$species[which(dataset_metadata$dataset == lab[i])]
# }
# 
# ldhot_trimmed_hard$species = factor(ldhot_trimmed_hard$dataset,
#                                     levels = levels(ldhot_trimmed_hard$dataset),
#                                     labels = lab)
# 
# ldhot_trimmed_hard$species = factor(ldhot_trimmed_hard$species,
#                                     levels = list_species_ordered,
#                                     labels = list_species_ordered)

# lab = levels(ldhot$dataset)
# for (i in 1:length(lab)) {
#   lab[i] = dataset_metadata$species[which(dataset_metadata$dataset == lab[i])]
# }
# 
# ldhot$species = factor(ldhot$dataset,
#                                     levels = levels(ldhot_trimmed_hard$dataset),
#                                     labels = lab)
# 
# ldhot$species = factor(ldhot$species,
#                                     levels = list_species_ordered,
#                                     labels = list_species_short)

# 
# p4a = ggplot(ldhot_trimmed_hard, aes(x = intensity, group = dataset, colour = dataset)) +
#   geom_density() +
#   xlim(0, 100) +
#   xlab("Hotspot intensity") +
#   theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
#         axis.title.x = element_text(color="black", size=fontsize),
#         axis.title.y = element_text(color="black", size=fontsize),
#         axis.text=element_text(size=fontsize, colour="black"),
#         axis.text.x=element_text(),
#         strip.text.x = element_text(size = fontsize-2),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.height = unit(2,"line"),
#         legend.key.width = unit(3,"line"),
#         legend.text=element_text(size=fontsize, face = "italic"),
#         legend.title=element_text(size=fontsize),
#         legend.position='none')
# 
# 
# p4b = ggplot(ldhot, aes(x = species, y = intensity, group = species, colour = species)) +
#   geom_violin(aes(fill = species)) +
#   geom_boxplot(alpha = 0.3) +
#   ylim(0, 200) +
#   ylab("Hotspot intensity") + xlab("") +
#   theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
#         axis.title.x = element_text(color="black", size=fontsize),
#         axis.title.y = element_text(color="black", size=fontsize),
#         axis.text=element_text(size=fontsize, colour="black"),
#         axis.text.x=element_text(angle = 90, face = "italic"),
#         strip.text.x = element_text(size = fontsize-2),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.height = unit(2,"line"),
#         legend.key.width = unit(3,"line"),
#         legend.text=element_text(size=fontsize, face = "italic"),
#         legend.title=element_text(size=fontsize),
#         legend.position='none')

# p4 = ggarrange(p4b, p4a, nrow = 2, labels = "")
# 
# p = ggarrange(p1, p2, p4b, p3, ncol = 2, nrow = 2, heights = c(4,4), widths = c(2,3), labels = "AUTO")

# p = ggpubr::ggarrange(p1, p2, p3, nrow = 3, heights = c(2,3,3), labels = "AUTO")

pa = ggpubr::ggarrange(p1, p3, nrow = 2, heights = c(3,3), labels = c("B", "C"))
p = ggpubr::ggarrange(p2, pa, ncol = 2, widths = c(3,2), labels = c("A", ""))

# p = ggarrange(pa, p3, nrow = 2, heights = c(4,4), labels = c("", "C")) 
p

ggsave(file = paste("Figure/Paper/Fig1.jpeg", sep = ""), plot = p, width = 15, height = 10, dpi = 300)

ggsave(file = paste("Figure/Paper/Fig1.jpeg", sep = ""), plot = p, width = 15, height = 10, dpi = 600)

