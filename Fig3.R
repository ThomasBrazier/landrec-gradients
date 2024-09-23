# Loading env
source("Source/init_knitr.R")
source("Source/set_ggplot_theme.R")

gc()

fontsize = 20

# Figure 3
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



meanUTR_length = data.frame(species = list_species)
meanUTR_length$utr5 = NA
meanUTR_length$utr3 = NA
for (i in 1:nrow(meanUTR_length)) {
  meanUTR_length$utr5[i] = mean(data_all$width[which(data_all$species == meanUTR_length$species[i] & data_all$feature == "utr5")], na.rm = TRUE)
  meanUTR_length$utr3[i] = mean(data_all$width[which(data_all$species == meanUTR_length$species[i] & data_all$feature == "utr3")], na.rm = TRUE)
}



fontsize = 20
dotsize = 0.2
linesize = 1

# Gradient ATG
# Recombination gradient in exons (bp)
interval = 200
upstream = 2000 # Size of the upstream region
genic = 2000 # Max size of the genic region

df_gradient = readRDS(file = paste("Data/Recombination/Gradient/medianRho.rds", sep = ""))

df_gradient = df_gradient[which(df_gradient$set %in% list_dataset),]

atg_tss = data.frame(set = list_dataset)
atg_tss$dist = NA
for (i in 1:nrow(atg_tss)) {
  atg_tss$dist[i] = - mean(df_gradient$atg_tss_dist[which(df_gradient$set == atg_tss$set[i])], na.rm = TRUE)
}


for (i in 1:nrow(df_gradient)) {
  genus = strsplit(df_gradient$set[i], split = "_")[[1]][1]
  sp = strsplit(df_gradient$set[i], split = "_")[[1]][2]
  df_gradient$species[i] = paste(genus, sp, sep = " ")
}


atg_tss = data.frame(species = list_species)
atg_tss$dist = NA
for (i in 1:nrow(atg_tss)) {
  atg_tss$dist[i] = - mean(df_gradient$atg_tss_dist[which(df_gradient$species == atg_tss$species[i])], na.rm = TRUE)
}

df_gradient_subset = df_gradient[which(df_gradient$hotoverlap == "both" & df_gradient$level == "Genic" & df_gradient$position == "ATG"),]
df_gradient_subset$condition = factor(df_gradient_subset$condition, levels = c("observed", "control"), labels = c("Observed", "Control"))
df_gradient_subset = df_gradient_subset[which(df_gradient_subset$species %in% list_species),]

# Rescale control
for (s in list_species) {
  # All
  obs = df_gradient_subset$medianRho[which(df_gradient_subset$species == s & df_gradient_subset$condition == "Observed")]
  m.obs = mean(obs, na.rm = TRUE)
  # sd.obs = sd(obs, na.rm = TRUE)
  control = df_gradient_subset$medianRho[which(df_gradient_subset$species == s & df_gradient_subset$condition == "Control")]
  m.control = mean(control, na.rm = TRUE)
  # sd.control = sd(control, na.rm = TRUE)
  # sd.control = 1
  # sd.obs = 1
  df_gradient_subset$medianRho[which(df_gradient_subset$species == s & df_gradient_subset$condition == "Control")] = ((control - m.control)) + m.obs
}


p1 = ggplot(data = df_gradient_subset, aes(x = idx/1000, y = medianRho, group = condition, linetype = condition)) +
  geom_line(size = linesize) +
  xlim(-upstream/1000, genic/1000) +
  xlab("Distance to the ATG (kb)") +
  ylab("Median recombination rate (ρ/kb)\n  ") +
  labs(linetype = "Condition") +
  # geom_smooth(method = "loess", se = T) +
  # geom_vline(aes(xintercept = - mean(atg_tss_dist)), linetype = "dashed") +
  facet_wrap(. ~ species, scales = "free_y", ncol = 1) +
  geom_vline(xintercept = 0, size = 1) +
  geom_vline(data = atg_tss, aes(xintercept = dist/1000), linetype = "dashed", size = 1, color = "#6E7889") +
  # geom_vline(data = meanUTR_length, aes(xintercept = utr5), linetype = "dashed", colour = "darkgrey", size = 1) +
  ggtitle("5' end") +
  facetted_pos_scales(y = list(species == "Arabidopsis thaliana" ~ scale_y_continuous(limits = c(0.8, 1.2)),
                               species == "Glycine max" ~ scale_y_continuous(limits = c(0.39, 0.56)),
                               species == "Populus tremula" ~ scale_y_continuous(limits = c(2, 5)))) +
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
        legend.position='none')
p1


# Gradient TTS
df_gradient_subset = df_gradient[which(df_gradient$hotoverlap == "both" & df_gradient$level == "Genic" & df_gradient$position == "TTS"),]
df_gradient_subset$condition = factor(df_gradient_subset$condition, levels = c("observed", "control"), labels = c("Observed", "Control"))
df_gradient_subset = df_gradient_subset[which(df_gradient_subset$species %in% list_species),]

# Rescale control
for (s in list_species) {
  # All
  obs = df_gradient_subset$medianRho[which(df_gradient_subset$species == s & df_gradient_subset$condition == "Observed")]
  m.obs = mean(obs, na.rm = TRUE)
  control = df_gradient_subset$medianRho[which(df_gradient_subset$species == s & df_gradient_subset$condition == "Control")]
  m.control = mean(control, na.rm = TRUE)
  df_gradient_subset$medianRho[which(df_gradient_subset$species == s & df_gradient_subset$condition == "Control")] = ((control - m.control)) + m.obs
}

p2 = ggplot(data = df_gradient_subset, aes(x = idx/1000, y = medianRho, group = condition, linetype = condition)) +
  geom_line(size = linesize) +
  xlim(-upstream/1000, genic/1000) +
  xlab("Distance to the TTS (kb)") +
  ylab(" ") +
  labs(linetype = "Condition") +
  # geom_smooth(method = "loess", se = T) +
  # geom_vline(aes(xintercept = - mean(atg_tss_dist)), linetype = "dashed") +
  facet_wrap(. ~ species, scales = "free_y", ncol = 1, nrow = 3) +
  geom_vline(xintercept = 0, size = 1) +
  geom_vline(data = meanUTR_length, aes(xintercept = -utr3/1000), linetype = "dashed", colour = "#6E7889", size = 1) +
  ggtitle("3' end") +
  facetted_pos_scales(y = list(species == "Arabidopsis thaliana" ~ scale_y_continuous(limits = c(0.8, 1.2)),
                               species == "Glycine max" ~ scale_y_continuous(limits = c(0.39, 0.56)),
                               species == "Populus tremula" ~ scale_y_continuous(limits = c(2, 5)))) +
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
p2

# 
# # df_gradient2 = readRDS(file = paste("Data/Recombination/Gradient/meanRho_hotoverlap.rds", sep = ""))
# # 
# # df_gradient2$hotoverlap = factor(df_gradient2$hotoverlap, levels = c("hot", "cold"), labels = c("Hot", "Cold"))
# # 
# # df_gradient2 = df_gradient2[which(df_gradient2$set %in% list_dataset),]
# 
# # Relative Rho
# df_subset = df_gradient[which(df_gradient$condition == "observed" & df_gradient$position == "TTS"),]
# 
# for (i in list_dataset) {
#   # df_subset$meanRho[which(df_subset$set == i)] = scale(df_subset$meanRho[which(df_subset$set == i)])
#   for (j in c("Hot", "Cold")) {
#       df_subset$meanRho[which(df_subset$set == i & df_subset$hotoverlap == j)] = df_subset$meanRho[which(df_subset$set == i & df_subset$hotoverlap == j)] - mean(df_subset$meanRho[which(df_subset$set == i & df_subset$hotoverlap == j)], na.rm = TRUE)
#   }
# }
# 
# 
# df_subset$species = df_subset$set
# for (i in 1:nrow(df_subset)) {
#   genus = strsplit(df_subset$species[i], split = "_")[[1]][1]
#   sp = strsplit(df_subset$species[i], split = "_")[[1]][2]
#   df_subset$species[i] = paste(genus, sp, sep = " ")
# }
# 
# df_subset = df_subset[which(df_subset$species %in% list_species),]


# Gradient Hot vs Cold ATG
df_gradient_subset = df_gradient[which(df_gradient$hotoverlap %in% c("hot", "cold") & df_gradient$level == "Genic" & df_gradient$position == "ATG"),]
df_gradient_subset$condition = factor(df_gradient_subset$condition, levels = c("observed", "control"), labels = c("Observed", "Control"))
df_gradient_subset$hotoverlap = factor(df_gradient_subset$hotoverlap, levels = c("hot", "cold"), labels = c("Hot", "Cold"))

df_gradient_subset = df_gradient_subset[which(df_gradient_subset$species %in% list_species),]



# All in one figure
p3 = ggplot(data = dist2hotspot_filtered_2[which(dist2hotspot_filtered_2$species %in% list_species),], aes(x = idx/10^3, y = geneCount.start/n.hotspots)) +
  # geom_line(data = dist2hotspot[which(dist2hotspot$species %in% list_species),], aes(x = idx/10^3, y = geneCount.start/n.hotspots), colour = brewer.pal(3, "Set1")[1]) +
  # geom_line(data = dist2hotspot[which(dist2hotspot$species %in% list_species),], aes(y = geneCount.start.control/n.hotspots), colour = brewer.pal(3, "Set1")[1], alpha = 0.4) +
  # geom_line(data = dist2hotspot_filtered_2[which(dist2hotspot$species %in% list_species),], aes(x = idx/10^3, y = geneCount.start/n.hotspots), colour = brewer.pal(3, "Set1")[2]) +
  # geom_line(data = dist2hotspot_filtered_2[which(dist2hotspot$species %in% list_species),], aes(y = geneCount.start.control/n.hotspots), colour = brewer.pal(3, "Set1")[2], alpha = 0.4) +
  geom_line(colour = brewer.pal(3, "Set1")[2], size = linesize) +
  geom_line(data = dist2hotspot_filtered_2[which(dist2hotspot_filtered_2$species %in% list_species),], aes(y = geneCount.start.control/n.hotspots), colour = brewer.pal(3, "Set1")[2], alpha = 0.4, size = linesize) +
  # geom_vline(xintercept = 0, color = "Grey", size = 1) +
  xlab("Distance to hotspot centre (kb)") + ylab("Overlapping TSS density\n  ") +
  facet_wrap(~ species, nrow = 3, scale = "free_y") +
  ggtitle("5' end") +
  scale_x_continuous(breaks = c(-5, -2.5, 0, 2.5, 5), labels = c("-5", "-2.5", "0", "2.5", "5")) +
  facetted_pos_scales(y = list(species == "Arabidopsis thaliana" ~ scale_y_continuous(limits = c(0.22, 0.40)),
                               species == "Glycine max" ~ scale_y_continuous(limits = c(0.049, 0.06)),
                               species == "Populus tremula" ~ scale_y_continuous(limits = c(0.06, 0.20))
            )) +
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
        legend.position='none')
p3


df_hotspot = bind_rows(dist2hotspot[which(dist2hotspot$species %in% list_species),],
                       dist2hotspot_filtered_2[which(dist2hotspot_filtered_2$species %in% list_species),],
                       dist2hotspot_filtered_4[which(dist2hotspot_filtered_4$species %in% list_species),])

df1 = dist2hotspot_filtered_2[which(dist2hotspot_filtered_2$species %in% list_species),]
df1$TTSdensity = df1$geneCount.end/df1$n.hotspots
df1$Condition = "Observed"
df2 = dist2hotspot_filtered_2[which(dist2hotspot_filtered_2$species %in% list_species),]
df2$TTSdensity = df2$geneCount.end.control/df2$n.hotspots
df2$Condition = "Control"

df = rbind(df1, df2)

df$Condition = factor(df$Condition,
levels = c("Observed", "Control"))

p4 = ggplot(data = df, aes(x = idx/10^3, y = TTSdensity, alpha = Condition)) +
  geom_line(colour = brewer.pal(3, "Set1")[2], size = linesize) +
  scale_alpha_manual(values = c(1, 0.4)) +
  # geom_line(colour = brewer.pal(3, "Set1")[2], size = linesize) +
  # geom_line(data = dist2hotspot_filtered_2[which(dist2hotspot_filtered_2$species %in% list_species),], aes(y = geneCount.end.control/n.hotspots), colour = brewer.pal(3, "Set1")[2], alpha = 0.4, size = linesize) +
  # geom_line(data = df_hotspot, aes(x = idx/10^3, y = geneCount.end/n.hotspots, group = legend, colour = legend)) +
  # geom_line(aes(y = geneCount.end.control/n.hotspots), alpha = 0.4) +
  xlab("Distance to hotspot centre (kb)") + ylab("Overlapping TTS density\n  ") +
  # geom_vline(xintercept = 0, color = "Grey", size = 1) +
  # scale_color_manual(name='Filtering strategy',
  #                    breaks=c('Unfiltered', 'Soft', 'Hard'),
  #                    values=c('Unfiltered'=brewer.pal(3, "Set1")[1],
  #                             'Soft'=brewer.pal(3, "Set1")[2],
  #                             'Hard'=brewer.pal(3, "Set1")[3])) +
  facet_wrap(~ species, nrow = 3, scale = "free_y") +
  ggtitle("3' end") +
  scale_x_continuous(breaks = c(-5, -2.5, 0, 2.5, 5), labels = c("-5", "-2.5", "0", "2.5", "5")) +
  facetted_pos_scales(y = list(species == "Arabidopsis thaliana" ~ scale_y_continuous(limits = c(0.22, 0.40)),
                               species == "Glycine max" ~ scale_y_continuous(limits = c(0.049, 0.06)),
                               species == "Populus tremula" ~ scale_y_continuous(limits = c(0.06, 0.20))
            )) +
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


p = ggpubr::ggarrange(p1, p2, p3, p4, ncol = 4, widths = c(3,4,3,4), labels = "AUTO")
# p = ggpubr::ggarrange(p1, p2, ncol = 2, widths = c(3,4), labels = "AUTO")
p


# ggsave(file = paste("Figure/Paper/Fig3.tiff", sep = ""), plot = p, width = 23, height = 12, dpi = 300, compression = "lzw")
ggsave(file = paste("Figure/Paper/Fig3.jpeg", sep = ""), plot = p, width = 23, height = 10, dpi = 300)
ggsave(file = paste("Figure/Paper/Fig3.jpeg", sep = ""), plot = p, width = 23, height = 10, dpi = 600)

