# Loading env
source("Source/init_knitr.R")
source("Source/set_ggplot_theme.R")



meanUTR_length = data.frame(species = list_species_ordered)
meanUTR_length$utr5 = NA
meanUTR_length$utr3 = NA
for (i in 1:nrow(meanUTR_length)) {
  meanUTR_length$utr5[i] = mean(data_all$width[which(data_all$species == meanUTR_length$species[i] & data_all$feature == "utr5")], na.rm = TRUE)
  meanUTR_length$utr3[i] = mean(data_all$width[which(data_all$species == meanUTR_length$species[i] & data_all$feature == "utr3")], na.rm = TRUE)
}



fontsize = 18
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


atg_tss = data.frame(species = list_species_ordered)
atg_tss$dist = NA
for (i in 1:nrow(atg_tss)) {
  atg_tss$dist[i] = - mean(df_gradient$atg_tss_dist[which(df_gradient$species == atg_tss$species[i])], na.rm = TRUE)
}

df_gradient_subset = df_gradient[which(df_gradient$hotoverlap == "both" & df_gradient$level == "Genic" & df_gradient$position == "ATG"),]
df_gradient_subset$condition = factor(df_gradient_subset$condition, levels = c("observed", "control"), labels = c("Observed", "Control"))
# df_gradient_subset = df_gradient_subset[which(df_gradient_subset$species %in% list_species),]

# Rescale control
for (s in list_species_ordered) {
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
  facet_wrap(. ~ species, scales = "free_y", ncol = 4) +
  geom_vline(xintercept = 0, size = 1) +
  geom_vline(data = atg_tss, aes(xintercept = dist/1000), linetype = "dashed", size = 1, color = "#6E7889") +
  # geom_vline(data = meanUTR_length, aes(xintercept = utr5), linetype = "dashed", colour = "darkgrey", size = 1) +
  ggtitle("5' end") +
  # facetted_pos_scales(y = list(species == "Arabidopsis thaliana" ~ scale_y_continuous(limits = c(0.8, 1.2)),
  #                              species == "Glycine max" ~ scale_y_continuous(limits = c(0.39, 0.56)),
  #                              species == "Populus tremula" ~ scale_y_continuous(limits = c(2, 5)))) +
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
# df_gradient_subset = df_gradient_subset[which(df_gradient_subset$species %in% list_species),]

# Rescale control
for (s in list_species_ordered) {
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
  facet_wrap(. ~ species, scales = "free_y", ncol = 4) +
  geom_vline(xintercept = 0, size = 1) +
  geom_vline(data = meanUTR_length, aes(xintercept = -utr3/1000), linetype = "dashed", colour = "#6E7889", size = 1) +
  ggtitle("3' end") +
  # facetted_pos_scales(y = list(species == "Arabidopsis thaliana" ~ scale_y_continuous(limits = c(0.8, 1.2)),
  #                              species == "Glycine max" ~ scale_y_continuous(limits = c(0.39, 0.56)),
  #                              species == "Populus tremula" ~ scale_y_continuous(limits = c(2, 5)))) +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        axis.text.x=element_text(),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize-4, face = "italic"),
        legend.title=element_text(size=fontsize-4),
        legend.position='right')
p2


p = ggpubr::ggarrange(p1, p2, ncol = 2, labels = "AUTO")
p

ggsave(file = paste("Figure/Paper/FigS4.jpeg", sep = ""), plot = p, width = 30, height = 10, dpi = 300)
