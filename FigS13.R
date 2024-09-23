# Loading env
source("Source/init_knitr.R")
source("Source/set_ggplot_theme.R")

fontsize = 18
dotsize = 0.2
linesize = 1


correlation_rho_rank = data.frame(species = unique(data_all$species))

correlation_rho_rank$obs.spearman = NA
correlation_rho_rank$obs.pvalue = NA
correlation_rho_rank$random.spearman = NA
correlation_rho_rank$random.pvalue = NA
correlation_rho_rank$simulgradient.spearman = NA
correlation_rho_rank$simulgradient.pvalue = NA
correlation_rho_rank$startafter.spearman = NA
correlation_rho_rank$startafter.pvalue = NA

for (i in 1:nrow(correlation_rho_rank)) {
  df = data_all[which(data_all$species == correlation_rho_rank$species[i] & data_all$feature == "CDS" & data_all$rank <= max.exons ),]
  # Empirical - Correlation mean Rho ~ rank
  cor = cor.test(df$weighted.mean.rho, df$rank, method = "spearman")
  correlation_rho_rank$obs.spearman[i] = cor$estimate
  correlation_rho_rank$obs.pvalue[i] = cor$p.value
  
  # Simulation. Random expectation. Resampling ranks among exons/introns
  df$rank_resample = sample(df$rank, replace = TRUE)
  cor = cor.test(df$weighted.mean.rho, df$rank_resample, method = "spearman")
  correlation_rho_rank$random.spearman[i] = cor$estimate
  correlation_rho_rank$random.pvalue[i] = cor$p.value
  
  # Simulation. Power to detect a true gradient if there is one. Re-ordering mean Rho per decreasing order (within genes).
  # WARNING Do not work with trimmed data.
  if (sum(!is.na(df$simulated.gradient)) > 10) {
    cor = cor.test(df$simulated.gradient, df$rank, method = "spearman")
    correlation_rho_rank$simulgradient.spearman[i] = cor$estimate
    correlation_rho_rank$simulgradient.pvalue[i] = cor$p.value
  }
  
  # Simulation. Robust to reject an artefactual gradient. Eliminate every Rho estimate not beginning in the actual window.
  cor = cor.test(df$weighted.mean.rho.startafter, df$rank, method = "spearman")
  correlation_rho_rank$startafter.spearman[i] = cor$estimate
  correlation_rho_rank$startafter.pvalue[i] = cor$p.value
  
}

rho_rank = aggregate(weighted.mean.rho.startafter ~ rank + species, data = data_all[which(data_all$feature == "CDS" & data_all$rank < max.exons & data_all$rank > 0),], mean)

rho_rank.control = aggregate(weighted.mean.rho.control ~ rank + species, data = data_all[which(data_all$feature == "CDS" & data_all$rank < max.exons & data_all$rank > 0),], mean)

# Rescale the control data
# https://stats.stackexchange.com/questions/46429/transform-data-to-desired-mean-and-standard-deviation
for (s in unique(rho_rank$species)) {
  obs = rho_rank$weighted.mean.rho.startafter[which(rho_rank$species == s)]
  m.obs = mean(obs, na.rm = TRUE)
  sd.obs = sd(obs, na.rm = TRUE)
  control = rho_rank.control$weighted.mean.rho.control[which(rho_rank.control$species == s)]
  m.control = mean(control, na.rm = TRUE)
  sd.control = sd(control, na.rm = TRUE)
  rho_rank.control$weighted.mean.rho.control[which(rho_rank.control$species == s)] = ((control - m.control)) + m.obs
}

p1 = ggplot(data = rho_rank[which(rho_rank$species %in% list_species),], aes(x = rank, y = weighted.mean.rho.startafter)) +
  geom_line() +
  geom_line(data = rho_rank.control[which(rho_rank.control$species %in% list_species),], aes(x = rank, y = weighted.mean.rho.control), colour = "Black", linetype = "dashed") +
  xlab("Rank") + ylab("Recombination rate (ρ/kb)") +
  facet_wrap(. ~ species, scales = "free", nrow = 3) +
  facetted_pos_scales(y = list(species == "Arabidopsis thaliana" ~ scale_y_continuous(breaks = c(1.4, 1.8)),
                               species == "Glycine max" ~ scale_y_continuous(breaks = c(1, 1.5, 2)),
                               species == "Populus tremula" ~ scale_y_continuous(breaks = c(2.5, 3, 3.5)))) +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        axis.text.x=element_text(),
        strip.text.x = element_text(size = 14, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize, face = "italic"),
        legend.title=element_text(size=fontsize),
        legend.position='none')
p1

df = correlation_rho_rank[,c(1:2, 4, 6, 8)]

df = df %>% gather(type, correlation, -species)
df$correlation = round(df$correlation, digits = 2)

df$type = as.factor(df$type)
df$type = factor(df$type, levels = c("random.spearman", "obs.spearman", "startafter.spearman", "simulgradient.spearman"),
                 labels = c("Random expectation", "Observed", "SNP included", "Max gradient (simulated)"))


p2 = ggplot(data = df, aes(x = species, y = type, fill = correlation)) +
  geom_tile() +
  # geom_text(aes(label = correlation), colour = "darkorange3", face = "bold", size = 4) +
  geom_text(aes(label = correlation), colour = "brown3", fontface = "bold", size = 4) +
  scale_fill_viridis_c() +
  xlab("") + ylab("") +
  coord_fixed() +
  labs(fill = "Correlation    ") +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        axis.text.x=element_text(angle = 90, face = "italic"),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(1.5,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=fontsize-6),
        legend.title=element_text(size=fontsize-6),
        legend.position='bottom')
p2

p = ggpubr::ggarrange(p1, p2, ncol = 2, widths = c(1,4), labels = "AUTO")

p

ggsave(file = paste("Figure/Paper/FigS13.jpeg", sep = ""), plot = p, width = 15, height = 6, dpi = 300, bg = "White")
