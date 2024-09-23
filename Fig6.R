# Loading env
source("Source/init_knitr.R")
source("Source/set_ggplot_theme.R")


fontsize = 18
dotsize = 0.2
linesize = 1.5


exons = data_all[which(data_all$feature == "CDS" & data_all$nb_exons <= max.exons),]
df_gradient_rank = aggregate(weighted.mean.rho ~ rank + nb_exons + species, exons, median)
df_gradient_rank$nb_exons = as.factor(df_gradient_rank$nb_exons)

# 
# df_gradient_rank = readRDS(file = paste("Data/Recombination/Gradient/meanRho_rank_nbexons.rds", sep = ""))
# 
# for (i in 1:nrow(df_gradient_rank)) {
#   genus = strsplit(df_gradient_rank$set[i], split = "_")[[1]][1]
#   sp = strsplit(df_gradient_rank$set[i], split = "_")[[1]][2]
#   df_gradient_rank$species[i] = paste(genus, sp, sep = " ")
# }
# alph = as.character(df_gradient_rank$species[-which(df_gradient_rank$nb_exons == "all")])
# alph = ifelse(as.character(alph) %in% list_species,
#                                0.3,
#                                1)

# df_gradient_rank$species = factor(df_gradient_rank$species,
#                                   levels = list_species_ordered,
#                                   labels = list_species_ordered)

df_gradient_rank$reference = ifelse(as.character(df_gradient_rank$species) %in% list_species, FALSE, TRUE)



p1 = ggplot(data = df_gradient_rank, aes(x = rank, y = weighted.mean.rho, group = nb_exons, colour = nb_exons, alpha = reference)) +
  scale_color_viridis_d() +
  geom_line() +
  geom_point() +
  # geom_line(data = df_gradient_rank, aes(x = rank, y = median.rho), colour = "black", size = 2) +
  xlim(1, 14) +
  xlab("CDS part rank") + ylab("Median recombination rate (ρ/kb)") +
  labs(colour = "# exons") +
  facet_wrap(. ~ species, scales = "free", nrow = 3) +
  scale_alpha_manual(values = c(0.4, 1)) +
  guides(alpha = "none") +
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
p1

df_avg = aggregate(weighted.mean.rho ~ rank + species, exons, median)
df_avg$reference = ifelse(as.character(df_avg$species) %in% list_species, FALSE, TRUE)

p1 = p1 + geom_line(data = df_avg, aes(x = rank, y = weighted.mean.rho, group = "black", colour = "black", alpha = reference), colour = "black", size = 1.5)
p1

# data_subset = data_all[which(data_all$species %in% c("Oryza sativa",
#                                                      "Triticum aestivum",
#                                                      "Spinacia oleracea",
#                                                      "Homo sapiens")),]
# 
# data_subset$species = factor(data_subset$species,
#                              levels = c("Oryza sativa", "Triticum aestivum", "Spinacia oleracea", "Homo sapiens"),
#                              labels = c("Oryza sativa", "Triticum aestivum", "Spinacia oleracea", "Homo sapiens"))
# 
# data_subset = data_subset[which(data_subset$feature == "CDS" & data_subset$nb_exons <= 14),]
# 
# list_genes = data_all$gene_id[which(data_all$feature == "gene" & data_all$nb_snp > 4)]
# 
# table(data_all$species[which(data_all$feature == "gene" & data_all$nb_snp > 4)])
# 
# table(data_all$species[which(data_all$feature == "gene" & data_all$nb_snp > 4)])/table(data_all$species[which(data_all$feature == "gene")])
# 
# data_subset2 = data_subset[which(data_subset$gene_id %in% list_genes),]
# df2 = aggregate(weighted.mean.rho ~ species + nb_exons + rank, data_subset2, mean)
# df2$nb_exons = as.factor(df2$nb_exons)
# 
# df2_avg = aggregate(weighted.mean.rho ~ species + rank, data_subset2, mean)
# 
# p2 = ggplot(data = df2, aes(x = rank, y = weighted.mean.rho, group = nb_exons, colour = nb_exons)) +
#   scale_color_viridis_d() +
#   geom_line() +
#   geom_point() +
#   geom_line(data = df2_avg, aes(x = rank, y = weighted.mean.rho, group = NA, colour = NA), colour = "black", size = 2) +
#   xlim(1, 14) +
#   xlab("CDS part rank") + ylab("Recombination rate (ρ/kb)\n ") +
#   labs(colour = "# exons") +
#   facet_wrap(. ~ species, scales = "free", nrow = 1) +
#   ggtitle("Only genes with at least five SNPs") +
#   theme(plot.title = element_text(color="black", size=fontsize +2, face="bold",hjust = 0.5),
#         axis.title.x = element_text(color="black", size=fontsize),
#         axis.title.y = element_text(color="black", size=fontsize),
#         axis.text=element_text(size=fontsize, colour="black"),
#         axis.text.x=element_text(),
#         title = element_text(color="black", size=fontsize + 2),
#         strip.text.x = element_text(size = fontsize, face = "italic"),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.height = unit(2,"line"),
#         legend.key.width = unit(3,"line"),
#         legend.text=element_text(size=fontsize - 2, face = "italic"),
#         legend.title=element_text(size=fontsize - 2),
#         legend.position='right')
# p2
# 
# 
# data_subset = data_all[which(data_all$species %in% c("Oryza sativa",
#                                                      "Triticum aestivum",
#                                                      "Spinacia oleracea",
#                                                      "Homo sapiens")),]
# 
# data_subset$species = factor(data_subset$species,
#                              levels = c("Oryza sativa", "Triticum aestivum", "Spinacia oleracea", "Homo sapiens"),
#                              labels = c("Oryza sativa", "Triticum aestivum", "Spinacia oleracea", "Homo sapiens"))
# 
# data_subset = data_subset[which(data_subset$feature == "gene" & data_subset$nb_exons <= 14),]
# 
# speciesRho = factor(data_subset$species)
# 
# lab = levels(speciesRho)
# for (i in 1:length(lab)) {
#   lab[i] = mean(data_subset$weighted.mean.rho[which(data_subset$species == lab[i])], na.rm = TRUE)
# }
# data_subset$speciesRho = as.numeric(as.character(factor(speciesRho,
#        levels = levels(speciesRho),
#        labels = lab)))
# 
# lab = levels(speciesRho)
# for (i in 1:length(lab)) {
#   lab[i] = quantile(data_subset$weighted.mean.rho[which(data_subset$species == lab[i])], 0.25, na.rm = TRUE)
# }
# data_subset$speciesRho.quantile25 = as.numeric(as.character(factor(speciesRho,
#        levels = levels(speciesRho),
#        labels = lab)))
# 
# 
# lab = levels(speciesRho)
# for (i in 1:length(lab)) {
#   lab[i] = quantile(data_subset$weighted.mean.rho[which(data_subset$species == lab[i])], 0.5, na.rm = TRUE)
# }
# data_subset$speciesRho.quantile50 = as.numeric(as.character(factor(speciesRho,
#        levels = levels(speciesRho),
#        labels = lab)))
# 
# 
# lab = levels(speciesRho)
# for (i in 1:length(lab)) {
#   lab[i] = quantile(data_subset$weighted.mean.rho[which(data_subset$species == lab[i])], 0.75, na.rm = TRUE)
# }
# data_subset$speciesRho.quantile75 = as.numeric(as.character(factor(speciesRho,
#        levels = levels(speciesRho),
#        labels = lab)))
# 
# data_subset$lowrec = ifelse(data_subset$weighted.mean.rho < data_subset$speciesRho.quantile75, TRUE, FALSE)
# 
# list_genes = data_subset$gene_id[which(data_subset$lowrec == TRUE)]
# 
# table(data_subset$species[which(data_subset$lowrec == TRUE)])
# 
# data_subset2 = data_all[which(data_all$gene_id %in% list_genes & data_all$feature == "CDS" & data_all$nb_exons <= 14),]
# 
# data_subset2 = data_subset2[which(data_subset2$species %in% c("Oryza sativa",
#                                                      "Triticum aestivum",
#                                                      "Spinacia oleracea",
#                                                      "Homo sapiens")),]
# 
# data_subset2$species = factor(data_subset2$species,
#                              levels = c("Oryza sativa", "Triticum aestivum", "Spinacia oleracea", "Homo sapiens"),
#                              labels = c("Oryza sativa", "Triticum aestivum", "Spinacia oleracea", "Homo sapiens"))
# 
# 
# 
# df2 = aggregate(weighted.mean.rho ~ species + nb_exons + rank, data_subset2, mean)
# df2$nb_exons = as.factor(df2$nb_exons)
# 
# df2_avg = aggregate(weighted.mean.rho ~ species + rank, data_subset2, mean)
# 
# p3 = ggplot(data = df2, aes(x = rank, y = weighted.mean.rho, group = nb_exons, colour = nb_exons)) +
#   scale_color_viridis_d() +
#   geom_line() +
#   geom_point() +
#   geom_line(data = df2_avg, aes(x = rank, y = weighted.mean.rho, group = NA, colour = NA), colour = "black", size = 2) +
#   xlim(1, 14) +
#   xlab("CDS part rank") + ylab("Recombination rate (ρ/kb)\n ") +
#   labs(colour = "# exons") +
#   facet_wrap(. ~ species, scales = "free", nrow = 1) +
#   ggtitle("Only the 75% of genes with the lowest recombination rate") +
#   theme(plot.title = element_text(color="black", size=fontsize +2, face="bold",hjust = 0.5),
#         axis.title.x = element_text(color="black", size=fontsize),
#         axis.title.y = element_text(color="black", size=fontsize),
#         axis.text=element_text(size=fontsize, colour="black"),
#         axis.text.x=element_text(),
#         strip.text.x = element_text(size = fontsize, face = "italic"),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.height = unit(2,"line"),
#         legend.key.width = unit(3,"line"),
#         legend.text=element_text(size=fontsize - 2, face = "italic"),
#         legend.title=element_text(size=fontsize - 2),
#         legend.position='right')
# p3
# 
# 
# df_subset = data_all[which(data_all$feature == "CDS" & data_all$nb_exons <= max.exons & data_all$species %in% c("Oryza sativa", "Triticum aestivum", "Spinacia oleracea", "Homo sapiens")),]
# 
# df_subset$species = factor(df_subset$species,
#                              levels = c("Oryza sativa", "Triticum aestivum", "Spinacia oleracea", "Homo sapiens"),
#                              labels = c("Oryza sativa", "Triticum aestivum", "Spinacia oleracea", "Homo sapiens"))
# 
# df_subset$hotoverlap = ifelse(df_subset$hotspot_overlap_intensity2 > 0, TRUE, FALSE)
# # df_subset$hotoverlap = ifelse(df_subset$hotspot_overlap_raw > 0, TRUE, FALSE)
# 
# hotspot_overlap = aggregate(hotoverlap ~ species + rank + nb_exons, df_subset, sum)
# 
# hotspot_overlap$nb_exons = factor(hotspot_overlap$nb_exons)
# 
# hotspot_n = aggregate(hotoverlap ~ species + rank + nb_exons, df_subset, length)
# hotspot_overlap$n_sample = hotspot_n$hotoverlap
# 
# p4 = ggplot(hotspot_overlap, aes(x = rank, y = hotoverlap/n_sample, group = nb_exons, colour = nb_exons)) +
#   geom_point() +
#   geom_line() +
#   scale_colour_viridis_d(option = "B") +
#   xlab("CDS part rank") + ylab("Hotspot overlap\n ") +
#   labs(colour = "# exons") +
#   ggtitle("Proportion of CDS overlapping a hotspot") +
#   facet_wrap(~ species, ncol = 4, scales = "free") +
#   xlim(1, 14) +
#   theme(plot.title = element_text(color="black", size=fontsize +2, face="bold",hjust = 0.5),
#         axis.title.x = element_text(color="black", size=fontsize),
#         axis.title.y = element_text(color="black", size=fontsize),
#         axis.text=element_text(size=fontsize, colour="black"),
#         axis.text.x=element_text(),
#         strip.text.x = element_text(size = fontsize, face = "italic"),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.height = unit(2,"line"),
#         legend.key.width = unit(3,"line"),
#         legend.text=element_text(size=fontsize - 2, face = "italic"),
#         legend.title=element_text(size=fontsize - 2),
#         legend.position='right')
# 



# p = ggarrange(p1, p2, p3, p4, ncol = 1, heights = c(4, 1.5, 1.5, 1.5), labels = "AUTO",
#               common.legend = TRUE, legend = "right") + bgcolor("White") + border(color = "White")
p1
ggsave(file = paste("Figure/Paper/Fig6.jpeg", sep = ""), plot = p1, width = 17, height = 11, dpi = 300)

ggsave(file = paste("Figure/Paper/Fig6.jpeg", sep = ""), plot = p1, width = 17, height = 11, dpi = 600)

# p = ggarrange(p2, p3, p4, ncol = 1, heights = c(1.5, 1.5, 1.5), labels = "AUTO",
#               common.legend = TRUE, legend = "right") + bgcolor("White") + border(color = "White")
# p
# 
# ggsave(file = paste("Figure/Paper/Fig5.jpeg", sep = ""), plot = p, width = 17, height = 11, dpi = 300)
