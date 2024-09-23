# Loading env
source("Source/init_knitr.R")
source("Source/set_ggplot_theme.R")


fontsize = 16
dotsize = 0.2
linesize = 1.5


df$nb_exons = as.character(df$nb_exons)
df = subset(data_all, data_all$feature %in% c("intron", "exon") & data_all$nb_exons <= max.exons & data_all$rank <= max.exons)

genes_with_intronsUTR = data_all$gene_id[which((data_all$feature == "gene" & data_all$intron3utr == TRUE) | (data_all$feature == "gene" & data_all$intron5utr == TRUE))]
df = df[!(df$gene_id %in% genes_with_intronsUTR),]
df = df[(df$species %in% list_species),]

rho = aggregate(weighted.mean.rho ~ species + rank + feature, df, median)
rho

p1 = ggplot(rho, aes(x = rank, y = weighted.mean.rho, group = feature, colour = feature)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ species, ncol = 3, scales = "free_y") +
  xlab(" \n Exon/intron rank") + ylab("Median ρ/kb") +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        axis.text.x=element_text(size=fontsize - 2),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize, face = "italic"),
        legend.title=element_text(size=fontsize),
        legend.position='right')
p1



# df_gradient_rank = readRDS(file = paste("Data/Recombination/Gradient/meanRho_rank.rds", sep = ""))

# for (i in 1:nrow(df_gradient_rank)) {
#   genus = strsplit(df_gradient_rank$set[i], split = "_")[[1]][1]
#   sp = strsplit(df_gradient_rank$set[i], split = "_")[[1]][2]
#   df_gradient_rank$species[i] = paste(genus, sp, sep = " ")
# }
# # df_gradient_rank = df_gradient_rank[which(df_gradient_rank$species %in% list_species),]

# df_gradient_rank = df_gradient_rank[which(df_gradient_rank$level %in% c("exon", "intron")),]
# df_gradient_rank$level = factor(df_gradient_rank$level, levels = c("exon", "intron"), labels = c("Exon", "Intron"))
# df_gradient_rank$condition = factor(df_gradient_rank$condition, levels = c("observed", "control"), labels = c("Observed", "Control"))


# df_gradient_rank$species = factor(df_gradient_rank$species,
#                                   levels = list_species_ordered)
  

# p1 = ggplot(data = df_gradient_rank, aes(x = rank, y = meanRho, colour = level, linetype = condition)) +
#   geom_line() +
#   xlim(1, 14) +
#   xlab("Exon/intron rank") + ylab("ρ/kb") +
#   labs(colour = "Feature", linetype = "Condition") +
#   facet_wrap(. ~ species, scales = "free", ncol = 4) +
#   scale_color_brewer(palette = "Dark2") +
#   theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
#         axis.title.x = element_text(color="black", size=fontsize),
#         axis.title.y = element_text(color="black", size=fontsize),
#         axis.text=element_text(size=fontsize, colour="black"),
#         axis.text.x=element_text(),
#         strip.text.x = element_text(size = fontsize, face = "italic"),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.height = unit(2,"line"),
#         legend.key.width = unit(3,"line"),
#         legend.text=element_text(size=fontsize, face = "italic"),
#         legend.title=element_text(size=fontsize),
#         legend.position='right')
# p1

ggsave(file = paste("Figure/Paper/FigS6.jpeg", sep = ""), plot = p1, width = 18, height = 12, dpi = 300)
