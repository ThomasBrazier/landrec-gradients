# Loading env
source("Source/init_knitr.R")
source("Source/set_ggplot_theme.R")


# Figure 4
fontsize = 20
dotsize = 0.2
linesize = 1.5


df_subset = data_all[which(data_all$feature == "CDS" & data_all$nb_exons <= max.exons & data_all$species %in% list_species),]

# df_subset$hotoverlap = ifelse(df_subset$hotspot_overlap_intensity2 > 0, TRUE, FALSE)
df_subset$hotoverlap = ifelse(df_subset$hotspot_overlap_raw > 0, TRUE, FALSE)

hotspot_overlap = aggregate(hotoverlap ~ species + rank + nb_exons, df_subset, sum)

hotspot_overlap$nb_exons = factor(hotspot_overlap$nb_exons)

hotspot_n = aggregate(hotoverlap ~ species + rank + nb_exons, df_subset, length)
hotspot_overlap$n_sample = hotspot_n$hotoverlap


df = subset(data_all, data_all$feature %in% c("intron", "exon") & data_all$nb_exons <= max.exons)
genes_with_intronsUTR = data_all$gene_id[which((data_all$feature == "gene" & data_all$intron3utr == TRUE) | (data_all$feature == "gene" & data_all$intron5utr == TRUE))]
df = df[!(df$gene_id %in% genes_with_intronsUTR),]
df = df[(df$species %in% list_species),]

df$category = factor(paste0(df$feature, "_", df$rank),
                        levels = paste0(rep(c("exon", "intron"), max.exons), "_", rep(1:max.exons, each = 2)))
df$f = ifelse(df$feature == "exon", "e", "i")
df$category_short =  factor(paste0(df$f, df$rank),
                        levels = paste0(rep(c("e", "i"), max.exons), rep(1:max.exons, each = 2)))

df = df[-which(df$category_short == "i14"),]

rho = aggregate(weighted.mean.rho ~ species + category_short, df, median)

rho$feature = NA
for (i in 1:nrow(rho)) {
  rho$feature[i] = ifelse(grepl("e[0-9]+", rho$category_short[i]), "exon", "intron")
}


p1 = ggplot(rho, aes(x = category_short, y = weighted.mean.rho)) +
  geom_line(aes(group = 1)) +
  geom_point(aes(colour = feature)) +
  facet_wrap(~ species, nrow = 3, scales = "free_y") +
  xlab(" \n Exon/intron rank") + ylab("Median ρ/kb") +
  labs(colour = "# exons") +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        axis.text.x=element_text(size=fontsize - 2, angle = 90),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize, face = "italic"),
        legend.title=element_text(size=fontsize),
        legend.position='none')
p1

df = subset(data_all, data_all$feature %in% c("downstream1kb", "upstream1kb", "downstream2kb", "upstream2kb", "downstream3kb", "upstream3kb", "CDS") & data_all$nb_exons <= max.exons)

df$category = as.character(df$feature)
df$category[which(df$feature %in% c("downstream1kb", "downstream2kb", "downstream3kb"))] = "Downstream"
df$category[which(df$feature %in% c("upstream1kb", "upstream2kb", "upstream3kb"))] = "Upstream"

df$rank[which(df$feature %in% c("downstream1kb", "downstream2kb", "downstream3kb"))] = 0
df$rank[which(df$feature %in% c("upstream1kb", "upstream2kb", "upstream3kb"))] = df$nb_exons[which(df$feature %in% c("upstream1kb", "upstream2kb", "upstream3kb"))] + 1

# genes_with_intronsUTR = data_all$gene_id[which((data_all$feature == "gene" & data_all$intron3utr == TRUE) | (data_all$feature == "gene" & data_all$intron5utr == TRUE))]
# df = df[!(df$gene_id %in% genes_with_intronsUTR),]

df = df[(df$species %in% list_species),]

rho = aggregate(weighted.mean.rho ~ species + category + rank + nb_exons, df, median)

rho$nb_exons = factor(rho$nb_exons)

rho$border = NA
# rho$border[which(rho$rank == 1)] = "Downstream"
rho$border[which(rho$category == "Downstream")] = "Downstream"
# rho$border[which(rho$rank == rho$nb_exons)] = "Upstream"
rho$border[which(rho$category == "Upstream")] = "Upstream"
rho$downstream_line = rho$border
rho$downstream_line[is.na(rho$downstream_line) & rho$rank == 1] = "Downstream"
rho$upstream_line = rho$border
rho$upstream_line[is.na(rho$upstream_line) & (rho$rank == rho$nb_exons)] = "Upstream"

p2 = ggplot(rho[which(rho$category %in% c("Upstream", "Downstream")),], aes(x = rank, y = weighted.mean.rho, group = nb_exons)) +
  geom_point() +
  geom_line(data = rho[which(rho$downstream_line == "Downstream"),], aes(x = rank, y = weighted.mean.rho, group = nb_exons)) +
  geom_line(data = rho[which(rho$upstream_line == "Upstream"),], aes(x = rank, y = weighted.mean.rho, group = nb_exons)) +
  facet_wrap(~ species, nrow = 3, scales = "free_y") +
  xlab("CDS part Rank") + ylab("ρ/kb") +
  # geom_point(data = rho[which(rho$category == "Upstream"),], aes(x = rank, y = weighted.mean.rho), colour = "Black") +
  # geom_point(data = rho[which(rho$category == "Downstream"),], aes(x = rank, y = weighted.mean.rho), colour = "Black") +
  geom_point(data = rho[-which(rho$category %in% c("Downstream", "Upstream")),], aes(x = rank, y = weighted.mean.rho, group = nb_exons, colour = nb_exons), alpha = 0.3) +
  geom_line(data = rho[-which(rho$category %in% c("Downstream", "Upstream")),], aes(x = rank, y = weighted.mean.rho, group = nb_exons, colour = nb_exons), alpha = 0.3) +
  scale_color_viridis_d() +
  labs(colour = "# exons") +
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
        legend.position='none')



# df_gradient_rank = readRDS(file = paste("Data/Recombination/Gradient/meanRho_rank_nbexons.rds", sep = ""))
exons = data_all[which(data_all$feature == "CDS" & data_all$nb_exons <= max.exons),]
subsetspecies = c("Arabidopsis thaliana", "Glycine max", "Populus tremula")
df_gradient_rank = aggregate(weighted.mean.rho ~ rank + nb_exons + species, exons[which(exons$species %in% subsetspecies),], median)
df_gradient_rank$nb_exons = as.factor(df_gradient_rank$nb_exons)


# for (i in 1:nrow(df_gradient_rank)) {
#   genus = strsplit(df_gradient_rank$set[i], split = "_")[[1]][1]
#   sp = strsplit(df_gradient_rank$set[i], split = "_")[[1]][2]
#   df_gradient_rank$species[i] = paste(genus, sp, sep = " ")
# }
# df_gradient_rank = df_gradient_rank[which(df_gradient_rank$species %in% c("Arabidopsis thaliana", "Populus tremula", "Glycine max")),]

p3 = ggplot(data = df_gradient_rank, aes(x = rank, y = weighted.mean.rho, group = nb_exons, colour = nb_exons)) +
  scale_color_viridis_d() +
  geom_line() +
  geom_point() +
  geom_line(data = df_gradient_rank[which(df_gradient_rank$nb_exons == "all"),], aes(x = rank, y = weighted.mean.rho), colour = "black", size = 2) +
  xlim(1, 14) +
  xlab("CDS part rank") + ylab("Median recombination rate (ρ/kb)") +
  labs(colour = "# exons") +
  facet_wrap(. ~ species, scales = "free", nrow = 3) +
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
df_avg = aggregate(weighted.mean.rho ~ rank + species, exons[which(exons$species %in% subsetspecies),], median)
p3 = p3 + geom_line(data = df_avg, aes(x = rank, y = weighted.mean.rho, group = "black", colour = "black"), colour = "black", size = 1.5)

p3

pa = ggpubr::ggarrange(p1, p2, nrow = 2, heights = c(3,3), labels = c("B", "C"))
p = ggpubr::ggarrange(p3, pa, ncol = 2, widths = c(8, 4), labels = c("A"))

p

# ggsave(file = paste("Figure/Paper/Fig4.tiff", sep = ""), plot = p, width = 16, height = 18, dpi = 300, compression = "lzw")
ggsave(file = paste("Figure/Paper/Fig4.jpeg", sep = ""), plot = p, width = 16, height = 18, dpi = 300)
ggsave(file = paste("Figure/Paper/Fig4.jpeg", sep = ""), plot = p, width = 18, height = 18, dpi = 600)
