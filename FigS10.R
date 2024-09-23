rm(list=ls(all=TRUE))
source("Source/init_knitr.R")
source("Source/set_ggplot_theme.R")

gc()


fontsize = 16
dotsize = 0.1
linesize = 0.7

gene = data_all[which(data_all$feature == "gene" & data_all$nb_exons < 15),]
exons = data_all[which(data_all$feature == "CDS" & data_all$nb_exons < 15 & data_all$rank < 15),]
introns = data_all[which(data_all$feature == "intron" & data_all$nb_exons < 15 & data_all$rank < 15),]



df_gradient_rank = aggregate(weighted.mean.rho ~ rank + nb_exons + species, exons, median)
df_gradient_rank$width = aggregate(width ~ rank + nb_exons + species, exons, mean)$width
df_gradient_rank$nb_exons = as.factor(df_gradient_rank$nb_exons)

df_gradient_rank$species_abbrev = factor(df_gradient_rank$species,
                            levels = list_species_ordered,
                            labels = list_species_abbreviations)


p1 = ggplot(data = df_gradient_rank, aes(x = rank, y = width, group = nb_exons, colour = nb_exons)) +
  scale_color_viridis_d() +
  geom_line() +
  geom_point() +
  xlim(1, 14) +
  xlab("CDS part rank") + ylab("Exon length (bp)") +
  labs(colour = "# exons") +
  facet_wrap(. ~ species_abbrev, scales = "free", nrow = 3) +
  ggtitle("Exon length") +
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
        legend.position='right')
p1

p1a = ggplot(data = df_gradient_rank[which(df_gradient_rank$species %in% c(list_species, "Homo sapiens")),], aes(x = rank, y = width, group = nb_exons, colour = nb_exons)) +
  scale_color_viridis_d() +
  geom_line() +
  geom_point() +
  xlim(1, 14) +
  xlab("CDS part rank") + ylab("Exon length (bp)") +
  labs(colour = "# exons") +
  facet_wrap(. ~ species, scales = "free", nrow = 1) +
  ggtitle("Exon length") +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        axis.text.x=element_text(),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(1,"line"),
        legend.key.width = unit(1,"line"),
        legend.text=element_text(size=fontsize - 6, face = "italic"),
        legend.title=element_text(size=fontsize - 6),
        legend.position='right')
p1a

# introns$nb_introns = introns$nb_exons - 1
# df_gradient_introns = aggregate(weighted.mean.rho ~ rank + nb_introns + species, introns, median)
introns$keep = (introns$rank < introns$nb_exons + 2)
introns_subset = introns[introns$keep,]

table(gene$intron5utr, gene$species)
table(gene$intron3utr, gene$species)

# Remove genes with introns in UTRs
blacklist = c(gene$gene_id[which(gene$intron5utr == "TRUE")],
                gene$gene_id[which(gene$intron3utr == "TRUE")])
length(blacklist)

introns_subset = introns[-which(introns$gene_id %in% blacklist),]
nrow(introns_subset)
introns_subset$keep = (introns_subset$rank < (introns_subset$nb_exons))
introns_subset = introns_subset[introns_subset$keep,]
nrow(introns_subset)



df_gradient_introns = aggregate(width ~ rank + nb_exons + species, introns_subset, mean)
df_gradient_introns$nb_exons = as.factor(df_gradient_introns$nb_exons)
df_gradient_introns$width[which(df_gradient_introns$species == "Citrullus lanatus")] = NA

# df_gradient_introns_avg = aggregate(width ~ rank + species, introns_subset, mean)
# df_gradient_introns_avg$width[which(df_gradient_introns_avg$species == "Citrullus lanatus")] = NA
# Sample size
table(introns_subset$nb_exons, introns_subset$rank, introns_subset$species)

df_gradient_introns$species_abbrev = factor(df_gradient_introns$species,
                            levels = list_species_ordered,
                            labels = list_species_abbreviations)
                            
p2 = ggplot(data = df_gradient_introns, aes(x = rank, y = width, group = nb_exons, colour = nb_exons)) +
  scale_color_viridis_d() +
  geom_line() +
  geom_point() +
  xlim(1, 14) +
  xlab("Intron rank") + ylab("Intron length (bp)") +
  labs(colour = "# exons") +
  facet_wrap(. ~ species_abbrev, scales = "free", nrow = 3) +
  ggtitle("Intron length") +
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
        legend.position='right')
p2



p2a = ggplot(data = df_gradient_introns[which(df_gradient_introns$species %in% c(list_species, "Homo sapiens")),], aes(x = rank, y = width, group = nb_exons, colour = nb_exons)) +
  scale_color_viridis_d() +
  geom_line() +
  geom_point() +
  xlim(1, 14) +
  xlab("Intron rank") + ylab("Intron length (bp)") +
  labs(colour = "# exons") +
  facet_wrap(. ~ species, scales = "free", nrow = 1) +
  ggtitle("Intron length") +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        axis.text.x=element_text(),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(1,"line"),
        legend.key.width = unit(1,"line"),
        legend.text=element_text(size=fontsize - 6, face = "italic"),
        legend.title=element_text(size=fontsize - 6),
        legend.position='right')
p2a


ggarrange(p1a, p2a, nrow = 2, labels = "AUTO", heights = c(2, 2))

# ggsave(file = paste("Figure/Paper/Fig7.jpeg", sep = ""), width = 12, height = 8, dpi = 300)


ggarrange(p1, p2, nrow = 2, labels = "AUTO")

ggsave(file = paste("Figure/Paper/FigS10.jpeg", sep = ""), width = 20, height = 12, dpi = 300)
