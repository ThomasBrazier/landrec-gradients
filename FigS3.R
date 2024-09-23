# Loading env
source("Source/init_knitr.R")
source("Source/set_ggplot_theme.R")


df_meanrho = readRDS(paste0("Data/Recombination/weightedmeanrho_per_feature.rds"))

df_meanrho = df_meanrho[-which(df_meanrho$feature %in% c("CDS", "promoter", "mRNA", "buffer_upstream", "buffer_downstream")),]

df_meanrho = df_meanrho[which(df_meanrho$species %in% list_species_ordered),]

# TODO Sort levels of genomic factors
df_meanrho$feature = factor(df_meanrho$feature, levels = c("intergenic", "intergenic_buffer", "gene", "mRNA", "upstream3kb", "upstream2kb", "upstream1kb", "promoter", "utr5", "CDS", "exon", "intron", "utr3", "downstream1kb", "downstream2kb", "downstream3kb"))

df_meanrho$Category = NA
df_meanrho$Category[which(df_meanrho$feature %in% c("utr5", "CDS", "exon", "intron", "utr3"))] = "Intra-genic"
df_meanrho$Category[which(df_meanrho$feature %in% c("upstream3kb", "upstream2kb", "upstream1kb", "downstream1kb", "downstream2kb", "downstream3kb"))] = "Flanking"
df_meanrho$Category[which(df_meanrho$feature %in% c("intergenic"))] = "Intergenic"
df_meanrho$Category[which(df_meanrho$feature %in% c("intergenic_buffer"))] = "Intergenic (buffer)"
df_meanrho$Category[which(df_meanrho$feature %in% c("gene"))] = "Genic"

df_meanrho$Category = factor(df_meanrho$Category, levels = c("Genic", "Intergenic", "Intergenic (buffer)", "Flanking", "Intra-genic"))


df_meanrho$feature = factor(df_meanrho$feature, labels = c("Intergenic", "Intergenic (buffer)", "Genic", "Up 3kb", "Up 2kb", "Up 1kb", "5' UTR", "Exon", "Intron", "3' UTR", "Down 1kb", "Down 2kb", "Down 3kb"))



df_meanrho$species = factor(df_meanrho$species,
                            levels = list_species_ordered)

fontsize = 14
dotsize = 0.2
linesize = 0.6

df_meanrho$Category = factor(df_meanrho$Category, levels = c("Genic", "Intra-genic", "Intergenic", "Intergenic (buffer)", "Flanking"))

# df_meanrho$species_label = gsub(" ", "\n", df_meanrho$species)
df_meanrho$species_label = factor(df_meanrho$species,
                            levels = list_species_ordered,
                            labels = list_species_abbreviations)


p1 = ggplot(data = df_meanrho[which(df_meanrho$Category %in% c("Genic", "Intergenic", "Intergenic (buffer)")),],
            aes(x = species_label, y = wmean, fill = Category)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(x = species_label, ymin = wmean_lower_95, ymax = wmean_upper_95, fill = Category),
                position=position_dodge(width=0.9)) +
  # scale_fill_manual(values = brewer.pal(5, "Dark2")[1:3]) +
  scale_fill_grey(labels = c("Genic", "Intergenic", "Intergenic (buffer)    ")) +
  xlab("Species") + ylab("Weighted mean ρ/kb") +
  facet_wrap(~ species_label, nrow = 1, scale = "free") +
  ggtitle("Genic vs Intergenic recombination rates") +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = fontsize-2, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize, face = "italic"),
        legend.title=element_text(size=fontsize),
        legend.position='right')
# p1


p2 = ggplot(data = df_meanrho[which(!(df_meanrho$Category %in% c("Genic", "Intergenic", "Intergenic (buffer)"))),],
            aes(x = feature, y = wmean, fill = Category)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(x = feature, ymin = wmean_lower_95, ymax = wmean_upper_95)) +
  # scale_fill_manual(values = brewer.pal(5, "Dark2")[4:5]) +
  scale_fill_grey() +
  xlab("") + ylab("Weighted mean ρ/kb") +
  facet_wrap(~ species, ncol = 4, scale = "free_y") +
  ggtitle("Recombination rates within genic regions") +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        axis.text.x=element_text(angle = 90),
        strip.text.x = element_text(size = fontsize-2, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize, face = "italic"),
        legend.title=element_text(size=fontsize),
        legend.position='right')


ggpubr::ggarrange(p1, p2, nrow = 2, labels = "AUTO", heights = c(1.5, 4))


ggsave(file = paste("Figure/Paper/FigS3.jpeg", sep = ""), width = 16, height = 10, dpi = 300)
