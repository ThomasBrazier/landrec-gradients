# Loading env
source("Source/init_knitr.R")
source("Source/set_ggplot_theme.R")


fontsize = 16
dotsize = 0.2
linesize = 1.5



# A. Rho estimates do not depend on SNP density or Pi
# 1. Gradient SNP Density
# 2. Gradient Rho/Pi
# "Arabidopsis_thaliana_1001genomes", "Populus_tremula_Liu2022", "Glycine_max_Yang2021"
pi_Athaliana = readRDS(paste0("Data/Genomic_landscapes/Polymorphism/gff_pi_Arabidopsis_thaliana_1001genomes.rds"))
pi_Ptremula = readRDS(paste0("Data/Genomic_landscapes/Polymorphism/gff_pi_Populus_tremula_Liu2022.rds"))
pi_Gmax = readRDS(paste0("Data/Genomic_landscapes/Polymorphism/gff_pi_Glycine_max_Yang2021.rds"))
pi_Athaliana$species = "Arabidopsis thaliana"
pi_Ptremula$species = "Populus tremula"
pi_Gmax$species = "Glycine max"
pi_data = rbind(pi_Athaliana, pi_Ptremula, pi_Gmax)

pi_CDS = pi_data[which(pi_data$feature == "CDS" & pi_data$rank > 0 & pi_data$nb_exons <= 14),]

pi_CDS$snp_density_bp = pi_CDS$n_snp/pi_CDS$width

pi_gradient = aggregate(pi_bp ~ nb_exons + rank + species, pi_CDS, mean)

rho_gradient = aggregate(weighted.mean.rho ~ nb_exons + rank + species, pi_CDS, median)

snp_gradient = aggregate(snp_density_bp ~ nb_exons + rank + species, pi_CDS, mean)

pirho_gradient = rho_gradient
pirho_gradient$pi_bp = pi_gradient$pi_bp
pirho_gradient$rho_pi = pirho_gradient$weighted.mean.rho/(pirho_gradient$pi_bp*1000)

pi_gradient$rho = rho_gradient$weighted.mean.rho
snp_gradient$rho = rho_gradient$weighted.mean.rho

rhosnp_gradient = rho_gradient
rhosnp_gradient$snp_density_bp = snp_gradient$snp_density_bp
rhosnp_gradient$rho_snp = rhosnp_gradient$weighted.mean.rho/(rhosnp_gradient$snp_density_bp*1000)


snp_gradient$nb_exons = as.factor(snp_gradient$nb_exons)
p1a = ggplot(data = snp_gradient, aes(x = rank, y = snp_density_bp, group = nb_exons, colour = nb_exons)) +
  geom_point() +
  geom_line() +
  xlab("CDS part rank") +
  scale_color_viridis_d() +
  ylab("SNP density (kb)") +
  labs(colour = "# exons") +
  facet_wrap(~ species, ncol = 3, scales = "free_y") +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize, face = "italic"),
        legend.title=element_text(size=fontsize),
        legend.position='none')
p1a


pirho_gradient$nb_exons = as.factor(pirho_gradient$nb_exons)
p1b = ggplot(data = pirho_gradient, aes(x = rank, y = rho_pi, group = nb_exons, colour = nb_exons)) +
  geom_point() +
  geom_line() +
  xlab("CDS part rank") +
  scale_color_viridis_d() +
  ylab("ρ/π") +
  labs(colour = "# exons") +
  facet_wrap(~ species, ncol = 3, scales = "free_y") +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize, face = "italic"),
        legend.title=element_text(size=fontsize),
        legend.position='bottom')
p1b

p1 = ggpubr::ggarrange(p1a, p1b, nrow = 2, labels = c("A", "B"), heights = c(2, 3))
p1



# B. We are interested in relative rates, variation along the genic sequence
# Gradients in relative rate
# Redo with relative recombination rate instead of absolute one
# exons = data_all[which(data_all$feature == "CDS" & data_all$nb_exons <= max.exons),]
# subsetspecies = c("Arabidopsis thaliana", "Glycine max", "Populus tremula")
# 
# df = aggregate(relative.meanrho ~ rank + nb_exons + species, exons[which(exons$species %in% subsetspecies),], median)
# exons$median.position = exons$dist_atg + exons$width/2
# df_exonPos = aggregate(median.position ~ rank + nb_exons + species, exons[which(exons$species %in% subsetspecies),], median)
# df$nb_exons = as.factor(df$nb_exons)
# df$median.position = df_exonPos$median.position
# 
# p2 = ggplot(df, aes(x = median.position, y = relative.meanrho, group = nb_exons, color = nb_exons)) +
#   geom_point() +
#   geom_line() +
#   xlab("Exon median position within the gene (distance to ATG)") + ylab("Relative rho/kb (exon rho/gene rho)") +
#   scale_color_viridis_d() +
#   facet_wrap(~ species, ncol = 3, scales = "free") +
#   theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
#         axis.title.x = element_text(color="black", size=fontsize),
#         axis.title.y = element_text(color="black", size=fontsize),
#         axis.text=element_text(size=fontsize, colour="black"),
#         axis.text.x=element_text(size=fontsize - 2, angle = 90),
#         strip.text.x = element_text(size = fontsize, face = "italic"),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.height = unit(2,"line"),
#         legend.key.width = unit(3,"line"),
#         legend.text=element_text(size=fontsize, face = "italic"),
#         legend.title=element_text(size=fontsize),
#         legend.position='none')
# p2
# 
# p2 = ggpubr::ggarrange(p2, labels = "C")



# C. Not an artefact of the LD-based method
# 1. Rowan_CO_hotspot_25kb_COcount CO count around hotspot centre
# 2. Rowan_CO_hotspot_5kb cM/Mb around hotspot centre
# 3. Mean CO count around TSS
# 4. Mean CO count around TTS
# 5. Mean CO count per exon (spanning more COs)

# Rowan et al data show the same pattern
rowan = read.table("Data/Rowan_2019/FileS2r1.csv", header = T, sep = ",")
rowan_range = GRanges(seqnames = rowan$chr, range = IRanges(start = rowan$block1.end.pos, end = rowan$block2.start.pos), strand = "*")
rowan_range

# Genome wide CO count in 200bp
chrsize = chromosome_metadata$chrsize.bp[which(chromosome_metadata$set == "Arabidopsis_thaliana_1001genomes")]
chrsize
intervals = lapply(chrsize, function(x) {seq(1, x, by = 200)})
nb_intervals = unlist(lapply(intervals, length))
chr = c(rep("1", nb_intervals[1]),
        rep("2", nb_intervals[2]),
        rep("3", nb_intervals[3]),
        rep("4", nb_intervals[4]),
        rep("5", nb_intervals[5]))
co_ranges = GRanges(seqnames = chr, range = IRanges(start = unlist(intervals), end = unlist(intervals) + 199), strand = "*")
overlaps = countOverlaps(co_ranges, rowan_range, type = "any")
summary(overlaps)
GwCO = mean(overlaps)

gff = readRDS("Data/Genomic_landscapes/Rho/gff_rho_Arabidopsis_thaliana_1001genomes.rds")
gff = gff[which(gff$feature == "CDS" & gff$nb_exons <= 14),]
gff_range = makeGRangesFromDataFrame(gff, keep.extra.columns = T)



# Count how many crossovers in each exon ----
# overlaps = countOverlaps(gff_range, rowan_range, type = "any")
# gff_range$rowan_CO = overlaps
# summary(gff_range$rowan_CO)

# # Aggregate mean CO count ----
# df = aggregate(rowan_CO ~ rank + nb_exons, as.data.frame(gff_range), mean)
# df_samplesize = aggregate(rowan_CO ~ rank + nb_exons, as.data.frame(gff_range), length)

# df$nb_exons = as.factor(df$nb_exons)

# # plot CO gradients
# p4a = ggplot2::ggplot(df, aes(x = rank, y = rowan_CO, group = nb_exons, color = nb_exons)) +
#   geom_point() +
#   geom_line() +
#   xlab("Exon rank") + ylab("CO count per exon") +
#   scale_color_viridis_d() +
#   labs(colour = "# exons") +
#   theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
#         axis.title.x = element_text(color="black", size=fontsize),
#         axis.title.y = element_text(color="black", size=fontsize),
#         axis.text=element_text(size=fontsize, colour="black"),
#         strip.text.x = element_text(size = fontsize, face = "italic"),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.height = unit(2,"line"),
#         legend.key.width = unit(3,"line"),
#         legend.text=element_text(size=fontsize, face = "italic"),
#         legend.title=element_text(size=fontsize),
#         legend.position='none')
# p4a






# CO hotspots ----
# (1) Are Rowan's COs hotspots of rho?
# (2) How many Rowan COs and LDhat hotspots are congruent/shared?
windows = 100
interval = seq(-5000, 5000, by = windows)
interval

# # replicate each row/CO to make interval around
# rowan_hotspots = rowan[rep(seq_len(nrow(rowan)), each = length(interval)), ]

# # Add index of position around CO centre
# rowan_hotspots$index = rep(interval, nrow(rowan))

# # Add physical position of the window (start and end)
# rowan_hotspots$start = rowan_hotspots$breakpoint.pos + rowan_hotspots$index - windows/2 -1
# rowan_hotspots$end = rowan_hotspots$breakpoint.pos + rowan_hotspots$index + windows/2
# colnames(rowan_hotspots)[1] = "chromosome"
# rowan_hotspots = makeGRangesFromDataFrame(rowan_hotspots, keep.extra.columns = T)

# # Estimate rho in windows
# ldmap = read.table(gzfile("Data/Recombination/LD/ldhat/Arabidopsis_thaliana_1001genomes.csv.gz"),
#                    header = T)
# ldmap = makeGRangesFromDataFrame(ldmap, keep.extra.columns = T)

# overlap = findOverlaps(rowan_hotspots, ldmap)
# overlap
# # rowan_hotspots[1]
# # ldmap[591]
# # overlap = overlap[1:100000,]
# meanRho = pbmclapply(1:length(rowan_hotspots),
#                      function(x) mean(ldmap$Mean_rho[subjectHits(overlap)[which(queryHits(overlap) == x)]]),
#                      mc.cores = 1)
# rowan_hotspots$meanRho = unlist(meanRho)

# SAVE and RELOAD
# saveRDS(rowan_hotspots, file = "Data/Rowan_2019/rowan_hotspots.rds")

rowan_hotspots = readRDS("Data/Rowan_2019/rowan_hotspots.rds")

# Plot rho ~ distance to CO centre
df = aggregate(meanRho ~ index, as.data.frame(rowan_hotspots), function(x) {mean(x, na.rm = T)})

p4a1 = ggplot(df, aes(x = index/1000, y = meanRho)) +
  geom_line() +
  xlab("Distance to \nCO centre (kb)") + ylab("ρ/kb") +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize, face = "italic"),
        legend.title=element_text(size=fontsize),
        legend.position='none')
p4a1




# CO count ~ distance to hotspot centre ----
# The cM/Mb map is not available, but I can count the number of COs
# Are LD-based hotspots truly CO hotspots?
windows = 500
interval = seq(-25000, 25000, by = windows)
interval

# LD-based Hotspots
# replicate each row/CO to make interval around
# Soft filtering
source("Source/read.ldhot.R")
ldhot = read.ldhot.all(max.length = 10000,
                       peak.rate = 0,
                       intensity = 0,
                       max.intensity = 10^8)
ldhot = ldhot[which(ldhot$dataset == "Arabidopsis_thaliana_1001genomes"),]
ld_hotspots = ldhot[rep(seq_len(nrow(ldhot)), each = length(interval)), ]

# Add index of position around hotspot centre
ld_hotspots$index = rep(interval, nrow(ldhot))

# Add physical position of the window (start and end)
ld_hotspots$start.hot = ld_hotspots$start
ld_hotspots$end.hot = ld_hotspots$end
ld_hotspots$start = ld_hotspots$start.hot + ld_hotspots$index - windows/2
ld_hotspots$end = ld_hotspots$start.hot + ld_hotspots$index + windows/2 - 1
ld_hotspots = makeGRangesFromDataFrame(ld_hotspots, keep.extra.columns = T)
ld_hotspots

# Count COs in windows
overlap = countOverlaps(ld_hotspots, rowan_range)
overlap
ld_hotspots$CO_count = as.numeric(overlap)

# Plot CO count ~ distance to LD-based hotspot centre
df = aggregate(CO_count ~ index, as.data.frame(ld_hotspots), function(x) {mean(x, na.rm = T)})

p4a2 = ggplot(df, aes(x = index/1000, y = CO_count)) +
  geom_line() +
  xlab("Distance to LD-based\nhotspot centre (kb)") + ylab("Mean CO count") +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize, face = "italic"),
        legend.title=element_text(size=fontsize),
        legend.position='none')
p4a2


# Spatial gradient around TSS/TTS ----
# CO count in sliding windows of 200 bp
# TSS ----
windows = 200
interval = seq(-20000, 20000, by = windows)

gff = readRDS("Data/Genomic_landscapes/Rho/gff_rho_Arabidopsis_thaliana_1001genomes.rds")
genes = gff[which(gff$feature == "mRNA" & gff$nb_exons <= 14 & gff$strand == "+"),]

# replicate each row/gene to make interval around
genes_range = genes[rep(seq_len(nrow(genes)), each = length(interval)), ]

# Add index of physical distance around TSS
genes_range$index = rep(interval, nrow(genes))

# Add physical position of the window (start and end)
genes_range$start.mRNA = genes_range$start
genes_range$end.mRNA = genes_range$end
genes_range$start = genes_range$start.mRNA + genes_range$index - windows/2 -1
genes_range$end = genes_range$start.mRNA + genes_range$index + windows/2

genes_range = makeGRangesFromDataFrame(genes_range, keep.extra.columns = T)
genes_range

# Count COs in windows
overlap = countOverlaps(genes_range, rowan_range)
overlap
genes_range$CO_count = as.numeric(overlap)
df = aggregate(CO_count ~ index, as.data.frame(genes_range), function(x) {mean(x, na.rm = T)})


# # Estimate CO count in windows
# overlap = findOverlaps(genes_range, rowan_range, type = "any",
#                        ignore.strand = T, minoverlap = windows*0.8)
# tab = table(queryHits(overlap))

# table(tab)
# hist(tab)

# genes_range$CO_count = NA
# genes_range$CO_count[as.numeric(names(tab))] = as.numeric(tab)

# # Plot rho ~ distance to CO centre
# df = aggregate(CO_count ~ index, as.data.frame(genes_range), sum)
# df_len = aggregate(CO_count ~ index, as.data.frame(genes_range), length)

# # Relative CO count
# df$CO_count_relative = df$CO_count/(GwCO*mean(df_len$CO_count))


p4b = ggplot(df, aes(x = index/1000, y = CO_count)) +
  geom_line() +
  geom_smooth(method = "loess", se = T) +
  # ylim(5700, 6100) +
  xlab("Distance to TSS (kb)") + ylab("Mean CO count") +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize, face = "italic"),
        legend.title=element_text(size=fontsize),
        legend.position='none')
p4b


# TTS ----
windows = 200
interval = seq(-20000, 20000, by = windows)

gff = readRDS("Data/Genomic_landscapes/Rho/gff_rho_Arabidopsis_thaliana_1001genomes.rds")
genes = gff[which(gff$feature == "mRNA" & gff$nb_exons <= 14 & gff$strand == "+"),]

# replicate each row/gene to make interval around
genes_range = genes[rep(seq_len(nrow(genes)), each = length(interval)), ]

# Add index of physical distance around TSS
genes_range$index = rep(interval, nrow(genes))

# Add physical position of the window (start and end)
genes_range$start.mRNA = genes_range$start
genes_range$end.mRNA = genes_range$end
genes_range$start = genes_range$end.mRNA + genes_range$index - windows/2 -1
genes_range$end = genes_range$end.mRNA + genes_range$index + windows/2

genes_range = makeGRangesFromDataFrame(genes_range, keep.extra.columns = T)
genes_range


# Count COs in windows
overlap = countOverlaps(genes_range, rowan_range)
overlap
genes_range$CO_count = as.numeric(overlap)
df = aggregate(CO_count ~ index, as.data.frame(genes_range), function(x) {mean(x, na.rm = T)})



# # Estimate CO count in windows
# overlap = findOverlaps(genes_range, rowan_range, type = "any",
#                        ignore.strand = T, minoverlap = windows*0.8)
# tab = table(queryHits(overlap))

# table(tab)
# hist(tab)

# genes_range$CO_count = NA
# genes_range$CO_count[as.numeric(names(tab))] = as.numeric(tab)

# # Plot rho ~ distance to CO centre
# df = aggregate(CO_count ~ index, as.data.frame(genes_range), sum)
# df_len = aggregate(CO_count ~ index, as.data.frame(genes_range), length)

# # Relative CO count
# df$CO_count_relative = df$CO_count/(GwCO*mean(df_len$CO_count))

p4c = ggplot(df, aes(x = index/1000, y = CO_count)) +
  geom_line() +
  geom_smooth(method = "loess", se = T) +
  # ylim(5700, 6100) +
  xlab("Distance to TTS (kb)") + ylab("Mean CO count") +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize, face = "italic"),
        legend.title=element_text(size=fontsize),
        legend.position='none')
p4c


p4 = ggpubr::ggarrange(p4a1, p4a2, p4b, p4c, ncol = 4, widths = c(1,1,1,1), labels = c("C", "D", "E", "F"))

p4


p = ggpubr::ggarrange(p1, p4, nrow = 2, heights = c(2.5,1))

p

ggsave(file = paste("Figure/Paper/Fig5.jpeg", sep = ""), plot = p, width = 14, height = 10, dpi = 300)
ggsave(file = paste("Figure/Paper/Fig5.jpeg", sep = ""), plot = p, width = 14, height = 12, dpi = 600)
