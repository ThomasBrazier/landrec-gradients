# Loading env
source("Source/init_knitr.R")
source("Source/set_ggplot_theme.R")


sumstats_species = read.xlsx(file = paste("Table/Table1.xls", sep = ""), sheetIndex = 1)


### Table 2
ldhot_trimmed_nofiltering = read.ldhot.all(max.length = 10^8,
                               peak.rate = 0,
                               intensity = 0,
                               max.intensity = 10^8)
ldhot_trimmed_soft = read.ldhot.all(max.length = 10^4,
                               peak.rate = 0,
                               intensity = 0,
                               max.intensity = 10^8)
ldhot_trimmed_hard = read.ldhot.all(max.length = 10^4,
                               peak.rate = 0,
                               intensity = 4,
                               max.intensity = 200)

#le nombre de hotspots détectés, la fraction du génome qu’ils occupent et la fraction de la recombinaison totale qui est concentrée dans ces hotspots (et faire ce tableau pour les différents critères de filtrage)
ds = sort(as.character(unique(ldhot_trimmed_hard$dataset)))
sp = ds
for (i in 1:length(sp)) {
  sp[i] = dataset_metadata$species[which(dataset_metadata$dataset == sp[i])]
}

n.hotspots.literature = c('8448', '', '', '', "25000-50000", '', '14125', '', '6000', '', '', '')
ref = c('Choi et al. 2013', '', '', '', 'Myers et al. 2005', '', 'Marand et al. 2019', '', 'Slavov et al. 2012', '', '', '')

df = data.frame(dataset = ds,
                species = sp,
                n.hotspots.raw = NA,
                density.hotspots.raw.Mb = NA,
                # fraction.genome.raw = NA,
                # fraction.recombination.raw = NA,
                n.hotspots.soft = NA,
                percent.hotspots.soft = NA,
                density.hotspots.soft.Mb = NA,
                # fraction.genome.soft = NA,
                # fraction.recombination.soft = NA,
                n.hotspots.hard = NA,
                percent.hotspots.hard = NA,
                density.hotspots.hard.Mb = NA,
                # fraction.genome.hard = NA,
                # fraction.recombination.hard = NA,
                n.hotspots.literature = n.hotspots.literature,
                ref = ref)

# ldhot_trimmed_nofiltering$rho = ldhot_trimmed_nofiltering$estimated_mean_rho * (ldhot_trimmed_nofiltering$length/10^3)
# ldhot_trimmed_soft$rho = ldhot_trimmed_soft$estimated_mean_rho * (ldhot_trimmed_soft$length/10^3)
# ldhot_trimmed_hard$rho = ldhot_trimmed_hard$estimated_mean_rho * (ldhot_trimmed_hard$length/10^3)
  

for (i in 1:nrow(df)) {
  # genome.size = sum(chromosome_metadata$chrsize.bp[which(chromosome_metadata$set == df$dataset[i] & chromosome_metadata$litteralname != "7D")])
  # map = read.table(gzfile(paste0("Data/Recombination/LD/ldhat/", df$dataset[i],".csv.gz")), header = TRUE)
  # genome.size = aggregate(end ~ chromosome, map, max)
  # genome.size = sum(genome.size$end)/10^6
  
  genome.size = sumstats_species$genome_length_Mb[which(sumstats_species$species == df$species[i])] * (1 - sumstats_species$genome_masked_percent[which(sumstats_species$species == df$species[i])]/100)
  
  # map$length = map$end - map$start
  # map$rho = map$Mean_rho * (map$length/10^3)
  # total.rho = sum(map$rho)

  
  df$n.hotspots.raw[i] = nrow(ldhot_trimmed_nofiltering[which(ldhot_trimmed_nofiltering$dataset == df$dataset[i]),])
  df$density.hotspots.raw.Mb[i] = nrow(ldhot_trimmed_nofiltering[which(ldhot_trimmed_nofiltering$dataset == df$dataset[i]),])/genome.size
  # df$fraction.genome.raw[i] = sum(ldhot_trimmed_nofiltering$length[which(ldhot_trimmed_nofiltering$dataset == df$dataset[i])], na.rm = TRUE)/genome.size
  # df$fraction.recombination.raw[i] = sum(ldhot_trimmed_nofiltering$rho[which(ldhot_trimmed_nofiltering$dataset == df$dataset[i])], na.rm = TRUE)/total.rho
  
  df$n.hotspots.soft[i] = nrow(ldhot_trimmed_soft[which(ldhot_trimmed_soft$dataset == df$dataset[i]),])
  df$percent.hotspots.soft[i] = nrow(ldhot_trimmed_soft[which(ldhot_trimmed_soft$dataset == df$dataset[i]),])/df$n.hotspots.raw[i] * 100
  df$density.hotspots.soft.Mb[i] = nrow(ldhot_trimmed_soft[which(ldhot_trimmed_soft$dataset == df$dataset[i]),])/genome.size
  # df$fraction.genome.soft[i] = sum(ldhot_trimmed_soft$length[which(ldhot_trimmed_soft$dbbataset == df$dataset[i])], na.rm = TRUE)/genome.size
  # df$fraction.recombination.soft[i] = sum(ldhot_trimmed_soft$rho[which(ldhot_trimmed_soft$dataset == df$dataset[i])], na.rm = TRUE)/total.rho
  
  df$n.hotspots.hard[i] = nrow(ldhot_trimmed_hard[which(ldhot_trimmed_hard$dataset == df$dataset[i]),])
  df$percent.hotspots.hard[i] = nrow(ldhot_trimmed_hard[which(ldhot_trimmed_hard$dataset == df$dataset[i]),])/df$n.hotspots.raw[i] * 100
  df$density.hotspots.hard.Mb[i] = nrow(ldhot_trimmed_hard[which(ldhot_trimmed_hard$dataset == df$dataset[i]),])/genome.size
  # df$fraction.genome.hard[i] = sum(ldhot_trimmed_hard$length[which(ldhot_trimmed_hard$dataset == df$dataset[i])], na.rm = TRUE)/genome.size
  # df$fraction.recombination.hard[i] = sum(ldhot_trimmed_hard$rho[which(ldhot_trimmed_hard$dataset == df$dataset[i])], na.rm = TRUE)/total.rho
}

df[,-1]
   
knitr::kable(df[,-1],
      format = "html",
      digits = 2,
      align = "c",
      caption = "Number of hotspots. Statistics are given for unfiltered, soft and hard filtering datasets. The number of hotspots already known in the literature is given.",
      col.names = c("Species",
                    "N hotspots (raw)",
                    "Hotspot density (Mb, raw)",
                    "N hotspots (soft)",
                    "% after soft filtering",
                    "Hotspot density (Mb, soft)",
                    "N hotspots (hard)",
                    "% after hard filtering",
                    "Hotspot density (Mb, hard)",
                    "N hotspots (literature)", "Reference"))

write.xlsx(x = df, file = paste("Table/Table2.xls", sep = ""), sheetName = "Summary statistics (hotspots)", row.names = FALSE, append = FALSE)



# The number of hotspots per species was significantly correlated with genome size (Spearman's rank correlation test, $\rho$ = 0.89, $\rho$ = 0.87, $\rho$ = 0.89 for unfiltered, soft and hard filtered datasets, p $<$ 0.001 for the three; Table \ref{table:TableHotspots}) but not significantly correlated with SNP density (Spearman's rank correlation test, $\rho$ = -0.26, p = 0.42, $\rho$ = -0.20, p = 0.51, $\rho$ = -0.24, p = 0.46 for the three datasets, respectively)
colnames(df)
colnames(sumstats_species)

for (i in 1:nrow(df)) {
  df$n_snps[i] = sumstats_species$n_snp[which(sumstats_species$species == df$species[i])]
  df$genome_length_Mb[i] = sumstats_species$genome_length_Mb[which(sumstats_species$species == df$species[i])]
  df$genome_masked_percent[i] = sumstats_species$genome_length_Mb[which(sumstats_species$species == df$species[i])] * (1 - sumstats_species$genome_masked_percent[which(sumstats_species$species == df$species[i])]/100)
}

df$snp_density = df$n_snps / df$genome_length_Mb
colnames(df)

# Nb hotspots ~ genome length
cor.test(df$n.hotspots.raw, df$genome_length_Mb)
cor.test(df$n.hotspots.soft, df$genome_length_Mb)
cor.test(df$n.hotspots.hard, df$genome_length_Mb)

cor.test(df$n.hotspots.raw, df$genome_masked_percent)
cor.test(df$n.hotspots.soft, df$genome_masked_percent)
cor.test(df$n.hotspots.hard, df$genome_masked_percent)

# Nb SNPs ~ Nb hotspots
cor.test(df$n.hotspots.raw, df$n_snps)
cor.test(df$n.hotspots.soft, df$n_snps)
cor.test(df$n.hotspots.hard, df$n_snps)

# SNP density ~ Hotspot density
cor.test(df$n.hotspots.raw, df$snp_density)
cor.test(df$n.hotspots.soft, df$snp_density)
cor.test(df$n.hotspots.hard, df$snp_density)
