# Loading env
source("Source/init_knitr.R")
source("Source/set_ggplot_theme.R")


### Table 1 - Species
df_ldmaplength = read.table(file = gzfile(paste("Data/Recombination/LD/ldhat/LD_widths.csv.gz", sep = "")), header = TRUE)

species = levels(as.factor(df_ldmaplength$set))
for (i in 1:length(species)) {
  genus = strsplit(species[i], split = "_")[[1]][1]
  sp = strsplit(species[i], split = "_")[[1]][2]
  species[i] = paste(genus, sp, sep = " ")
}
df_ldmaplength$species = factor(df_ldmaplength$set, levels = levels(as.factor(df_ldmaplength$set)),
                                labels = species)


df_ldmaplength = data.frame(species = df_ldmaplength$species,
                            feature = "ldmap",
                            width = df_ldmaplength$length)

# references
sumstats_species = data.frame(species = sort(unique(as.character(data_all$species))))
sumstats_species$species
sumstats_species$reference = c("The 1001 Genomes Consortium",
                               "Zhang et al. 2021",
                               "Guo et al. 2019",
                               "Yang et al. 2021",
                               "Sudmant et al. 2015",
                               "Sun et al. 2020",
                               "Wang et al. 2018",
                               "Wu et al. 2020",
                               "Liu et al. 2022",
                               "Lozano et al. 2021",
                               "Cai et al. 2021",
                               "Zhou et al. 2020")

sumstats_species$accession = NA
for (i in 1:nrow(sumstats_species)) {
  s = unique(data_all$set[which(data_all$species == sumstats_species$species[i])])
  sumstats_species$accession[i] = unique(chromosome_metadata$accession[which(chromosome_metadata$set == s)])
}

# SNP density (SNP/kb)
# min/max, mean
sumstats_species$n_snp = NA
sumstats_species$genome_length_Mb = NA
sumstats_species$genome_masked_percent = NA
sumstats_species$snp_density_kb = NA
for (i in 1:nrow(sumstats_species)) {
  sumstats_species$n_snp[i] = sum(df_ldmaplength$species == sumstats_species$species[i])
  sumstats_species$genome_length_Mb[i] = sum(df_ldmaplength$width[which(df_ldmaplength$species == sumstats_species$species[i])], na.rm = TRUE)/10^6
  sumstats_species$snp_density_kb[i] = sumstats_species$n_snp[i]/(sumstats_species$genome_length_Mb[i]*10^3)
}

# Number of chromosomes per species
sumstats_species$n_chromosome = NA
for (i in 1:nrow(sumstats_species)) {
  sumstats_species$n_chromosome[i] = length(unique(data_all$chromosome[which(data_all$species == sumstats_species$species[i])]))
}


sumstats_species$n_diploid_genomes_sampled = c(40, 40, 40, 40, 40, 37, 40, 40, 36, 32, 40, 40)


sumstats_species$mating_system = c("Selfer", "Outcrosser", "Mixed", "Selfer",
                                   "Outcrosser", "Outcrosser", "Selfer", "Selfer",
                                   "Outcrosser", "Selfer", "Outcrosser", "Selfer")

# Mean Rho/bp
sumstats_species$mean_rho_kb = NA
sumstats_species$median_rho_kb = NA
for (i in 1:nrow(sumstats_species)) {
  s = unique(data_all$set[which(data_all$species == sumstats_species$species[i])])
  mapfile = paste("Data/Recombination/LD/ldhat/", s, ".csv.gz", sep = "")
  ldmap = read.table(gzfile(mapfile), header = T)
  ldmap$length = ldmap$end - ldmap$start
  mask_length = sum(ldmap$length[which(is.na(ldmap$Mean_rho))], na.rm = TRUE) / 10^6
  sumstats_species$genome_masked_percent[i] = round((mask_length/sumstats_species$genome_length_Mb[i]) * 100, digits = 1)
  sumstats_species$mean_rho_kb[i] = weighted.mean(ldmap$Mean_rho, (ldmap$end - ldmap$start), na.rm = TRUE)
  sumstats_species$median_rho_kb[i] = median(ldmap$Mean_rho, na.rm = TRUE)
}

# Theta
ldmap_metadata = readODS::read_ods(paste("Data/Planning/LDRecMaps_TODO_RUN2.ods", sep =""), sheet = 1)
sumstats_species$theta_bp = NA
for (i in 1:nrow(sumstats_species)) {
  sumstats_species$theta_bp[i] = unique(ldmap_metadata$theta.estimate[which(ldmap_metadata$set == unique(data_all$set[which(data_all$species == sumstats_species$species[i])]))])
}

# Ratio Rho/Theta
sumstats_species$ratio.rho.theta = (sumstats_species$mean_rho_kb/10^3)/sumstats_species$theta_bp

sumstats_species

# Save Table 1
write.xlsx(x = sumstats_species, file = paste("Table/Table1.xls", sep = ""), sheetName = "Summary statistics (species)", row.names = FALSE, append = FALSE)
