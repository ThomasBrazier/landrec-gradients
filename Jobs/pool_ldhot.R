#!/usr/bin/env Rscript
source("Source/init.R")

# # List of dataset
# marey_data = read.table(paste("Data/Recombination/marey_dataset.csv", sep = ""), header = TRUE, sep = "\t")
# chromosome_metadata = read.table(paste0("Data/Genome/genome_chromosome_metadata.csv"),
#                                  header = TRUE, sep = "\t")
# 
# 
chromosome_metadata = chromosome_metadata[which(!is.na(chromosome_metadata$ldmapname)),]

chromosome_metadata = chromosome_metadata[which(chromosome_metadata$set %in% list_dataset),]


# datasets = c("Malus_domestica_Sun2020", "Solanum_lycopersicum_Gupta2020", "Sorghum_bicolor_Lozano2021",
#          "Populus_tremula_Liu2022", "Glycine_max_Yang2021",
#          "Phaseolus_vulgaris_Wu2020", "Camellia_sinensis_Zhang2021", "Citrullus_lanatus_Guo2019",
#          "Spinacia_oleracea_Cai2021", "Arabidopsis_thaliana_1001genomes",
#          "Oryza_sativa_Wang2018", "Triticum_aestivum_Zhou2020")
# 
# source("Source/read.ldhot.R")
# 
# params = read_yaml("parameters.yaml")
# peakrate = params$peakrate
# hotspot.maxlength = params$hotspot.maxlength
# hotspot.intensity = params$hotspot.intensity

#peakrate = 1

# source("Source/read.ldhot.R")
# dataset = "Populus_tremula_Liu2022"
# chromosome = "chr1"
# ldhot = read.ldhot(dataset,
#                    chromosome,
#            max.length = 1000000,
#            peak.rate = 0,
#            detailed.results = T,
#            intensity = 0)


# Raw data
cat("Raw data - Unfiltered hotspots\n")
hotspots = lapply(1:nrow(chromosome_metadata),
                  function(x) {read.ldhot(chromosome_metadata$set[x],
                                          chromosome_metadata$ldmapname[x],
                                          max.length = 1000000,
                                          peak.rate = 0,
                                          detailed.results = T,
                                          intensity = 0)})
hotspots = rbindlist(hotspots)

write.table(hotspots, file = gzfile("Data/Recombination/LD/ldhotspots_raw.csv.gz"), col.names = T,
            row.names = F, sep = "\t", quote = F)

# cat("Hotspots filtered\n")

# Deprecated since read.ldhot.all()

# Hotspot filtering is applied in read.ldhot()
# (1) TODO Merge overlapping hotspots
# (2) Keep only hotspots smaller than 10 kb
# hotspots = lapply(1:nrow(chromosome_metadata),
#               function(x) {read.ldhot(chromosome_metadata$set[x],
#                                       chromosome_metadata$ldmapname[x],
#                                       max.length = hotspot.maxlength,
#                                       peak.rate = peakrate,
#                                       detailed.results = T,
#                                       intensity = 2)})
# hotspots = rbindlist(hotspots)
# 
# write.table(hotspots, file = gzfile("Data/Recombination/LD/ldhotspots_filtered_intensity2.csv.gz"), col.names = T,
#             row.names = F, sep = "\t", quote = F)
# 
# 
# hotspots = lapply(1:nrow(chromosome_metadata),
#                   function(x) {read.ldhot(chromosome_metadata$set[x],
#                                           chromosome_metadata$ldmapname[x],
#                                           max.length = hotspot.maxlength,
#                                           peak.rate = peakrate,
#                                           detailed.results = T,
#                                           intensity = 4)})
# hotspots = rbindlist(hotspots)
# 
# write.table(hotspots, file = gzfile("Data/Recombination/LD/ldhotspots_filtered_intensity4.csv.gz"), col.names = T,
#             row.names = F, sep = "\t", quote = F)

gc()


