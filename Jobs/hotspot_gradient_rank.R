#!/usr/bin/env Rscript

########################################################################## #
#     JOB - Dataset Assembly
########################################################################## #

# This script take the results of LDhat and LDhot
# as well as parsed GFF, in order to compute rho gradients

#============================================================================#
# LOADING ENVIRONMENT ----
#============================================================================#
source("Source/init.R")


#============================================================================#
# Loading variables & objects ----
#============================================================================#
# Get the directory of the file & set working directory
# wd=dirname(rstudioapi::getSourceEditorContext()$path)
# wd=gsub("/Source", "", wd)
# setwd(wd)
# setwd("~/Academic/landrec-gradients")
ncores = 16

# List of dataset
# marey_data = read.table(paste("Data/Recombination/marey_dataset.csv", sep = ""), header = TRUE, sep = "\t")
# chromosome_metadata = read.table(paste0("Data/Genome/genome_chromosome_metadata.csv"),
#                                  header = TRUE, sep = "\t")
# 
# 
# chromosome_metadata = chromosome_metadata[which(!is.na(chromosome_metadata$ldmapname)),]
# 
# sets = c("Malus_domestica_Sun2020", "Solanum_lycopersicum_Gupta2020", "Sorghum_bicolor_Lozano2021",
#          "Populus_tremula_Liu2022", "Glycine_max_Yang2021",
#          "Phaseolus_vulgaris_Wu2020", "Camellia_sinensis_Zhang2021", "Citrullus_lanatus_Guo2019",
#          "Spinacia_oleracea_Cai2021", "Arabidopsis_thaliana_1001genomes",
#          "Oryza_sativa_Wang2018", "Triticum_aestivum_Zhou2020")
# # "Zea_mays_QiSun2018"
# # "Brassica_napus_Wu2019"
# sets2 = c("Arabidopsis_thaliana_1001genomes", "Populus_tremula_Liu2022", "Glycine_max_Yang2021")


# max.exons = 14


### Raw hotspots ----
# Genes grouped by number of exons
rm(df_gradient_rank_nbexons)
rm(gff_recombination)
for (s in list_dataset) {
  cat(s, "\n")
  gff_recombination = data_all[data_all$set == s,]
  # gff_recombination = readRDS(paste("Data/Recombination/Gradient/gff_rho_", s, ".rds", sep = ""))
  gff_recombination$weighted.mean.rho.control = as.numeric(gff_recombination$weighted.mean.rho.control)
  # Filter protein coding genes only
  if ("protein_coding" %in% gff_recombination$gene_biotype) {
    pcgenes = gff_recombination$id[which(gff_recombination$gene_biotype == "protein_coding")]
    gff_recombination = gff_recombination[which(gff_recombination$id %in% pcgenes | gff_recombination$gene_id %in% pcgenes),]
  }
  
  # Filter genes
  # By gene length:
  # - less than 10 kb
  # - less than 15 exons
  # max.exons = 14
  
  gff_recombination$weighted.mean.rho = as.numeric(gff_recombination$weighted.mean.rho)
  gff_recombination$weighted.mean.rho.control = as.numeric(gff_recombination$weighted.mean.rho.control)

  # Subset
  gff_recombination = gff_recombination[which(gff_recombination$feature %in% c("CDS", "exon", "intron") & 
                                                gff_recombination$rank <= max.exons &
                                                gff_recombination$rank > 0),]
  
  hotspot_rank_exons = aggregate(hotspot_overlap_raw ~ rank, data = gff_recombination[which(gff_recombination$feature == "CDS" & gff_recombination$rank <= max.exons & gff_recombination$rank > 0 & gff_recombination$nb_exons <= max.exons),], sum)
  gff_recombination$rank.control = sample(gff_recombination$rank, replace = FALSE)
  hotspot_rank_exons.control = aggregate(hotspot_overlap_raw ~ rank.control, data = gff_recombination[which(gff_recombination$feature == "CDS" & gff_recombination$rank <= max.exons & gff_recombination$rank > 0 & gff_recombination$nb_exons <= max.exons),], sum)
  colnames(hotspot_rank_exons.control)[1] = "rank"
  hotspot_rank_exons_perclass = aggregate(hotspot_overlap_raw ~ rank + nb_exons, data = gff_recombination[which(gff_recombination$feature == "CDS" & gff_recombination$rank <= max.exons & gff_recombination$rank > 0 & gff_recombination$nb_exons <= max.exons),], sum)
  
  hotspot_rank_exons_perclass = hotspot_rank_exons_perclass[which(hotspot_rank_exons_perclass$nb_exons <= max.exons),]
  hotspot_rank_exons_perclass$n_sample = NA
  for (i in 1:nrow(hotspot_rank_exons_perclass)) {
    hotspot_rank_exons_perclass$n_sample[i] = nrow(gff_recombination[which(gff_recombination$feature == "CDS" & gff_recombination$rank == hotspot_rank_exons_perclass$rank[i] & gff_recombination$nb_exons == hotspot_rank_exons_perclass$nb_exons[i]),])
  }
  
  hotspot_rank_exons$nb_exons = "all"
  for (i in 1:nrow(hotspot_rank_exons)) {
    hotspot_rank_exons$n_sample[i] = nrow(gff_recombination[which(gff_recombination$feature == "CDS" & gff_recombination$rank == hotspot_rank_exons$rank[i]),])
  }
  hotspot_rank_exons.control$nb_exons = "control"
  for (i in 1:nrow(hotspot_rank_exons.control)) {
    hotspot_rank_exons.control$n_sample[i] = nrow(gff_recombination[which(gff_recombination$feature == "CDS" & gff_recombination$rank.control == hotspot_rank_exons$rank[i]),])
  }
  # hotspot_rank_exons$n_sample = nrow(gff_recombination[which(gff_recombination$feature == "CDS" & gff_recombination$rank <= max.exons & gff_recombination$rank > 0 & gff_recombination$nb_exons <= max.exons),])
  
  hotspot_rank_exons_perclass = rbind.fill(hotspot_rank_exons_perclass,
                                           hotspot_rank_exons,
                                           hotspot_rank_exons.control)
  hotspot_rank_exons_perclass$nb_exons = factor(hotspot_rank_exons_perclass$nb_exons, levels = c(1:15, "all", "control"))
  
  hotspot_rank_exons_perclass$set = s
  
  if (exists("df_gradient_rank_nbexons")) {
    df_gradient_rank_nbexons = rbind(df_gradient_rank_nbexons, hotspot_rank_exons_perclass)
    rm(hotspot_rank_exons_perclass)
  } else {
    df_gradient_rank_nbexons = hotspot_rank_exons_perclass
    rm(hotspot_rank_exons_perclass)
  }
}
# df_gradient_rank_nbexons

saveRDS(df_gradient_rank_nbexons, file = paste("Data/Recombination/Hotspot/hotspot_gradient_rank.rds", sep = ""))
rm(df_gradient_rank_nbexons)

### Trimmed hotspots ----
# Genes grouped by number of exons
rm(df_gradient_rank_nbexons)
rm(gff_recombination)
for (s in list_dataset) {
  cat(s, "\n")
  # gff_recombination = readRDS(paste("Data/Recombination/Gradient/gff_rho_", s, ".rds", sep = ""))
  gff_recombination = data_all[data_all$set == s,]
  gff_recombination$weighted.mean.rho.control = as.numeric(gff_recombination$weighted.mean.rho.control)
  # Filter protein coding genes only
  if ("protein_coding" %in% gff_recombination$gene_biotype) {
    pcgenes = gff_recombination$id[which(gff_recombination$gene_biotype == "protein_coding")]
    gff_recombination = gff_recombination[which(gff_recombination$id %in% pcgenes | gff_recombination$gene_id %in% pcgenes),]
  }
  
  # Filter genes
  # By gene length:
  # - less than 10 kb
  # - less than 15 exons
  # max.exons = 14
  
  gff_recombination$weighted.mean.rho.control = as.numeric(gff_recombination$weighted.mean.rho.control)
  
  
  # Subset
  gff_recombination = gff_recombination[which(gff_recombination$feature %in% c("CDS", "exon", "intron") & 
                                                gff_recombination$rank <= max.exons &
                                                gff_recombination$rank > 0),]
  
  hotspot_rank_exons = aggregate(hotspot_overlap_intensity4 ~ rank, data = gff_recombination[which(gff_recombination$feature == "CDS" & gff_recombination$rank <= max.exons & gff_recombination$rank > 0),], sum)
  
  hotspot_rank_exons_perclass = aggregate(hotspot_overlap_intensity4 ~ rank + nb_exons, data = gff_recombination[which(gff_recombination$feature == "CDS" & gff_recombination$rank <= max.exons & gff_recombination$rank > 0),], sum)
  
  hotspot_rank_exons_perclass = hotspot_rank_exons_perclass[which(hotspot_rank_exons_perclass$nb_exons <= max.exons),]
  hotspot_rank_exons_perclass$n_sample = NA
  for (i in 1:nrow(hotspot_rank_exons_perclass)) {
    hotspot_rank_exons_perclass$n_sample[i] = nrow(gff_recombination[which(gff_recombination$feature == "CDS" & gff_recombination$rank == hotspot_rank_exons_perclass$rank[i] & gff_recombination$nb_exons == hotspot_rank_exons_perclass$nb_exons[i]),])
  }
  
  hotspot_rank_exons$nb_exons = "all"
  for (i in 1:nrow(hotspot_rank_exons)) {
    hotspot_rank_exons$n_sample[i] = nrow(gff_recombination[which(gff_recombination$feature == "CDS" & gff_recombination$rank == hotspot_rank_exons$rank[i]),])
  }
  # hotspot_rank_exons$n_sample = nrow(gff_recombination[which(gff_recombination$feature == "CDS" & gff_recombination$rank <= max.exons & gff_recombination$rank > 0 & gff_recombination$nb_exons <= max.exons),])
  
  
  hotspot_rank_exons_perclass = rbind(hotspot_rank_exons_perclass,
                                      hotspot_rank_exons)
  hotspot_rank_exons_perclass$nb_exons = factor(hotspot_rank_exons_perclass$nb_exons, levels = c(1:15, "all"))
  
  hotspot_rank_exons_perclass$set = s
  
  if (exists("df_gradient_rank_nbexons")) {
    df_gradient_rank_nbexons = rbind(df_gradient_rank_nbexons, hotspot_rank_exons_perclass)
    rm(hotspot_rank_exons_perclass)
  } else {
    df_gradient_rank_nbexons = hotspot_rank_exons_perclass
    rm(hotspot_rank_exons_perclass)
  }
}
# df_gradient_rank_nbexons

saveRDS(df_gradient_rank_nbexons, file = paste("Data/Recombination/Hotspot/hotspot_gradient_rank_filtered.rds", sep = ""))
rm(df_gradient_rank_nbexons)
