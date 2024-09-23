#!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# if (length(args) != 2) {
#   stop("<wd> and <ncores> must be supplied.\n", call.=FALSE)
# }
# wd = args[1]
# setwd(wd)
# ncpus = args[2]

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
# ncores = min(ncpus, detectCores())
ncores = 8

# List of dataset
# marey_data = read.table(paste("Data/Recombination/marey_dataset.csv", sep = ""), header = TRUE, sep = "\t")
# chromosome_metadata = read.table(paste0("Data/Genome/genome_chromosome_metadata.csv"),
#                                  header = TRUE, sep = "\t")
# 
# 
# chromosome_metadata = chromosome_metadata[which(!is.na(chromosome_metadata$ldmapname)),]

# sets = c("Malus_sieversii_Sun2020", "Sorghum_bicolor_Lozano2021",
#          "Populus_tremula_Liu2022", "Glycine_max_Yang2021", "Homo_sapiens_Sudmant2015",
#          "Phaseolus_vulgaris_Wu2020", "Camellia_sinensis_Zhang2021", "Citrullus_lanatus_Guo2019",
#          "Spinacia_oleracea_Cai2021", "Arabidopsis_thaliana_1001genomes",
#          "Oryza_sativa_Wang2018", "Triticum_aestivum_Zhou2020")
# "Zea_mays_QiSun2018"
# "Brassica_napus_Wu2019"
# sets2 = c("Arabidopsis_thaliana_1001genomes", "Populus_tremula_Liu2022", "Glycine_max_Yang2021")


#============================================================================#
# Import data from cluster ----
#============================================================================#



#============================================================================#
# Pool LD maps ----
#============================================================================#
cat("======================================\n")
cat("Pool LD maps...\n")
# source("Jobs/pool_ldmaps.R")


#============================================================================#
# Pool LD hotspots ----
#============================================================================#
cat("======================================\n")
cat("Pool LD hotspots...\n")
# source("Jobs/pool_ldhot.R")



#============================================================================#
# GFF Rho ----
#============================================================================#
#  TODO Pool GFF rho per species
cat("======================================\n")
cat("GFF rho\n")
# source("Jobs/pool_gffrho.R")



#============================================================================#
# Number of exons per gene ----
#============================================================================#
cat("======================================\n")
cat("Number of exons per gene\n")

df_rho = readRDS(file = "Data/Recombination/Gradient/gff_rho_all.rds")

species = list_dataset
for (i in 1:length(species)) {
  cat(species[i], "\n")
  genus = strsplit(species[i], split = "_")[[1]][1]
  sp = strsplit(species[i], split = "_")[[1]][2]
  species[i] = paste(genus, sp, sep = " ")
}
nb_exon = data.frame(set = rep(list_dataset, each = 20),
                     species = rep(species, each = 20),
                     nb_exons = 1:20,
                     nb_genes = NA)

for (i in 1:nrow(nb_exon)) {
  nb_exon$nb_genes[i] = sum(df_rho$nb_exons[which(df_rho$set == nb_exon$set[i])] == nb_exon$nb_exons[i], na.rm = TRUE)
}

saveRDS(nb_exon, file = "Data/Genomic_landscapes/nb_exons.rds")


#============================================================================#
# Relative Rho/kb ~ distance to hotspot center +- 5 kb ----
#============================================================================#
cat("======================================\n")
cat("Relative Rho/kb ~ distance to hotspot center\n")

# source("Jobs/rho_hotspot_centre.R")




#============================================================================#
# Rho gradient ----
#============================================================================#
cat("======================================\n")
cat("Rho gradient\n")

# source("Jobs/pool_rho_gradient_5kb.R")

# source("Jobs/pool_rho_gradient_25kb.R")


#============================================================================#
# Gradient per rank ----
#============================================================================#
cat("======================================\n")
cat("Gradient per rank")


rm(df_gradient_rank)
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
  
  gff_recombination$weighted.mean.rho.control = as.numeric(gff_recombination$weighted.mean.rho.control)
  
  
  # Subset
  gff_recombination = gff_recombination[which(gff_recombination$feature %in% c("CDS", "exon", "intron") & 
                                                gff_recombination$rank <= max.exons &
                                                gff_recombination$rank > 0),]
  
  # Random control - resampling ranks
  gff_recombination$rank_resample = NA
  set.seed(42)
  gff_recombination$rank_resample[which(gff_recombination$rank > 0 & gff_recombination$feature == "CDS")] = sample(gff_recombination$rank[which(gff_recombination$rank > 0) & gff_recombination$feature == "CDS"],
                                                                                                                   replace = TRUE)
  gff_recombination$rank_resample[which(gff_recombination$rank > 0 & gff_recombination$feature == "exon")] = sample(gff_recombination$rank[which(gff_recombination$rank > 0) & gff_recombination$feature == "exon"],
                                                                                                                    replace = TRUE)
  gff_recombination$rank_resample[which(gff_recombination$rank > 0 & gff_recombination$feature == "intron")] = sample(gff_recombination$rank[which(gff_recombination$rank > 0) & gff_recombination$feature == "intron"],
                                                                                                                      replace = TRUE)
  
  rho_rank = aggregate(weighted.mean.rho ~ rank, data = gff_recombination[which(gff_recombination$feature == "CDS" &
                                                                                  gff_recombination$rank <= max.exons &
                                                                                  gff_recombination$rank > 0),], mean)
  
  rho_rank.control = aggregate(weighted.mean.rho ~ rank_resample, data = gff_recombination[which(gff_recombination$feature == "CDS" &
                                                                                                   gff_recombination$rank <= max.exons &
                                                                                                   gff_recombination$rank > 0),], mean)
  
  rho_rank_exons = aggregate(weighted.mean.rho ~ rank, data = gff_recombination[which(gff_recombination$feature == "exon" &
                                                                                        gff_recombination$rank <= max.exons &
                                                                                        gff_recombination$rank > 0),], mean)
  
  rho_rank_exons.control = aggregate(weighted.mean.rho ~ rank_resample, data = gff_recombination[which(gff_recombination$feature == "exon" &
                                                                                                         gff_recombination$rank <= max.exons &
                                                                                                         gff_recombination$rank > 0),], mean)
  
  
  rho_rank_introns = aggregate(weighted.mean.rho ~ rank, data = gff_recombination[which(gff_recombination$feature == "intron" &
                                                                                          gff_recombination$rank <= max.exons &
                                                                                          gff_recombination$rank > 0),], mean)
  
  rho_rank_introns.control = aggregate(weighted.mean.rho ~ rank_resample, data = gff_recombination[which(gff_recombination$feature == "intron" &
                                                                                                           gff_recombination$rank <= max.exons &
                                                                                                           gff_recombination$rank > 0),], mean)
  colnames(rho_rank) = c("rank", "meanRho")
  colnames(rho_rank_exons) = c("rank", "meanRho")
  colnames(rho_rank_introns) = c("rank", "meanRho")
  colnames(rho_rank.control) = c("rank", "meanRho")
  colnames(rho_rank_exons.control) = c("rank", "meanRho")
  colnames(rho_rank_introns.control) = c("rank", "meanRho")
  
  rho_rank$condition = "observed"
  rho_rank_exons$condition = "observed"
  rho_rank_introns$condition = "observed"
  rho_rank.control$condition = "control"
  rho_rank_exons.control$condition = "control"
  rho_rank_introns.control$condition = "control"
  
  
  rho_rank = rbind(rho_rank, rho_rank.control)
  rho_rank_exons = rbind(rho_rank_exons, rho_rank_exons.control) 
  rho_rank_introns = rbind(rho_rank_introns, rho_rank_introns.control)
  
  rho_rank$level = "CDS" 
  rho_rank_exons$level = "exon" 
  rho_rank_introns$level = "intron"
  
  rho_rank = rbind(rho_rank, rho_rank_exons, rho_rank_introns)
  
  rho_rank$set = s
  
  if (exists("df_gradient_rank")) {
    df_gradient_rank = rbind(df_gradient_rank, rho_rank)
    rm(rho_rank)
  } else {
    df_gradient_rank = rho_rank
    rm(rho_rank)
  }
}
# df_gradient_rank

saveRDS(df_gradient_rank, file = paste("Data/Recombination/Gradient/meanRho_rank.rds", sep = ""))
rm(df_gradient_rank)


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
  
  gff_recombination$weighted.mean.rho.control = as.numeric(gff_recombination$weighted.mean.rho.control)
  
  # Subset
  gff_recombination = gff_recombination[which(gff_recombination$feature %in% c("CDS", "exon", "intron") & 
                                                gff_recombination$rank <= max.exons &
                                                gff_recombination$rank > 0),]
  
  if (!(s %in% c("Brassica_napus_Wu2019"))) {
    rho_rank_exons = aggregate(weighted.mean.rho ~ rank, data = gff_recombination[which(gff_recombination$feature == "CDS" & gff_recombination$rank <= max.exons & gff_recombination$rank > 0 & gff_recombination$nb_exons <= max.exons),], mean)
    
    rho_rank_exons_perclass = aggregate(weighted.mean.rho ~ rank + nb_exons, data = gff_recombination[which(gff_recombination$feature == "CDS" & gff_recombination$rank <= max.exons & gff_recombination$rank > 0 & gff_recombination$nb_exons <= max.exons),], mean)
    
    rho_rank_exons_perclass = rho_rank_exons_perclass[which(rho_rank_exons_perclass$nb_exons < 15),]
    
    rho_rank_exons$nb_exons = "all"
    rho_rank_exons_perclass = rbind(rho_rank_exons_perclass,
                                    rho_rank_exons)
    
    rho_rank_exons_perclass$nb_exons = factor(rho_rank_exons_perclass$nb_exons, levels = c(1:15, "all"))
  } else {
    rho_rank_exons_perclass = aggregate(weighted.mean.rho ~ rank, data = gff_recombination[which(gff_recombination$feature == "CDS" & gff_recombination$rank <= max.exons & gff_recombination$rank > 0),], mean)
    rho_rank_exons_perclass$nb_exons = "all"
  }
  
  rho_rank_exons_perclass$set = s
  
  if (exists("df_gradient_rank_nbexons")) {
    df_gradient_rank_nbexons = rbind(df_gradient_rank_nbexons, rho_rank_exons_perclass)
    rm(rho_rank_exons_perclass)
  } else {
    df_gradient_rank_nbexons = rho_rank_exons_perclass
    rm(rho_rank_exons_perclass)
  }
}
# df_gradient_rank_nbexons

saveRDS(df_gradient_rank_nbexons, file = paste("Data/Recombination/Gradient/meanRho_rank_nbexons.rds", sep = ""))
rm(df_gradient_rank_nbexons)



#============================================================================#
# Hotspot gradient per rank ----
#============================================================================#
cat("======================================\n")
cat("Hotspot gradient per rank\n")

# source("Jobs/hotspot_gradient_rank.R")



#============================================================================#
# Mean Rho in genomic features ----
#============================================================================#
cat("======================================\n")
cat("Mean rho in genomic features\n")

# source("Jobs/meanrho_per_feature.R")



#============================================================================#
# Weighted Mean Rho per bp in genomic features ----
#============================================================================#
cat("======================================\n")
cat("Mean rho weighted per sequence length in genomic features\n")

# source("Jobs/weightedmeanrho_per_feature.R")

#============================================================================#
# Summary statistics ----
#============================================================================#


# Nb genes
# Gene length bp
# Gene length nb exons
# Mean exon length bp
# Mean intron length bp
# Mean utr5 length
# Mean utr3 length

# Mean Rho·
# sd Rho
# Mean Rho exon 1
# sd Rho exon1
# Spearman rho gene ~ length bp
# Spearman rho gene ~ nb of exons
# Spearman rho CDS ~ rank

# Nb of hotspots - raw
# Nb of hotspots - trimmed
# Hotspot mean size
# Hotspot size sd
# Hotspot intensity (min/max)


#============================================================================#
# End of script ----
#============================================================================#
 
