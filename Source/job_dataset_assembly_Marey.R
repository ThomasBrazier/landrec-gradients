#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


########################################################################## #
#     JOB - Dataset Assembly
########################################################################## #

#============================================================================#
# LOADING ENVIRONMENT ----
#============================================================================#
library(readODS)
library(ggplot2)
library(pbmcapply)
library(viridis)
library(ggpubr)
library(tidyr)
library(tidyverse)
library(gridExtra)
library(purrr)


# TODO keys to merge dataset (gff, LD maps) are
# Gene id, chromosome, positions

# TODO Weighted mean Rho in a given window (chrom, start, end)

# TODO Datasets
# - pooled LD maps
# - Pooled Hotspots
# - Species GFF with Rho
# - Species Genes only with Rho
# - pooled GFFs with Rho
# - pooled Genes only with Rho

# TODO GFF Features
# - introns in UTR
# - mono-exonic
# - gene overlap hotspot/hotspot position
# - Number of SNPs
# - Number of exons of the gene
# - ID of the gene

# TODO Summary statistics
# - nb SNPs, SNP density
# - number/proportion of genes overlapping SNPs
# - SNP overlapping TSS (± 100 bp) ?
# - SNP overlapping ATG (± 100 bp) ?



#============================================================================#
# Loading variables & objects ----
#============================================================================#
# Get the directory of the file & set working directory
# wd=dirname(rstudioapi::getSourceEditorContext()$path)
# wd=gsub("/Source", "", wd)
# setwd(wd)

ncpus = 8
ncores = min(ncpus, detectCores())


#============================================================================#
# Import data from cluster ----
#============================================================================#



#============================================================================#
# Marey dataset ----
#============================================================================#

# Pooled GFFs for Marey maps ----
# List of dataset
marey_data = read.table(paste("Data/Recombination/marey_dataset.csv", sep = ""), header = TRUE, sep = "\t")
chromosome_metadata = read.table(paste0("Data/Genome/genome_chromosome_metadata.csv"),
                                 header = TRUE, sep = "\t")
# List of dataset
set_list = marey_data$set
# Discard species without GC distribution
set_list = set_list[which(!(set_list %in% c("Nelumbo_nucifera_Gui2018", "Gossypium_hirsutum_Zhang2019", "Dioscorea_alata_Cormier2019")))]

# Pool all species
rm(gff_Marey_pooled)
for (i in 1:length(set_list)) {
  cat("\n", set_list[i], "\n")
  gfffile = paste("Data/Genomic_landscapes/GC_genes/", set_list[i], ".csv.gz", sep = "")
  gff = read.table(gzfile(gfffile), header = T, sep = "\t")
  # TODO Bugfix: Fin de fichier (EOF) dans une chaîne de caractères entre guillemnts le nombre d'objets lus n'est pas un multiple du nombre de colonne
  # Filter protein coding genes only
  if ("protein_coding" %in% gff$gene_biotype) {
    pcgenes = gff$id[which(gff$gene_biotype == "protein_coding")]
    gff = gff[which(gff$id %in% pcgenes | gff$parent %in% pcgenes | gff$parent %in%
                    gff$id[which(gff$parent %in% pcgenes)]),]
  }
  gff$set = set_list[i]
  # Filter only chromosomes
  gff = gff[which(gff$seqname %in% chromosome_metadata$annotname[which(chromosome_metadata$set == set_list[i])]),]
  gff = gff[, !(colnames(gff) %in% c("index", "gene_biotype", "gene_length", "exon_length", "intron_length", "gene_nbexons"))]
  if (exists("gff_Marey_pooled")) {
    gff_Marey_pooled = rbind(gff_Marey_pooled, gff)
  } else {
    gff_Marey_pooled = gff
  }
}
# Save multi-species GFF dataset (all species with a Marey map and annotation)
save(gff_Marey_pooled, file = paste0("Data/gff_Marey_pooled.Rda"))


# Pooled 100kb marey landscapes ----
for (i in 1:length(set_list)) {
  gc100kb = read.table(gzfile(paste("Data/Genomic_landscapes/GC_100kb/", set_list[i], ".csv.gz", sep = "")), header = TRUE)
  gc100kb$set = set_list[i]
  if (exists("Marey_landscapes100kb")) {
    Marey_landscapes100kb = rbind(Marey_landscapes100kb, gc100kb)
  } else {
    Marey_landscapes100kb = gc100kb
  }
}
save(Marey_landscapes100kb, file = paste0("Data/Marey_landscapes100kb.Rda"))

# Gradient of GC in Marey species ----
GC_gradients = data.frame(set = rep(unique(gff_Marey_pooled$set), each = 20))
GC_gradients$rank = rep(1:20, length(unique(gff_Marey_pooled$set)))

GC_gradients$gc = NA
GC_gradients$gc3 = NA
for (i in 1:nrow(df)) {
  GC_gradients$gc[i] = mean(gff_Marey_pooled$gc[which(gff_Marey_pooled$feature == "CDS" & gff_Marey_pooled$set == df$set[i] & gff_Marey_pooled$rank == df$rank[i])], na.rm = TRUE)
  GC_gradients$gc3[i] = mean(gff_Marey_pooled$gc3[which(gff_Marey_pooled$feature == "CDS" & gff_Marey_pooled$set == df$set[i] & gff_Marey_pooled$rank == df$rank[i])], na.rm = TRUE)
}

save(GC_gradients, file = paste0("Data/GC_gradients.Rda"))


#============================================================================#
# GFF dataset ----
#============================================================================#

genids = pbmclapply(1:nrow(gff_recombination), function(x){genid(x)})
length(unlist(genids))

gff_recombination$nb_exons = NA


#============================================================================#
# LD dataset ----
#============================================================================#

# Pooled LD maps ----
# List of LD dataset
set_list = c("Oryza_sativa_Wang2018", "Citrullus_lanatus_Guo2019", "Brassica_napus_Wu2019", "Phaseolus_vulgaris_Wu2020",
             "Solanum_lycopersicum_Gupta2020", "Malus_domestica_Sun2020", "Sorghum_bicolor_Lozano2021",
             "Camellia_sinensis_Zhang2021", "Populus_tremula_Liu2022",
             "Arabidopsis_thaliana_1001genomes", "Glycine_max_Yang2021")
set_list

bpen = 5
if (exists("LD_landscapes")) {
  rm(LD_landscapes)
}

for (i in 1:length(set_list)) {
  # cat("\n", set_list[i], "\n")
  chr_list = chromosome_metadata$ldmapname[which(chromosome_metadata$set == set_list[i])]
  # cat("Chromosomes: ", chr_list, "\n")
  # Load the map only if the chromosome has been processed
  for (c in chr_list) {
    if (file.exists(paste("Data/Recombination/LD/ldhat/", set_list[i], ".", c, ".bpen", bpen,".res.txt.gz", sep = ""))) {
      ldmap = read.table(gzfile(paste("Data/Recombination/LD/ldhat/", set_list[i], ".", c, ".bpen", bpen,".res.txt.gz", sep = "")),
                         header = TRUE)
      ldmap$set = set_list[i]
      ldmap$chromosome = c
      # -------------- LD map quality ----------------------
      # Put positions in bp
      ldmap$Loci = ldmap$Loci * 1000
      # Remove negative positions
      ldmap = ldmap[-which(ldmap$Loci < 0),]
      ldmap$start = c(0, head(ldmap$Loci, -1))
      ldmap$end = (ldmap$Loci - 1)
      if (exists("LD_landscapes")) {
        LD_landscapes = rbind(LD_landscapes, ldmap)
      } else {
        LD_landscapes = ldmap
      }
    }
  }
}
save(LD_landscapes, file = paste0("Data/LD_landscapes.Rda"))


# Pooled GFFs for LD maps + estimates (e.g. Rho, GC) ----
# Discard species without genic GC estimates
set_list = set_list[which(!(set_list %in% c("Brassica_napus_Wu2019")))]
# Pool all species
if (exists("gff_LD_pooled")) {
  rm(gff_LD_pooled)
}
for (i in 1:length(set_list)) {
  # cat("\n", set_list[i], "\n")
  gfffile = paste("Data/Genomic_landscapes/GC_genes/", set_list[i], ".csv.gz", sep = "")
  gff = read.table(gzfile(gfffile), header = T, sep = "\t")
  # TODO Bugfix: Fin de fichier (EOF) dans une chaîne de caractères entre guillemnts le nombre d'objets lus n'est pas un multiple du nombre de colonne
  # Filter protein coding genes only
  if ("protein_coding" %in% gff$gene_biotype) {
    pcgenes = gff$id[which(gff$gene_biotype == "protein_coding")]
    gff = gff[which(gff$id %in% pcgenes | gff$parent %in% pcgenes | gff$parent %in%
                      gff$id[which(gff$parent %in% pcgenes)]),]
  }
  gff$set = set_list[i]
  # Filter only chromosomes
  gff = gff[which(gff$seqname %in% chromosome_metadata$annotname[which(chromosome_metadata$set == set_list[i])]),]
  gff = gff[, !(colnames(gff) %in% c("index", "gene_biotype", "gene_length", "exon_length", "intron_length", "gene_nbexons"))]
  
  if (exists("gff_LD_pooled")) {
    gff_LD_pooled = rbind(gff_LD_pooled, gff)
  } else {
    gff_LD_pooled = gff
  }
}
# Save multi-species GFF dataset (all species with a Marey map and annotation)
save(gff_LD_pooled, file = paste0("Data/gff_LD_pooled.Rda"))



# Rho in genes ----
load(file = paste0("Data/gff_LD_pooled.Rda"))
Genic_Rho = gff_LD_pooled[which(gff_LD_pooled$feature == "gene"),]

Genic_Rho$rho = NA
for (i in 1:nrow(genes)) {
  chr = chromosome_metadata$ldmapname[which(chromosome_metadata$annotname == Genic_Rho$seqname[i] & chromosome_metadata$set == Genic_Rho$set[i])]
  Genic_Rho$rho[i] = mean(LD_landscapes$Mean_rho[which(LD_landscapes$set == Genic_Rho$set[i] &
                                                     LD_landscapes$chromosome == chr &
                                                     LD_landscapes$start >= Genic_Rho$start[i] &
                                                     LD_landscapes$end <= Genic_Rho$end[i])], na.rm = TRUE)
}
save(genes, file = paste0("Data/Genic_Rho.Rda"))




idgenes = unique(gff_recombination$id[which(gff_recombination$feature == "gene")])
# idgenes = idgenes[1:10]
# list_idgenes = pbmclapply(idgenes, function(x){which(grepl(gsub("gene-", "", x), gff_recombination$parent) | grepl(gsub("gene-", "", x), gff_recombination$id))})

gff_recombination$gene_id = NA

# pb = txtProgressBar(min = 1, max = length(idgenes))
# for (i in 1:length(idgenes)) {
#   gff_recombination$nb_exons[list_idgenes[[i]]] = idgenes[i]
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# # apply(1:length(idgenes), function(x){gff_recombination$nb_exons[list_idgenes[[x]]] = idgenes[x]})
idgenes_modified = gsub("gene-", "", idgenes)

pb = txtProgressBar(min = 1, max = nrow(gff_recombination), style = 3)
for (i in 1:nrow(gff_recombination)) {
  if (gff_recombination$feature[i] == "gene") {
    gff_recombination$gene_id[i] = gff_recombination$id[i]
  } else {
    gff_recombination$gene_id[i] = idgenes[which(unlist(lapply(idgenes_modified, function(x){grepl(x, gff_recombination$parent[i])})))]
  }
  setTxtProgressBar(pb, i)
}
close(pb)

genid = function(x) {
  if (gff_recombination$feature[i] == "gene") {
    gene_id = gff_recombination$id[i]
  } else {
    gene_id = idgenes[which(unlist(lapply(idgenes_modified, function(x){grepl(x, gff_recombination$parent[i])})))]
  }
  return(gene_id)
}





#============================================================================#
# Summary statistics ----
#============================================================================#





#============================================================================#
# End of script ----
#============================================================================#
