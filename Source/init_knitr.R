marey_data = read.table("Data/Recombination/marey_dataset.csv",
                        header = TRUE,
                        sep = "\t")
chromosome_metadata = read.table("Data/Genome/genome_chromosome_metadata.csv",
                                 header = TRUE,
                                 sep = "\t")

# Loading libraries and scripts
packages = c("rstudioapi","tidyr","ggplot2","MareyMap","stringr","gdata","tidyverse","taxizedb","dplyr","rentrez","rBLAST","Biostrings","purrr","ggpubr","reshape2","RColorBrewer","ggplotify","pals","cowplot","lme4","sjPlot","sjmisc","phyr","kableExtra","car","ape","plyr","lmerTest", "viridis", "pbmcapply", "MuMIn", "ggtree", "treeio", "ggnewscale", "dplyr", "data.table", "xlsx", "gridExtra", "hrbrthemes", "parallel", "readr","rcartocolor", "ggh4x", "GenomicRanges", "yaml")
# lapply(packages, function(x) install.packages(as.character(x)))
lapply(packages, function(x) require(as.character(x), character.only = T))

source(paste0("Source/read.ldmap.R"))
source(paste0("Source/read.ldhot.R"))

list_dataset = c("Malus_sieversii_Sun2020", "Sorghum_bicolor_Lozano2021",
                 "Populus_tremula_Liu2022", "Glycine_max_Yang2021", "Homo_sapiens_Sudmant2015",
                 "Phaseolus_vulgaris_Wu2020", "Camellia_sinensis_Zhang2021", "Citrullus_lanatus_Guo2019",
                 "Spinacia_oleracea_Cai2021", "Arabidopsis_thaliana_1001genomes",
                 "Oryza_sativa_Wang2018", "Triticum_aestivum_Zhou2020")


list_species_ordered = c("Arabidopsis thaliana", "Glycine max", "Populus tremula", "Homo sapiens",
                         "Camellia sinensis", "Citrullus lanatus", "Malus sieversii", "Oryza sativa",
                         "Phaseolus vulgaris", "Sorghum bicolor", "Spinacia oleracea", "Triticum aestivum")
list_species = c("Arabidopsis thaliana", "Glycine max", "Populus tremula")

list_species_short = c("A. thaliana", "G. max", "P. tremula", "H. sapiens",
                         "C. sinensis", "C. lanatus", "M. sieversii", "O. sativa",
                         "P. vulgaris", "S. bicolor", "S. oleracea", "T. aestivum")

list_species_abbreviations = c("A. thal.", "G. max.", "P. trem.", "H. sap.",
                         "C. sin.", "C. lan.", "M. sie.", "O. sat.",
                         "P. vul.", "S. bic.", "S. ole.", "T. aest.")


params = read_yaml("parameters.yaml")
bpen = params$bpen
max.exons = params$max.exons
hotspot.maxlength = params$hotspot.maxlength
peakrate = params$peakrate

dataset_metadata = readODS::read_ods("Data/Metadata/dataset.ods")


# data_all = readRDS(file = paste0("Data/Recombination/Gradient/gff_rho_all.rds"))
# data_all = data_all[which(data_all$set %in% list_dataset),]
# # data_all = data_all[which(data_all$feature %in% c("gene", "mRNA", "exon", "intron", "CDS", "utr3", "utr5", "promoter", "intergenic", "upstream1kb", "upstream2kb", "upstream3kb", "downstream1kb", "downstream2kb", "downstream3kb")),]
# data_all = data_all[which(data_all$feature %in% c("gene", "exon", "intron", "CDS", "utr3", "utr5", "intergenic", "upstream1kb", "upstream2kb", "upstream3kb", "downstream1kb", "downstream2kb", "downstream3kb")),]
# data_all$feature = factor(data_all$feature, levels = c("intergenic", "gene", "upstream1kb", "upstream2kb", "upstream3kb", "utr5", "exon", "CDS", "intron", "utr3", "downstream1kb", "downstream2kb", "downstream3kb"))
if (file.exists(paste0("Data/Recombination/Gradient/gff_rho_all.rds"))) {
  data_all = readRDS(file = paste0("Data/Recombination/Gradient/gff_rho_all.rds"))
  data_all = data_all[which(data_all$set %in% list_dataset),]
  # data_all = data_all[which(data_all$feature %in% c("gene", "mRNA", "exon", "intron", "CDS", "utr3", "utr5", "promoter", "intergenic", "upstream1kb", "upstream2kb", "upstream3kb", "downstream1kb", "downstream2kb", "downstream3kb")),]
  data_all = data_all[which(data_all$feature %in% c("gene", "exon", "intron", "CDS", "utr3", "utr5", "intergenic", "intergenic_buffer", "buffer_upstream", "upstream1kb", "upstream2kb", "upstream3kb", "downstream1kb", "downstream2kb", "downstream3kb", "buffer_downstream")),]
  data_all$feature = factor(data_all$feature, levels = c("intergenic", "intergenic_buffer", "gene", "buffer_upstream", "upstream1kb", "upstream2kb", "upstream3kb", "utr5", "exon", "CDS", "intron", "utr3", "downstream1kb", "downstream2kb", "downstream3kb", "buffer_downstream"))
}

# Remove discarded chromosomes based on LDhot results
ldhot = read.table(file = gzfile("Data/Recombination/LD/ldhotspots_raw.csv.gz"), header = T)
ldhot$dataset = as.factor(ldhot$dataset)
ldhot = hotspot_filter(ldhot, 10^8, 0, 0, 10^6)

# Remove discarded chromosomes in GFF
data_all = data_all[(data_all$chromosome %in% ldhot$seqnames & data_all$set %in% ldhot$dataset),]  

data_all$species = factor(data_all$species,
                          levels = list_species_ordered)

gc()
