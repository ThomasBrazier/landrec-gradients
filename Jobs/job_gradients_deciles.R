library(GenomicRanges)
library(ggplot2)
library(pbmcapply)
library(dplyr)
library(tidyr)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GenomicRanges")

list_dataset = c("Malus_sieversii_Sun2020", "Sorghum_bicolor_Lozano2021",
                 "Populus_tremula_Liu2022", "Glycine_max_Yang2021",
                 "Phaseolus_vulgaris_Wu2020", "Camellia_sinensis_Zhang2021", "Citrullus_lanatus_Guo2019",
                 "Spinacia_oleracea_Cai2021", "Arabidopsis_thaliana_1001genomes",
                 "Oryza_sativa_Wang2018", "Triticum_aestivum_Zhou2020", "Homo_sapiens_Sudmant2015")

ncores = 8
# ncores = 1

allmaps = read.table(gzfile("Data/Recombination/LD/ldhat/LD_maps.csv.gz"),
           header = TRUE, sep = "\t")

# rm(allmaps)
# ncores = 1
# ds = "Arabidopsis_thaliana_1001genomes"
rm(deciles)
n_bins = 10


for (ds in list_dataset) {
  cat(ds, "\n")
  ldmap = allmaps[which(allmaps$set == ds),]
  ldmap_ranges = makeGRangesFromDataFrame(ldmap, keep.extra.columns = TRUE)
  
  ldmap_ranges
  
  gff = readRDS(paste0("Data/Genomic_landscapes/Rho/gff_rho_", ds, ".rds"))
  
  columns = c("chromosome", "start", "end",                           
              "width", "strand", "source",                        
              "feature", "score", "frame",                         
              "attribute", "id", "parent",                        
              "name", "gene_biotype", "rank", "nb_exons", 
              "mean.rho", "median.rho", "mean.rho.control", "median.rho.control",
              "weighted.mean.rho", "weighted.mean.rho.control", "mean.rho.startafter",           
              "weighted.mean.rho.startafter", "relative.meanrho", "relative.wmeanrho")
  
  genes = gff[which(gff$feature == "gene" & gff$nb_exons <= 14),]
  
  genes = genes[,columns]
  
  # replicate each row/gene to make interval around
  gene_range = genes[rep(seq_len(nrow(genes)), each = n_bins), ]
  
  # Add index of physical distance around TSS
  # Consider the strand + or -
  gene_range$index = NA
  gene_range$index = rep(1:n_bins, nrow(gene_range)/10)
  
  gene_range_forward = gene_range[which(gene_range$strand == "+"),]
  gene_range_forward$start.gene = gene_range_forward$start
  gene_range_forward$end.gene = gene_range_forward$end
  gene_range_forward$start = gene_range_forward$start.gene + (gene_range_forward$index - 1) * round(gene_range_forward$width / n_bins, digits = 0)
  gene_range_forward$end = gene_range_forward$start.gene + (gene_range_forward$index) * round(gene_range_forward$width / n_bins, digits = 0)
  
  gene_range_reverse = gene_range[which(gene_range$strand == "-"),]
  gene_range_reverse$start.gene = gene_range_reverse$start
  gene_range_reverse$end.gene = gene_range_reverse$end
  gene_range_reverse$start = gene_range_reverse$end.gene - (gene_range_reverse$index) * round(gene_range_reverse$width / n_bins, digits = 0)
  gene_range_reverse$end = gene_range_reverse$end.gene - (gene_range_reverse$index - 1) * round(gene_range_reverse$width / n_bins, digits = 0)
  
  gene_range = rbind(gene_range_forward, gene_range_reverse)
  
  # gene_range$end = ifelse(gene_range$index == 10, gene_range$end.gene, gene_range$end)
  
  gene_range = gene_range[(gene_range$end > (gene_range$start - 1)),]
  
  gene_range = makeGRangesFromDataFrame(gene_range, keep.extra.columns = TRUE)
  gene_range
  
  
  # Estimate mean rho for each bin
  cat("Mean Rho\n")
  hits = findOverlaps(gene_range, ldmap_ranges)
  query = queryHits(hits)
  subj = subjectHits(hits)
  mean.rho = unlist(pbmclapply(1:length(gene_range), function(x) {if (x %in% query) {if (length(subj[which(query == x)]) > 0) {mean(ldmap_ranges$Mean_rho[subj[which(query == x)]], na.rm = TRUE)} else {NA}} else {NA}}, mc.cores = ncores))
  gene_range$mean.rho = mean.rho
  
  
  # hits = findOverlaps(gene_range, ldmap_ranges, type = c("within"))
  # query = queryHits(hits)
  # subj = subjectHits(hits)
  # mean.rho = unlist(pbmclapply(1:length(gene_range), function(x) {if (x %in% query) {if (length(subj[which(query == x)]) > 0) {mean(ldmap_ranges$Mean_rho[subj[which(query == x)]], na.rm = TRUE)} else {NA}} else {NA}}, mc.cores = ncores))
  # gene_range$mean.rho.within = mean.rho.within
  
  
  gene_range$dataset = ds
  
  if (!exists("deciles")) {
    deciles = as.data.frame(gene_range)
  } else {
    deciles = rbind(deciles, as.data.frame(gene_range))
  }
  
}


saveRDS(deciles, "Data/Recombination/Gradient/gradients_deciles.rds")
