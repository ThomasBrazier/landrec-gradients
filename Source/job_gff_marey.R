#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


########################################################################## #
#     JOB - Infer gene ID and number of exons in Marey GFFs
########################################################################## #
# test if there is at least one argument: if not, return an error
if (length(args) != 1) {
  stop("One argument must be supplied: <dataset>.", call.=FALSE)
} else if (length(args) == 1 ) {
  set = args[1]
}

# Configure here your dataset manually
# set = "Arabidopsis_thaliana_1001genomes"


#============================================================================#
# LOADING ENVIRONMENT ----
#============================================================================#
library(readODS)
library(pbmcapply)
library(tidyr)
library(tidyverse)
library(data.table)
library(GenomicRanges)


#============================================================================#
# Loading variables & objects ----
#============================================================================#
# Get the directory of the file & set working directory
# wd=dirname(rstudioapi::getSourceEditorContext()$path)
# wd=gsub("/Source", "", wd)
# setwd(wd)


ncpus = 8
ncores = min(ncpus, detectCores())

cat(set, "\n")

# List of dataset
marey_data = read.table(paste("Data/Recombination/marey_dataset.csv", sep = ""), header = TRUE, sep = "\t")
chromosome_metadata = read.table(paste0("Data/Genome/genome_chromosome_metadata.csv"),
                                 header = TRUE, sep = "\t")

# List of chromosomes in the LD map
# map_list = chromosome_metadata$ldmapname[which(chromosome_metadata$set == set)]
# map_list = map_list[which(!is.na(map_list))]
# map_list = chromosome
# map_list
# # List of chromosomes in the gff
# chr_list = chromosome_metadata$annotname[which(chromosome_metadata$set == set &
#                                                  chromosome_metadata$ldmapname == chromosome)]
# chr_list = chr_list[which(!is.na(chr_list))]
# chr_list

#============================================================================#
# GFF dataset ----
#============================================================================#

cat("Import gff...\n")
gfffile = paste("Data/Genomic_landscapes/GC_genes/", set, ".csv.gz", sep = "")
# gff = read.table(gzfile(gfffile), header = T, sep = "\t")
# gfffile = paste(wd, "/Data/Genomic_landscapes/GFF_parsed/", set, "_parsed.csv.gz", sep = "")

# Make the dataset - Pooled chromosomes
gff = read.table(gzfile(gfffile), header = T, sep = "\t", quote = "", 
                 na.strings = c(NA, "NaN", ""), fill = TRUE)
# gff = gff[which(gff$seqname %in% chr_list),]
# Filter protein coding genes only
# Do not filter when no protein coding information - all genes considered as protein coding
if ("protein_coding" %in% gff$gene_biotype) {
  pcgenes = gff$id[which(gff$gene_biotype == "protein_coding")]
  gff = gff[which(gff$id %in% pcgenes | gff$parent %in% pcgenes | gff$parent %in% gff$id[which(gff$parent %in% pcgenes)]),]
}

# The parsed GFF where recombination rates has been added
# load(file = paste("Output/gff_", set, ".Rda", sep = ""))

# A subset for developping
# gff = gff[1:10000,]
cat("Make GFF genomic ranges\n")
gff = gff[which(gff$end > (gff$start - 1)),]
gff = gff[which(gff$strand %in% c("+", "-", "*")),]

gff_ranges = makeGRangesFromDataFrame(gff, keep.extra.columns=TRUE)



# 1. Get gene ID for each feature ----
cat("Gene IDs\n")

idgenes = unique(gff_ranges$id[which(gff_ranges$feature == "gene")])

# DEBUG
# idgenes = idgenes[1:100]


# Get children for each gene
# TODO Optimize, grepl is slow
children = function(id) {
  # Take a gene id in argument
  # Return all children of a gene and the gene itself
  
  # Remove prefixes
  # id_modified = gsub("gene-", "", id)
  # child = gff_ranges$id[which(grepl(id_modified, gff_ranges$parent))]
  # More simple method
  ch = gff_ranges$id[which(gff_ranges$parent == id)]
  ch = ch[!is.na(ch)] # Take care of NA values
  child = gff_ranges$id[which(gff_ranges$id == id | gff_ranges$parent == id | gff_ranges$parent %in% ch)]
  
  if (length(child) > 0) {
    child = data.frame(gene = id, children = child)
    return(child)
  } else {
    return(NULL)
  }
}
list_ids = pbmclapply(idgenes, function(x){children(x)})
list_ids = rbindlist(list_ids)
list_ids = list_ids[which(!grepl("rna", list_ids$gene)),]

# Get gene id for each feature
genid = function(x) {
  if (gff_ranges$feature[x] == "gene") {
    gene_id = gff_ranges$id[x]
  } else {
    gene_id = list_ids$gene[which(list_ids$children == gff_ranges$id[x])]
  }
  gene_id = unique(gene_id)
  if (length(gene_id) != 1) {
    gene_id = NA
  }
  return(gene_id)
}

# For C. lanatus, no need to infer gene id
if (set == "Citrullus_lanatus_Guo2019") {
  gff_ranges$gene_id = list_ids$gene
} else {
  genids = pbmclapply(1:length(gff_ranges), function(x){genid(x)})
  
  if (length(unlist(genids)) != length(gff_ranges)) {
    warning("Parsing gene ids introduced missing rows")
  }
  gff_ranges$gene_id = unlist(genids)
}


# 2. Get number of exons per gene ----
cat("Number of exons\n")

# idgenes = unique(gff_ranges$id[which(gff_ranges$feature == "gene")])

# Get number of exons for each gene
# !! Some annotations have only CDS
nbexons = function(id) {
  nbex = max(gff_ranges$rank[which(gff_ranges$gene_id == id & gff_ranges$feature == "CDS")])
  nbex = data.frame(gene = id, nb_exons = nbex)
  return(nbex)
}
nexons = pbmclapply(1:length(gff_ranges), function(x){nbexons(gff_ranges$gene_id[x])})
nexons = rbindlist(nexons)
# nexons = nexons[-grep("rna", nexons$gene),]
nexons$nb_exons[which(nexons$nb_exons == -Inf)] = NA

if (nrow(nexons) != length(gff_ranges)) {
  warning("Parsing number of exons introduced missing rows")
}
gff_ranges$nb_exons = nexons$nb_exons



# Save ----
cat("Save ", paste("Output/gc_genes_marey_", set, ".rds", sep = ""),"...\n")
saveRDS(gff_ranges, file = paste("Output/gc_genes_marey_", set, ".rds", sep = ""))

# df = as.data.frame(gff_ranges, row.names = NULL)
# 
# cat("Save .csv file\n")
# write.table(df, file = gzfile(paste("Output/gc_genes_marey_", set, ".csv.gz", sep = "")),
#             row.names = F, col.names = T, sep = "\t", quote = F)




#============================================================================#
# End of script ----
#============================================================================#
