#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

########################################################################## #
#     JOB - GFF with population-scaled recombination rates
########################################################################## #
# test if there is at least one argument: if not, return an error
if (length(args) != 1) {
  stop("One argument must be supplied. <dataset>.", call.=FALSE)
} else if (length(args) == 1 ) {
  dataset = args[1]
}

#============================================================================#
# LOADING ENVIRONMENT ----
#============================================================================#
library(readODS)
library(pbmcapply)
library(tidyr)
library(tidyverse)
library(data.table)
library(GenomicRanges)

source('Source/read.ldmap.R')
source('Source/read.ldhot.R')
source('Source/read.gff.R')

ncpus = 8
ncores = min(ncpus, detectCores())

cat(dataset, "\n")

# List of dataset
marey_data = read.table(paste("Data/Recombination/marey_dataset.csv", sep = ""), header = TRUE, sep = "\t")
chromosome_metadata = read.table(paste0("Data/Genome/genome_chromosome_metadata.csv"),
                                 header = TRUE, sep = "\t")


cat("==============================================\n")
cat("Import gff...\n")
gfffile = paste("Data/Genomic_landscapes/GFF_parsed/", dataset, "_parsed.csv.gz", sep = "")

# Make the dataset - Pooled chromosomes
# gff = read.table(gzfile(gfffile), header = T, sep = "\t", quote = "",
#                  na.strings = c(NA, "NaN", ""), fill = TRUE)
gff = read.gff(gfffile)


# Filter protein coding genes only
# Do not filter when no protein coding information - all genes considered as protein coding
if ("protein_coding" %in% gff$gene_biotype) {
  pcgenes = gff$id[which(gff$gene_biotype == "protein_coding")]
  gff = gff[which(gff$id %in% pcgenes | gff$parent %in% pcgenes | gff$parent %in% gff$id[which(gff$parent %in% pcgenes)]),]
}

cat("==============================================\n")
cat("Make GFF genomic ranges\n")
summary(gff)

# Remove regions with start < end
# Prevent the following error:
# Error in .width_as_unnamed_integer(width, msg = "an end that is greater or equal to its start minus one") : 
# each range must have an end that is greater or equal to its start minus one
gff = gff[which((gff$start - 1) < gff$end),]

gff_ranges = makeGRangesFromDataFrame(gff,
                                      keep.extra.columns = TRUE,
                                      seqnames.field = "seqname",
                                      start.field = "start",
                                      end.field = "end")

gff_ranges

#============================================================================#
# 5. Add intergenic regions ----
cat("==============================================\n")
cat("Intergenic regions\n")

if ('gene' %in% unique(gff_ranges$feature)) {
  genes = gff_ranges[which(gff_ranges$feature == "gene")]
} else {
  if ('mRNA' %in% unique(gff_ranges$feature)) {
    genes = gff_ranges[which(gff_ranges$feature == "mRNA")]
  }
}

chr_sizes = chromosome_metadata$chrsize.bp[which(chromosome_metadata$set == dataset)]

chr_ranges = makeGRangesFromDataFrame(data.frame(seqnames = as.character(levels(seqnames(gff_ranges))),
                                                 start = rep(1, length(chr_sizes)),
                                                 end = as.numeric(chr_sizes), strand = "*"))

collapsed_genes = reduce(genes)
strand(collapsed_genes) = "*"
intergenic = setdiff(chr_ranges, collapsed_genes)
intergenic$source = NA
intergenic$feature = "intergenic"
intergenic$score = NA
intergenic$frame = NA
intergenic$attribute = NA
intergenic$id = paste0("intergenic-", seq(1, length(intergenic)))
intergenic$parent = NA
intergenic$name = NA
intergenic$gene_biotype = NA
intergenic$rank = NA

cat("==============================================\n")
intergenic


# Buffer +-3kb around genes
buffer = 3
intergenic_buffer = intergenic
start(intergenic_buffer) = start(intergenic_buffer) - buffer*10^3
end(intergenic_buffer) = end(intergenic_buffer) + buffer*10^3
intergenic_buffer$feature = "intergenic_buffer"

cat("==============================================\n")
intergenic_buffer


flanking_upstream = genes
start_pos = start(genes)
end_pos = end(genes)
start(flanking_upstream)[which(strand(flanking_upstream) == "+")] = start_pos[which(strand(flanking_upstream) == "+")] - buffer*10^3
start(flanking_upstream)[which(strand(flanking_upstream) == "-")] = end_pos[which(strand(flanking_upstream) == "-")] + 1
end(flanking_upstream)[which(strand(flanking_upstream) == "+")] = start_pos[which(strand(flanking_upstream) == "+")] - 1
end(flanking_upstream)[which(strand(flanking_upstream) == "-")] = end_pos[which(strand(flanking_upstream) == "-")] + buffer*10^3
flanking_upstream$feature = "flanking_upstream"

cat("==============================================\n")
flanking_upstream



flanking_downstream = genes
start_pos = start(genes)
end_pos = end(genes)
start(flanking_downstream)[which(strand(flanking_downstream) == "+")] = end_pos[which(strand(flanking_downstream) == "+")] + 1
start(flanking_downstream)[which(strand(flanking_downstream) == "-")] = start_pos[which(strand(flanking_downstream) == "-")] - buffer*10^3
end(flanking_downstream)[which(strand(flanking_downstream) == "+")] = end_pos[which(strand(flanking_downstream) == "+")] + buffer*10^3
end(flanking_downstream)[which(strand(flanking_downstream) == "-")] = start_pos[which(strand(flanking_downstream) == "-")] - 1
flanking_downstream$feature = "flanking_downstream"

cat("==============================================\n")
flanking_downstream

cat("==============================================\n")
cat("Concatenate results\n")
#gff_ranges = c(gff_ranges, intergenic, intergenic_buffer5kb, intergenic_buffer10kb)

gff_ranges = as.data.frame(gff_ranges)
intergenic = as.data.frame(intergenic)
intergenic_buffer = as.data.frame(intergenic_buffer)
flanking_upstream = as.data.frame(flanking_upstream)
flanking_downstream = as.data.frame(flanking_downstream)

cat("==============================================\n")
summary(gff_ranges)
cat("==============================================\n")
summary(intergenic)

df = rbind(gff_ranges, intergenic, intergenic_buffer, flanking_upstream, flanking_downstream)

colnames(df)[which(colnames(df) == "seqnames")] = "seqname"


cat("==============================================\n")
cat("Flag flanking overlapping a gene or another flanking (redundant)\n")
cat("==============================================\n")


gff_ranges = makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)

gff_genes = gff_ranges[which(gff_ranges$feature == "gene")]

isOverlappingGene = countOverlaps(gff_ranges, gff_genes, minoverlap = 2)
# table(isOverlappingGene)
# table(isOverlappingGene[which(gff_ranges$feature == "upstream3kb")])

isOverlappingGene = ifelse(isOverlappingGene > 0, TRUE, FALSE)
# table(isOverlappingGene)

gff_flanking = gff_ranges[which(gff_ranges$feature %in% c("upstream3kb",
                                                          "upstream2kb",
                                                          "upstream1kb",
                                                          "downstream1kb",
                                                          "downstream2kb",
                                                          "downstream3kb"))]

isOverlappingFlanking = countOverlaps(gff_ranges, gff_flanking, minoverlap = 2)
isOverlappingFlanking = isOverlappingFlanking - 1
# table(isOverlappingFlanking)
# table(isOverlappingFlanking[which(gff_ranges$feature == "upstream3kb")])

isOverlappingFlanking = ifelse(isOverlappingFlanking > 0, TRUE, FALSE)
# table(isOverlappingFlanking)


gff_ranges$isOverlappingGene = isOverlappingGene
gff_ranges$isOverlappingFlanking = isOverlappingFlanking

df = as.data.frame(gff_ranges)

colnames(df)[1] = "seqname"

cat("==============================================\n")
cat("Translate chromosome names\n")
cat("==============================================\n")
annotname = chromosome_metadata$annotname[chromosome_metadata$set == dataset]
ldmapname = chromosome_metadata$ldmapname[chromosome_metadata$set == dataset]

df$seqname = factor(df$seqname, levels = annotname, labels = ldmapname)

gfffile = paste("Data/Genomic_landscapes/GFF_parsed/", dataset, ".csv.gz", sep = "")

write.table(df, file = gzfile(gfffile),
            row.names = F, col.names = T, sep = "\t", quote = F)
