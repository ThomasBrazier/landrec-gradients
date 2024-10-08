#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

########################################################################## #
#     JOB - LDhat Rho gradients
########################################################################## #

#============================================================================#
# ARGUMENTS ----
#============================================================================#
# test if there is at least one argument: if not, return an error
if (length(args) != 2) {
 stop("Two arguments must be supplied: <dataset> and <chromosome>.\n", call.=FALSE)
} else if (length(args) == 2) {
 dataset = args[1]
 chromosome = args[2]
}
# dataset = args[1]
#============================================================================#
# LOADING ENVIRONMENT ----
#============================================================================#
library(readODS)
library(pbmcapply)
library(tidyr)
library(tidyverse)
library(GenomicRanges)

source("Source/init.R")

ncores = 16
# ncores = min(ncores, detectCores())

#============================================================================#
# Get the list of chromosomes to process ----
#============================================================================#
cat("Get the list of chromosomes to process\n")
chromosome_metadata = read.table("Data/Genome/genome_chromosome_metadata.csv",
                                 header = TRUE, sep = "\t")

# List of chromosomes in the LD map
# map_list = chromosome_metadata$ldmapname[which(chromosome_metadata$set == dataset)]
# chr_list = chromosome_metadata$annotname[which(chromosome_metadata$set == dataset)]
if (chromosome == "all") {
 map_list = chromosome_metadata$ldmapname[which(chromosome_metadata$set == dataset)]
 chr_list = chromosome_metadata$annotname[which(chromosome_metadata$set == dataset)]
} else {
 map_list = chromosome_metadata$ldmapname[which(chromosome_metadata$set == dataset
                                                & chromosome_metadata$ldmapname == chromosome)]
 chr_list = chromosome_metadata$annotname[which(chromosome_metadata$set == dataset
                                                & chromosome_metadata$ldmapname == chromosome)]
}
map_list = map_list[which(!is.na(map_list))]
map_list
# List of chromosomes in the gff
chr_list = chr_list[which(!is.na(chr_list))]
chr_list
gff_list = chr_list

#============================================================================#
# Make the dataset - Pooled chromosomes ----
#============================================================================#

cat("Make the dataset - Pooled chromosomes\n")

cat("Import gff...\n")
gfffile = paste("Data/Genomic_landscapes/GFF_parsed/", dataset, "_parsed.csv.gz", sep = "")
# gff = read.table(gzfile(gfffile), header = T, sep = "\t")

# Make the dataset - Pooled chromosomes
gff = read.table(gzfile(gfffile), header = T, sep = "\t", quote = "", 
                 na.strings = c(NA, "NaN", ""))
# unique(gff$seqname)
gff = gff[which(gff$seqname %in% chr_list),]
# Filter protein coding genes only
# Do not filter when no protein coding information - all genes considered as protein coding
cat("Keep only protein-coding genes\n")
if ("protein_coding" %in% gff$gene_biotype) {
  pcgenes = gff$id[which(gff$gene_biotype == "protein_coding")]
  gff = gff[which(gff$id %in% pcgenes | gff$parent %in% pcgenes | gff$parent %in% gff$id[which(gff$parent %in% pcgenes)]),]
}
# Take care of empty strings -> convert to NA
gff$id[which(gff$id == "")] = NA
gff$parent[which(gff$parent == "")] = NA


# Pool LD maps
source("Source/read.ldmap.R")

ldmap_trimmed = read.ldmap(dataset, chr = map_list)
# Trimming, if necessary


#============================================================================#
# Mask non variable regions ----
#============================================================================#

# # Detect groups
# grp = seq_consecutive(ldmap_trimmed$Mean_rho)
# ldmap_trimmed$group = grp
# 
# # Get the genomic size of each group
# grp_start = aggregate(start ~ grp, ldmap_trimmed, min)
# grp_end = aggregate(end ~ grp, ldmap_trimmed, max)
# grp_rho = aggregate(Mean_rho ~ grp, ldmap_trimmed, mean)
# grp_sdrho = aggregate(Mean_rho ~ grp, ldmap_trimmed, sd)
# chr = aggregate(chromosome ~ grp, ldmap_trimmed, unique)                    
# 
# grp_length = grp_end - grp_start
# 
# g = data.frame(grp = unique(grp))
# g$grp_length = grp_length$end
# g$start = grp_start$start
# g$end = grp_end$end
# g$Mean_rho = grp_rho$Mean_rho
# g$SDMean_rho = grp_sdrho$Mean_rho
# g$chromosome = chr$chromosome
# 
# 
# msk_kb = 20
# mask = ldmap_trimmed$group %in% g$grp[which(g$grp_length > (msk_kb * 10^3))]
# 
# ldmap_trimmed$mask = mask
# ldmap_trimmed$Mean_rho[mask] = NA
# 


ldmap_ranges = ldmap_trimmed
# Change chromosome names to GFF annotation names
for (chr in unique(ldmap_ranges$chromosome)) {
  ldmap_ranges$chromosome[which(ldmap_ranges$chromosome == chr)] = chromosome_metadata$annotname[which(chromosome_metadata$ldmapname == chr & chromosome_metadata$set == dataset)]
}
ldmap_ranges = makeGRangesFromDataFrame(ldmap_ranges, keep.extra.columns = TRUE)


#============================================================================#
# Import GFF ----
#============================================================================#

# Reduce GFF to seqnames in LD maps
gff = gff[which(gff$seqname %in% seqnames(ldmap_ranges)),]

# Error in .width_as_unnamed_integer(width, msg = "an end that is greater or equal to its start minus one") : 
#   each range must have an end that is greater or equal to its start minus
# one
# Calls: makeGRangesFromDataFrame ... .new_IRanges_from_start_end -> .width_as_unnamed_integer
if (sum(gff$end < (gff$start - 1)) > 0) {
  gff = gff[-which(gff$end < (gff$start - 1)),]
}

#============================================================================#
# Recombination gradient in exons (bp) ----
#============================================================================#


# Distance from ATG start codon (i.e. start of CDS rank 1)

cat("Recombination gradient in exons (bp)\n")
interval = 200
upstream = -5000 # Size of the upstream region
downstream = 5000 # Max size of the downstream region

# Get each gene ID
cat("Splitting by gene\n")
# Reduce analyses to sequences in LD maps
ids = as.factor(gff$id[which(gff$feature == "gene" & gff$seqname %in% as.character(levels(seqnames(ldmap_ranges))))])
ids = split(ids, as.factor(ids))

sample.gene = function(id) {
  # Sample a gene from a given ID
  children = gff$id[which(gff$parent == id)]
  children = children[!is.na(children)] # Take care of NA values
  g = gff[which((gff$id == id | gff$parent == id | gff$parent %in% children)),]
  return(g)
}
# sample.gene(as.character(ids[[1]]))
# DEBUG
# ids = ids[1:100]
cat("....... Sample genes\n")
genes = pbmclapply(ids,
                   function(x){sample.gene(x)},
                   mc.cores = ncores)
# genes = lapply(1:length(ids), function(x){sample.gene(as.character(ids[[x]]))})
# DONE Bugfix not all parts of the gene are sampled

# Keep only genes with CDS
genes = genes[unlist(lapply(genes, function(x) {"CDS" %in% x$feature}))]

# For each gene, get a range and a random control ----

cat("Get position of TSS\n")
# Get position of first ATG ----
# TODO BUGFIX problem of wrong direction for strand "-"
# TODO take end instead of start for ATG
gene_tss = pbmclapply(genes,
                      function(x){if (sum(x$feature == "CDS" & x$rank == 1) > 0) {ifelse("+" %in% x$strand,
                                                                                                      min(x$start[which(x$feature == "gene" | x$feature == "mRNA")]),
                                                                                                      max(x$end[which(x$feature == "gene" | x$feature == "mRNA")]))}
  else {NA}}, mc.cores = ncores)

strand = pbmclapply(genes,
                    function(x){unique(as.character(x$strand))},
                    mc.cores = ncores)
seqname = lapply(genes,
                 function(x){unique(x$seqname)})

cat("Get distance to ATG\n")
# Get distance to ATG ----
atg = pbmclapply(genes,
                 function(x){if (sum(x$feature == "CDS" & x$rank == 1) > 0) {ifelse("+" %in% x$strand,
                                                                                     min(x$start[which(x$feature == "CDS" & x$rank == 1)]),
                                                                                     max(x$end[which(x$feature == "CDS" & x$rank == 1)]))}
  else {NA}}, mc.cores = ncores)

dist_atg = abs(unlist(gene_tss) - unlist(atg))


cat("Keep only first part of the gene\n")
# p.sample = 0.8 # Gene proportion to sample
# in order to avoid effects of gene end

# Get size of each gene ----
gene_size = unlist(pbmclapply(genes,
                              function(x){x$end[x$feature == "gene"] - x$start[x$feature == "gene"] + 1},
                              mc.cores = ncores))
# Calculate size to cut
# gene_size = gene_size * p.sample
# gene_size = unlist(pbmclapply(gene_size, function(x) {min(downstream, x)}, mc.cores = ncores))


cat("Make ranges\n")
# Make ranges ----
idx = seq(from = upstream, to = downstream, by = abs(interval))

# One specific range per gene
# idx = seq(from = upstream, to = downstream, by = abs(interval))
# idx = pbmclapply(gene_size, function(x) {seq(from = upstream, to = x, by = abs(interval))}, mc.cores = ncores)

make.range = function(x) {
  # Make a -upstream:+downstream interval in which to sample windows
  # Tak an index of an atg as argument
  if ("+" == unique(as.character(strand[x]))) {
    start = as.numeric(gene_tss[x]) + as.numeric(idx)
    range = GRanges(seqnames = as.character(seqname[x]),
                    ranges = IRanges(start = start,
                                     width = abs(interval)),
                    strand = as.character(strand[x]),
                    idx = as.numeric(idx),
                    gene_id = names(genes)[x],
                    gene_size = as.numeric(gene_size[x]),
                    gene_tss = as.character(gene_tss[x]),
                    dist_atg = dist_atg[x])
  } else {
    # TODO take end instead of start for ATG
    start = as.numeric(gene_tss[x]) - abs(interval) + 1 - as.numeric(idx)
    end = as.numeric(gene_tss[x]) - as.numeric(idx)
    range = GRanges(seqnames = as.character(seqname[x]),
                    ranges = IRanges(start = start, end = end),
                    strand = as.character(strand[x]),
                    idx = as.numeric(idx),
                    gene_id = names(genes)[x],
                    gene_size = as.numeric(gene_size[x]),
                    gene_tss = as.character(gene_tss[x]),
                    dist_atg = dist_atg[x])
  }
  return(range)
}
# View(as.data.frame(make.range(2)))
# gene_ranges = lapply(1:length(gene_tss), function(x){make.range(x)})
gene_ranges = pbmclapply(1:length(gene_tss),
                         function(x){make.range(x)},
                         mc.cores = ncores)

cat("Random ranges\n")
# Random ranges ----
seqlengths = pbmclapply(1:length(genes),
                        function(x){chromosome_metadata$chrsize.bp[which(chromosome_metadata$set == dataset & chromosome_metadata$annotname == unique(genes[[x]]$seqname))]},
                        mc.cores = ncores)
random_pos = pbmclapply(seqlengths,
                        function(x){sample(5001:(x - 5001), size = 1, replace = TRUE)},
                        mc.cores = ncores)

random.range = function(x) {
  if ("+" == unique(as.character(strand[x]))) {
    start = as.numeric(random_pos[x]) + as.numeric(idx)
    range = GRanges(seqnames = as.character(seqname[x]),
                    ranges = IRanges(start = start,
                                     width = abs(interval)), strand = as.character(strand[x]))
  } else {
    start = as.numeric(random_pos[x]) - abs(interval) + 1 - as.numeric(idx)
    end = as.numeric(random_pos[x]) - as.numeric(idx)
    range = GRanges(seqnames = as.character(seqname[x]),
                    ranges = IRanges(start = start, end = end),
                    strand = as.character(strand[x]))
  }
  return(range)
}

random_ranges = pbmclapply(1:length(random_pos),
                           function(x){random.range(x)},
                           mc.cores = ncores)

# pbmclapply(1:length(random_pos), function(x){if ("+" == unique(as.character(strand[x]))) {
#   GRanges(seqnames = as.character(seqname[x]),
#           ranges = IRanges(start = as.numeric(random_pos[x]) + idx_forward,
#                            width = abs(interval)), strand = as.character(strand[x]))
# } else {
#   GRanges(seqnames = as.character(seqname[x]),
#           ranges = IRanges(start = as.numeric(random_pos[x]) - (abs(interval) - 1) + idx_backward,
#                            width = abs(interval)), strand = as.character(strand[x]))
# }
# }, mc.cores = ncores)


# Go from gene-level to interval-level
gene_ranges = GRangesList(gene_ranges)
gene_ranges = unlist(gene_ranges)

random_ranges = GRangesList(random_ranges)
random_ranges = unlist(random_ranges)

cat("Estimate Rho\n")
# Get Rho for each interval ----
# Estimating mean Rho ----
index = seq(1, length(gene_ranges))
res = split(gene_ranges, as.factor(index))
hits = findOverlaps(res, ldmap_ranges)
hits2 = split(hits, as.factor(queryHits(hits)))
mean.rho = unlist(lapply(1:length(res),
                         function(x) {if (x %in% names(hits2)) {if (length(subjectHits(hits2[[which(names(hits2) == x)]])) > 0) {mean(ldmap_ranges$Mean_rho[subjectHits(hits2[[which(names(hits2) == x)]])], na.rm = TRUE)} else {NA}} else {NA}}))

index = seq(1, length(random_ranges))
res = split(random_ranges, as.factor(index))
hits = findOverlaps(res, ldmap_ranges)
hits2 = split(hits, as.factor(queryHits(hits)))
mean.rho.control = unlist(lapply(1:length(res), function(x) {if (x %in% names(hits2)) {if (length(subjectHits(hits2[[which(names(hits2) == x)]])) > 0) {mean(ldmap_ranges$Mean_rho[subjectHits(hits2[[which(names(hits2) == x)]])], na.rm = TRUE)} else {NA}} else {NA}}))

gene_ranges$mean.rho = mean.rho
gene_ranges$mean.rho.control = mean.rho.control


# Weighted mean Rho overlap
wmeanrho = function(x) {
  # x is a GRanges object of size 1
  hits = subjectHits(findOverlaps(x, ldmap_ranges))
  rho = ldmap_ranges$Mean_rho[hits]
  w = width(restrict(ldmap_ranges[hits], start(x), end(x)))
  wmean = weighted.mean(rho, w)
  return(wmean)
}

cat("Weighted mean rho\n")
# Mean Rho
# x = gene_ranges[1]
wmean = pbmclapply(1:length(gene_ranges),
                   function(x) {wmeanrho(gene_ranges[x])},
                   mc.cores = ncores)
gene_ranges$wmean.rho = unlist(wmean)

wmean.control = pbmclapply(1:length(gene_ranges),
                           function(x) {wmeanrho(random_ranges[x])},
                           mc.cores = ncores)
gene_ranges$wmean.rho.control = unlist(wmean.control)




cat("Check genes upstream\n")
# Check if the position upstream overlap another gene that the one sampled ----
# Discard if overlap
gene_boundaries = makeGRangesFromDataFrame(gff[which(gff$feature == 'gene'),],
                                           keep.extra.columns = T)
hits = findOverlaps(gene_ranges, gene_boundaries)
gene_ranges$gene_overlap = NA
gene_ranges$gene_overlap[queryHits(hits)] = TRUE
gene_ranges$gene_overlap[-queryHits(hits)] = FALSE

gene_ranges$gene_upstream_overlap = (gene_ranges$gene_overlap & (gene_ranges$idx < 0))

hits = countOverlaps(gene_ranges, gene_boundaries)
gene_ranges$nb_gene = hits

# DEPRECATED AFTER OPTIMIZATION
# list_ranges = split(gene_ranges, as.factor(1:length(gene_ranges)))
# TRUE if upstream position (i.e. negative idx) overlap a gene before
# gene_boundaries = makeGRangesFromDataFrame(gff[which(gff$feature == 'gene'),],
#                                            keep.extra.columns = T)
# gene_boundaries$gene_id = gff$id[which(gff$feature == 'gene')]
# gene_overlap = pbmclapply(list_ranges, function(x) {ifelse(x$idx < 0, countOverlaps(x, gene_boundaries[-which(gene_boundaries$gene_id == x$gene_id)]) > 0, NA)},
#                           mc.cores = ncores)
# gene_ranges$gene_overlap = unlist(gene_overlap)
# # Get number of genes on the position ----
# nb_gene = pbmclapply(list_ranges, function(x) {countOverlaps(x, gene_boundaries)},
#                           mc.cores = ncores)
# gene_ranges$nb_gene = unlist(nb_gene)


cat("Check exon/intron overlap\n")
# Check if I am in exons or introns for each pos ----
# Overlapping at least p % of the interval
p_overlap = 0.7
exon_boundaries = makeGRangesFromDataFrame(gff[which(gff$feature == 'CDS'),])
exon_overlap = (countOverlaps(gene_ranges, exon_boundaries, minoverlap = p_overlap*interval) > 0)
gene_ranges$exon_overlap = exon_overlap

intron_boundaries = makeGRangesFromDataFrame(gff[which(gff$feature == 'intron'),])
intron_overlap = (countOverlaps(gene_ranges, intron_boundaries, minoverlap = p_overlap*interval) > 0)
gene_ranges$intron_overlap = intron_overlap

# 
# 
# # How much it overlap
# size_overlap = function(x, y) {
#   # x is a GRanges object of size 1
#   # y is a GRanges object in which to search overlaps, e.g. introns or exons
#   overlap = findOverlaps(x, y)
#   if (length(overlap) > 0) {
#     query = queryHits(overlap)
#     subject = subjectHits(overlap)
#     wi = (intersect(x, y[subject]))
#     wi = sum(width(wi), na.rm = TRUE)
#   } else {
#     wi = NA
#   }
#   return(wi)
# }
# cat("Exon overlap\n")
# exon_boundaries = makeGRangesFromDataFrame(gff[which(gff$feature == 'CDS'),])
# ex_over = pbmclapply(1:length(gene_ranges), function(x) {size_overlap(gene_ranges[x], exon_boundaries)}, mc.cores = ncores)
# gene_ranges$exon_overlap = unlist(ex_over)
# 
# cat("Intron overlap\n")
# intron_boundaries = makeGRangesFromDataFrame(gff[which(gff$feature == 'intron'),])
# int_over = pbmclapply(1:length(gene_ranges), function(x) {size_overlap(gene_ranges[x], intron_boundaries)}, mc.cores = ncores)
# gene_ranges$intron_overlap = unlist(int_over)
# 
# 
# # DEBUG
# # df = as.data.frame(gene_ranges)
# cat("Weighted mean Rho for exons\n")
# # Weighted mean Rho overlap
# wmeanrho_overlap = function(x, y) {
#   # x is a GRanges object of size 1
#   # y is a GRanges object in which to search overlaps, e.g. introns or exons
#   overlap = findOverlaps(x, y)
#   if (length(overlap) > 0) {
#     query = queryHits(overlap)
#     subject = subjectHits(overlap)
#     wi = intersect(x, y[subject])
#     hits = subjectHits(findOverlaps(wi, ldmap_ranges))
#     rho = ldmap_ranges$Mean_rho[hits]
#     w = width(restrict(ldmap_ranges[hits], start(x), end(x)))
#     wmean = weighted.mean(rho, w)
#   } else {
#     wmean = NA
#   }
#   return(wmean)
# }
# # Mean Rho for exons
# # x = gene_ranges[32]
# # y = exon_boundaries
# exon_boundaries = makeGRangesFromDataFrame(gff[which(gff$feature == 'CDS'),])
# ex_over = pbmclapply(1:length(gene_ranges), function(x) {wmeanrho_overlap(gene_ranges[x], exon_boundaries)}, mc.cores = ncores)
# gene_ranges$wmean_exon = unlist(ex_over)
# 
# # DEBUG
# # for (x  in 1:1000) {
# #   wmeanrho_overlap(gene_ranges[x], exon_boundaries)
# # }
# 
# cat("Weighted mean Rho for introns\n")
# # Mean Rho for introns
# intron_boundaries = makeGRangesFromDataFrame(gff[which(gff$feature == 'intron'),])
# int_over = pbmclapply(1:length(gene_ranges), function(x) {wmeanrho_overlap(gene_ranges[x], intron_boundaries)}, mc.cores = ncores)
# gene_ranges$wmean_intron = unlist(int_over)



# Hotspot count ----
cat("Hotspot count\n")
source("Source/read.ldhot.R")

ldhot.raw = read.ldhot.all(max.length = 10^8,
                           peak.rate = 0,
                           intensity = 0,
                           max.intensity = 10^8)
ldhot.filtered.2 = read.ldhot.all(max.length = 10^4,
                                  peak.rate = 0,
                                  intensity = 0,
                                  max.intensity = 10^8)
ldhot.filtered.4 = read.ldhot.all(max.length = 10^4,
                                  peak.rate = 0,
                                  intensity = 4,
                                  max.intensity = 200)

# Reduce GFF to seqnames in LD maps
ldhot.raw = ldhot.raw[which(ldhot.raw$seqnames %in% map_list & ldhot.raw$dataset == dataset),]
ldhot.raw$seqnames = gff_list
ldhot.raw = makeGRangesFromDataFrame(ldhot.raw,
                                     keep.extra.columns=TRUE)

ldhot.filtered.2 = ldhot.filtered.2[which(ldhot.filtered.2$seqnames %in% map_list & ldhot.filtered.2$dataset == dataset),]
ldhot.filtered.2$seqnames = gff_list
ldhot.filtered.2 = makeGRangesFromDataFrame(ldhot.filtered.2,
                                            keep.extra.columns=TRUE)

ldhot.filtered.4 = ldhot.filtered.4[which(ldhot.filtered.4$seqnames %in% map_list & ldhot.filtered.4$dataset == dataset),]
ldhot.filtered.4$seqnames = gff_list
ldhot.filtered.4 = makeGRangesFromDataFrame(ldhot.filtered.4,
                                            keep.extra.columns=TRUE)

# How many hotspots overlap the position
hotOverlap = countOverlaps(gene_ranges, ldhot.raw)
table(hotOverlap)
gene_ranges$hotspot_overlap.raw = hotOverlap

hotOverlap = countOverlaps(gene_ranges, ldhot.filtered.2)
table(hotOverlap)
gene_ranges$hotspot_overlap.filtered.2 = hotOverlap

hotOverlap = countOverlaps(gene_ranges, ldhot.filtered.4)
table(hotOverlap)
gene_ranges$hotspot_overlap.filtered.4 = hotOverlap

# How many hotspots midpoints overlap the position
hotspots_midpoint = ldhot.raw
end(hotspots_midpoint) = start(hotspots_midpoint) + 1
hotOverlap = countOverlaps(gene_ranges, hotspots_midpoint)
table(hotOverlap)
gene_ranges$hotspot_count_midpoint.raw = hotOverlap

hotspots_midpoint = ldhot.filtered.2
end(hotspots_midpoint) = start(hotspots_midpoint) + 1
hotOverlap = countOverlaps(gene_ranges, hotspots_midpoint)
table(hotOverlap)
gene_ranges$hotspot_count_midpoint.filtered.2 = hotOverlap

hotspots_midpoint = ldhot.filtered.4
end(hotspots_midpoint) = start(hotspots_midpoint) + 1
hotOverlap = countOverlaps(gene_ranges, hotspots_midpoint)
table(hotOverlap)
gene_ranges$hotspot_count_midpoint.filtered.4 = hotOverlap

# How many hotspots on the chromosome
gene_ranges$n_hotspots_chromosome.raw = length(ldhot.raw)
gene_ranges$n_hotspots_chromosome.filtered.2 = length(ldhot.filtered.2)
gene_ranges$n_hotspots_chromosome.filtered.4 = length(ldhot.filtered.4)


# SAVE ----
cat("Save\n")

df_pos = as.data.frame(gene_ranges)

if (chromosome == "all") {
  save(df_pos, file = paste("Output/RhoGradient_5kbTSS_", dataset, ".Rda", sep = ""))
} else {
  save(df_pos, file = paste("Output/RhoGradient_5kbTSS_", dataset, "_", chromosome, ".Rda", sep = ""))
}

gc()

#============================================================================#
# End of script ----
#============================================================================#
