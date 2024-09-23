#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


########################################################################## #
#     JOB - GFF with population-scaled recombination rates
########################################################################## #
# test if there is at least one argument: if not, return an error
if (length(args) != 2) {
  stop("Two arguments must be supplied. <dataset> and <chromosome>.", call.=FALSE)
} else if (length(args) ==2 ) {
  dataset = args[1]
  chromosome = args[2]
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

# TODO keys to merge dataset (gff, LD maps) are
# dataset, Gene id, chromosome, positions

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

# source('Source/get_r.R')
source('Source/read.ldmap.R')
source('Source/read.ldhot.R')
source('Source/read.gff.R')


ncores = 16
# ncores = min(ncpus, detectCores())

cat(dataset, "\n")

# List of dataset
marey_data = read.table(paste("Data/Recombination/marey_dataset.csv", sep = ""), header = TRUE, sep = "\t")
chromosome_metadata = read.table(paste0("Data/Genome/genome_chromosome_metadata.csv"),
                                 header = TRUE, sep = "\t")

# List of chromosomes in the LD map
# map_list = chromosome_metadata$ldmapname[which(chromosome_metadata$set == set)]
# map_list = map_list[which(!is.na(map_list))]
map_list = chromosome
map_list
# List of chromosomes in the gff
# chr_list = chromosome_metadata$annotname[which(chromosome_metadata$set == dataset &
#                                                  chromosome_metadata$ldmapname == chromosome)]
# chr_list = chr_list[which(!is.na(chr_list))]
# chr_list

#============================================================================#
# GFF dataset ----
#============================================================================#

cat("Import gff...\n")
gfffile = paste("Data/Genomic_landscapes/GFF_parsed/", dataset, ".csv.gz", sep = "")
# gff = read.table(gzfile(gfffile), header = T, sep = "\t")
# gfffile = paste(wd, "/Data/Genomic_landscapes/GFF_parsed/", dataset, "_parsed.csv.gz", sep = "")

# Make the dataset - Pooled chromosomes
# gff = read.table(gzfile(gfffile), header = T, sep = "\t", quote = "",
#                  na.strings = c(NA, "NaN", ""), fill = TRUE)
gff = read.gff(gfffile)

gff_nongenic = gff[which(gff$feature %in% c("intergenic", "intergenic_buffer")),]


if (sum(gff$feature == "lnc_RNA") > 0) {
  gff = gff[-which(gff$feature == "lnc_RNA"),]  
}


# Filter protein coding genes only
# Do not filter when no protein coding information - all genes considered as protein coding
if ("protein_coding" %in% gff$gene_biotype) {
  pcgenes = gff$id[which(gff$gene_biotype == "protein_coding")]
  gff = gff[which(gff$id %in% pcgenes | gff$parent %in% pcgenes | gff$parent %in% gff$id[which(gff$parent %in% pcgenes)]),]
} else {
  gff = gff[-which(gff$feature %in% c("intergenic", "intergenic_buffer")),]
}

# The parsed GFF where recombination rates has been added
# load(file = paste("Output/gff_", dataset, ".Rda", sep = ""))

# LD MAPS ----
# Import LD maps
cat("Import LD maps\n")

bpen = 5
ldmap = read.ldmap(dataset, map_list)

# Change chromosome names to GFF annotation names
# for (chr in unique(ldmap$chromosome)) {
#   # cat(chr, '\n')
#   ldmap$chromosome[which(ldmap$chromosome == chr)] = chromosome_metadata$annotname[which(chromosome_metadata$ldmapname == chr & chromosome_metadata$set == dataset)]
# }
# unique(ldmap$chromosome)
ldmap_ranges = makeGRangesFromDataFrame(ldmap, keep.extra.columns = TRUE)


# TODO Reduce GFF to seqnames in LD maps
gff = gff[which(gff$seqname %in% unique(ldmap$chromosome)),]


# Take care of improper coordinates
# e.g. start > (end - 1)
gff = gff[which(gff$start < (gff$end - 1)),]

# gff = gff[-which(grepl("stream", gff$feature)),]

# TODO GFF Features
# - introns in UTR
# - mono-exonic
# - gene overlap hotspot/hotspot position
# - Number of SNPs
# - Number of exons of the gene
# - ID of the gene
# - Mean Rho/weighted mean Rho

# A subset for developping
# gff = gff[1:10000,]
cat("Make GFF genomic ranges\n")
gff_ranges = makeGRangesFromDataFrame(gff, keep.extra.columns=TRUE)


# LD hotspots ----
cat("Import LD hotspots\n")
# Import LD hotspots
# Pool LD maps
bpen = 5
# Raw data
# ldhot = lapply(1:nrow(chromosome_metadata),
#                   function(x) {read.ldhot(dataset,
#                                           chromosome,
#                                           max.length = 1000000,
#                                           peak.rate = 1)})
# ldhot.filtered = lapply(1:nrow(chromosome_metadata),
#                function(x) {read.ldhot(dataset,
#                                        chromosome,
#                                        max.length = 10000,
#                                        peak.rate = 1)})
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


ldhot.raw = makeGRangesFromDataFrame(ldhot.raw,
                                            keep.extra.columns=TRUE)

midpoint = (start(ldhot.raw) + end(ldhot.raw))/2
ldhot.raw$midpoint = midpoint

ldhot.filtered.2 = makeGRangesFromDataFrame(ldhot.filtered.2,
                                            keep.extra.columns=TRUE)

midpoint = (start(ldhot.filtered.2) + end(ldhot.filtered.2))/2
ldhot.filtered.2$midpoint = midpoint

ldhot.filtered.4 = makeGRangesFromDataFrame(ldhot.filtered.4,
                                            keep.extra.columns=TRUE)

midpoint = (start(ldhot.filtered.4) + end(ldhot.filtered.4))/2
ldhot.filtered.4$midpoint = midpoint







#============================================================================#
# Retrieve promoter regions with GenomicRanges ----
# cat("Promoters\n")

# rm(gff)
# TODO Add sequence length
# chr_size = data.frame(seqnames = levels(seqnames(gff_ranges)))
# for (i in 1:nrow(chr_size)) {
#   chr_size$seqlength[i] = chromosome_metadata$chrsize.bp[which(chromosome_metadata$set == set &
#                                                                  chromosome_metadata$annotname == chr_size$seqnames[i])]
# }
# seqlengths(gff_ranges) = chr_size$seqlength
# 
# 
# if ('gene' %in% unique(gff_ranges$feature)) {
#   prom = promoters(gff_ranges[which(gff_ranges$feature == 'gene')], upstream = 2000, downstream = 0)
#   prom = trim(prom)
#   prom$feature = 'promoter'
#   prom$parent = prom$id
#   gff_ranges = c(gff_ranges, prom)
#   rm(prom)
# }



#============================================================================#
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
if (dataset == "Citrullus_lanatus_Guo2019") {
  gff_ranges$gene_id = list_ids$gene
} else {
  genids = pbmclapply(1:length(gff_ranges), function(x){genid(x)})

  if (length(unlist(genids)) != length(gff_ranges)) {
    warning("Parsing gene ids introduced missing rows")
  }
  gff_ranges$gene_id = unlist(genids)
}


#============================================================================#
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



#============================================================================#
# 3. Assess if intron in 5'/3' UTRs (TRUE/FALSE) ----
cat("Introns in UTRs\n")

intron_utr = function(x) {
  if ((gff_ranges$feature[x] == "gene") | (!("gene" %in% unique(gff_ranges$feature)) & (gff_ranges$feature[x] == "mRNA"))) {
    gene = gff_ranges[which(gff_ranges$gene_id == gff_ranges$gene_id[x]),]
    if (unique(as.character(strand(gene[gene$feature == "gene"])))  == "+") {
      intron5 = min(start(gene)[gene$feature == "intron"]) < min(start(gene)[gene$feature == "CDS"])
      intron3 = max(start(gene)[gene$feature == "intron"]) > max(start(gene)[gene$feature == "CDS"])
    } else {
      intron5 = max(start(gene)[gene$feature == "intron"]) > max(start(gene)[gene$feature == "CDS"])
      intron3 = min(start(gene)[gene$feature == "intron"]) < min(start(gene)[gene$feature == "CDS"])
    }
    introns = data.frame(intron5 = intron5, intron3 = intron3)
  } else {
    introns = data.frame(intron5 = NA, intron3 = NA)
  }
  return(introns)
}

# DEBUG
# for (i in 1:length(gff_ranges)) {
#   cat(i, "\n")
#   intron_utr(i)
# }

if (dataset != "Citrullus_lanatus_Guo2019") {
  intronsutr = pbmclapply(1:length(gff_ranges), function(x){intron_utr(x)})
  intronsutr = rbindlist(intronsutr)
  if (nrow(intronsutr) != length(gff_ranges)) {
    warning("Parsing introns in UTRs introduced missing rows")
  }
  gff_ranges$intron5utr = intronsutr$intron5
  gff_ranges$intron3utr = intronsutr$intron3
} else {
  gff_ranges$intron5utr = NA
  gff_ranges$intron3utr = NA
}





#============================================================================#
# 4. Distance to closest hotspot (bp) ----
# Take the midpoint of hotspots

# hotpos = makeGRangesFromDataFrame(data.frame(seqnames = as.character(seqnames(ldhot.filtered.2)),
#                                              start = ldhot.filtered.2$midpoint,
#                                              end = (ldhot.filtered.2$midpoint + 1)))
# nearest_hotspot = distanceToNearest(gff_ranges, hotpos, ignore.strand = TRUE)

# gff_ranges$nearest_hotspot_bp = NA
# gff_ranges$nearest_hotspot_bp[nearest_hotspot@from] = mcols(nearest_hotspot)$distance
# 
# 
# # 4.b. Distance to closest hotspot (bp) ----
# # Take the midpoint of hotspots
# # Take the start of each Grange
# # Consider only distance from start position of the feature
# gff_ranges_copy = gff_ranges
# end(gff_ranges_copy[which(strand(gff_ranges_copy) == "+")]) = start(gff_ranges_copy[which(strand(gff_ranges_copy) == "+")]) + 1
# start(gff_ranges_copy[which(strand(gff_ranges_copy) == "-")]) = end(gff_ranges_copy[which(strand(gff_ranges_copy) == "-")]) - 1
# nearest_hotspot = distanceToNearest(gff_ranges_copy, hotpos, ignore.strand = TRUE)
# 
# gff_ranges$nearest_hotspot_start_bp = NA
# gff_ranges$nearest_hotspot_start_bp[nearest_hotspot@from] = mcols(nearest_hotspot)$distance
# 
# 
# # 4.c. Distance to closest hotspot (bp) ----
# # Take the midpoint of hotspots
# # Take the end of each Grange
# # Consider only distance from start position of the feature
# gff_ranges_copy = gff_ranges
# start(gff_ranges_copy[which(strand(gff_ranges_copy) == "+")]) = end(gff_ranges_copy[which(strand(gff_ranges_copy) == "+")]) - 1
# end(gff_ranges_copy[which(strand(gff_ranges_copy) == "-")]) = start(gff_ranges_copy[which(strand(gff_ranges_copy) == "-")]) + 1
# nearest_hotspot = distanceToNearest(gff_ranges_copy, hotpos, ignore.strand = TRUE)
# 
# gff_ranges$nearest_hotspot_end_bp = NA
# gff_ranges$nearest_hotspot_end_bp[nearest_hotspot@from] = mcols(nearest_hotspot)$distance
# 

# #============================================================================#
# # 5. Add intergenic regions ----
cat("Intergenic regions\n")

# gff = rbind(gff, gff_nongenic)
intergenic = makeGRangesFromDataFrame(gff_nongenic,
                                      keep.extra.columns = TRUE)

#
# if ('gene' %in% unique(gff_ranges$feature)) {
#   genes = gff_ranges[which(gff_ranges$feature == "gene")]
# } else {
#   if ('mRNA' %in% unique(gff_ranges$feature)) {
#     genes = gff_ranges[which(gff_ranges$feature == "mRNA")]
#   }
# }
#
# chr_ranges = makeGRangesFromDataFrame(data.frame(seqnames = as.character(levels(seqnames(gff_ranges))),
#                                                  start = rep(1, length(seqlengths(gff_ranges))),
#                                                  end = as.numeric(seqlengths(gff_ranges)), strand = "*"))
# collapsed_genes = reduce(genes)
# strand(collapsed_genes) = "*"
# intergenic = setdiff(chr_ranges, collapsed_genes)
# intergenic$id = paste0("intergenic-", seq(1, length(intergenic)))
# intergenic$feature = "intergenic"
#
gff_ranges = c(gff_ranges, intergenic)



# # 6. Number of SNPs ----
# cat("Number of SNPs\n")
#
# # SNPs are breakpoints in the LD map, i.e. each start position is a SNP
#
# # Make SNP set
# snps = data.frame(seqnames = seqnames(ldmap_ranges),
#                   start = start(ldmap_ranges),
#                   end = start(ldmap_ranges) + 1)
# snps = makeGRangesFromDataFrame(snps)
# # SNP count is simply the number of overlaps
# snpcount = countOverlaps(gff_ranges, snps, type = "any")
# # gff_ranges[1]
# # snps[23:27]
#
# cat("Summary of SNP count\n")
# cat("-----------------------\n")
# table(snpcount)
# cat("-----------------------\n")
# cat("\n")
#
# gff_ranges$snp_count = snpcount


#============================================================================#
# 7. Hotspot overlaps ----
cat("Hotspot overlap\n")
# Trimmed hotspots
# ldhot.filtered = makeGRangesFromDataFrame(ldhot.filtered,
#                                  keep.extra.columns=TRUE)
hotoverlap = countOverlaps(gff_ranges, ldhot.filtered.2, type = "any")

if (length(hotoverlap) != length(gff_ranges)) {
  warning("Parsing hotspot overlap introduced missing rows")
}
gff_ranges$hotspot_overlap_intensity2 = hotoverlap


hotoverlap = countOverlaps(gff_ranges, ldhot.filtered.4, type = "any")

if (length(hotoverlap) != length(gff_ranges)) {
  warning("Parsing hotspot overlap introduced missing rows")
}
gff_ranges$hotspot_overlap_intensity4 = hotoverlap


hotoverlap = countOverlaps(gff_ranges, ldhot.raw, type = "any")

if (length(hotoverlap) != length(gff_ranges)) {
  warning("Parsing hotspot overlap introduced missing rows")
}
gff_ranges$hotspot_overlap_raw = hotoverlap



# 4b. Hotspot coverage - How much of the sequence is covered by a hostpot
# Size of intersect
hits = findOverlaps(gff_ranges, ldhot.raw, ignore.strand = TRUE)
grl = extractList(ldhot.raw, as(hits, "List"))
hotcoverage = width(pintersect(gff_ranges, grl, ignore.strand = TRUE))
hotcoverage[which(unlist(lapply(hotcoverage, length)) == 0)] = 0
hotcoverage = lapply(hotcoverage, function(x){mean(unlist(x))})
hotcoverage = unlist(hotcoverage)

if (length(hotcoverage) != length(gff_ranges)) {
  warning("Parsing hotspot coverage introduced missing rows")
}
gff_ranges$hotspot_coverage_bp_raw = hotcoverage


hits = findOverlaps(gff_ranges, ldhot.filtered.2, ignore.strand = TRUE)
grl = extractList(ldhot.filtered.2, as(hits, "List"))
hotcoverage = width(pintersect(gff_ranges, grl, ignore.strand = TRUE))
hotcoverage[which(unlist(lapply(hotcoverage, length)) == 0)] = 0
hotcoverage = lapply(hotcoverage, function(x){mean(unlist(x))})
hotcoverage = unlist(hotcoverage)

if (length(hotcoverage) != length(gff_ranges)) {
  warning("Parsing hotspot coverage introduced missing rows")
}
gff_ranges$hotspot_coverage_bp_intensity2 = hotcoverage


hits = findOverlaps(gff_ranges, ldhot.filtered.4, ignore.strand = TRUE)
grl = extractList(ldhot.filtered.4, as(hits, "List"))
hotcoverage = width(pintersect(gff_ranges, grl, ignore.strand = TRUE))
hotcoverage[which(unlist(lapply(hotcoverage, length)) == 0)] = 0
hotcoverage = lapply(hotcoverage, function(x){mean(unlist(x))})
hotcoverage = unlist(hotcoverage)

if (length(hotcoverage) != length(gff_ranges)) {
  warning("Parsing hotspot coverage introduced missing rows")
}
gff_ranges$hotspot_coverage_bp_intensity4 = hotcoverage




#============================================================================#
# 8. Weighted mean Rho ----
#============================================================================#

cat("Recombination rates\n")

# DEBUG
# gff_ranges = gff_ranges[1:10000]


# Sampled ranges ----
sample_ranges = gff_ranges

# Random ranges for control ----
set.seed(42)
reshuffling = sample(10000:(max(start(sample_ranges)) - 10000), length(sample_ranges))
random_ranges = GRanges(seqnames = seqnames(sample_ranges), strand = strand(sample_ranges),
        ranges = IRanges(start = reshuffling, width = width(sample_ranges)))


# 9. Estimating mean & median Rho ----
cat("Mean Rho\n")
# index = seq(1, length(sample_ranges))
# res = split(sample_ranges, as.factor(index))
hits = findOverlaps(sample_ranges, ldmap_ranges)
# hits2 = split(hits, as.factor(queryHits(hits)))
query = queryHits(hits)
subj = subjectHits(hits)

mean.rho = unlist(pbmclapply(1:length(sample_ranges), function(x) {if (x %in% query) {if (length(subj[which(query == x)]) > 0) {mean(ldmap_ranges$Mean_rho[subj[which(query == x)]], na.rm = TRUE)} else {NA}} else {NA}}, mc.cores = ncores))
# mean.rho = unlist(pbmclapply(1:length(res), function(x) {if (x %in% names(hits2)) {if (length(subjectHits(hits2[[which(names(hits2) == x)]])) > 0) {mean(ldmap_ranges$Mean_rho[subjectHits(hits2[[which(names(hits2) == x)]])], na.rm = TRUE)} else {NA}} else {NA}}, mc.cores = ncores))
median.rho = unlist(pbmclapply(1:length(sample_ranges), function(x) {if (x %in% query) {if (length(subj[which(query == x)]) > 0) {median(ldmap_ranges$Mean_rho[subj[which(query == x)]], na.rm = TRUE)} else {NA}} else {NA}}, mc.cores = ncores))

cat("Weighted mean\n")
# Restrict weights to sample range
# Weighted mean Rho overlap
wmeanrho = function(x) {
  # x is a GRanges object of size 1
  hits = subjectHits(findOverlaps(x, ldmap_ranges))
  rho = ldmap_ranges$Mean_rho[hits]
  w = width(restrict(ldmap_ranges[hits], start(x), end(x)))
  wmean = weighted.mean(rho, w)
  if (is.numeric(wmean) & length(wmean) == 1 & !is.nan(wmean)) {
    return(wmean)
  } else {
    return(NA)
  }
}
wmean.rho = pbmclapply(1:length(gff_ranges), function(x) {wmeanrho(gff_ranges[x])})

# w = pbmclapply(1:length(res), function(x) {if (x %in% query) {if (length(subj[which(query == x)]) > 0) {ldmap_ranges[subj[which(query == x)]]} else {NA}} else {NA}}, mc.cores = ncores)
# w2 = pbmclapply(1:length(w), function(x) {if(length(w[[x]]) == 0 | length(res[[x]]) == 0 | sum(is.na(w[[x]])) > 0 | sum(is.na(res[[x]])) > 0) {NA} else {pintersect(w[[x]], res[[x]])}})
# w3 = pbmclapply(1:length(w2), function(x) {if(length(w2[[x]]) == 0 | sum(is.na(w2[[x]])) > 0 ) {NA} else {width(w2[[x]])}})
# w3 = pbmclapply(1:length(w2), function(x) {width(w2[[x]])})

# w[175847:181096]
# for (x in 175847:181096) {
#   if(length(w[[x]]) == 0 | length(res[[x]]) == 0 | sum(is.na(w[[x]])) > 0 | sum(is.na(res[[x]])) > 0) {NA} else {pintersect(w[[x]], res[[x]])}
# }
# x = 179219

# wmean.rho = unlist(pbmclapply(1:length(res), function(x) {if (x %in% query) {if (length(subj[which(query == x)]) > 0) {weighted.mean(x = ldmap_ranges$Mean_rho[subj[which(query == x)]], w = w3[[x]], na.rm = TRUE)} else {NA}} else {NA}}, mc.cores = ncores))

# rm(w)
# rm(w2)
# rm(w3)

# Mean Rho fro control
cat("Control\n")
# index = seq(1, length(random_ranges))
# res = split(random_ranges, as.factor(index))
hits = findOverlaps(random_ranges, ldmap_ranges)
# hits2 = split(hits, as.factor(queryHits(hits)))
query = queryHits(hits)
subj = subjectHits(hits)

mean.rho.control = unlist(pbmclapply(1:length(random_ranges), function(x) {if (x %in% query) {if (length(subj[which(query == x)]) > 0) {mean(ldmap_ranges$Mean_rho[subj[which(query == x)]], na.rm = TRUE)} else {NA}} else {NA}}, mc.cores = ncores))
median.rho.control = unlist(pbmclapply(1:length(random_ranges), function(x) {if (x %in% query) {if (length(subj[which(query == x)]) > 0) {median(ldmap_ranges$Mean_rho[subj[which(query == x)]], na.rm = TRUE)} else {NA}} else {NA}}, mc.cores = ncores))

# Weigthed mean Rho for control
cat("Weighted mean control\n")
# w = pbmclapply(1:length(res), function(x) {if (x %in% query) {if (length(subj[which(query == x)]) > 0) {ldmap_ranges[subj[which(query == x)]]} else {NA}} else {NA}}, mc.cores = ncores)
# w2 = pbmclapply(1:length(w), function(x) {if (length(w[[x]]) > 0 & length(res[[x]])) {pintersect(w[[x]], res[[x]])} else {NA}})
# w3 = pbmclapply(1:length(w2), function(x) {width(w2[[x]])})

wmean.rho.control = unlist(pbmclapply(1:length(random_ranges), function(x) {wmeanrho(random_ranges[x])}))
# wmean.rho.control = unlist(pbmclapply(1:length(res), function(x) {if (x %in% query) {if (length(subj[which(query == x)]) > 0) {weighted.mean(x = ldmap_ranges$Mean_rho[subj[which(query == x)]], w = w3[[x]], na.rm = TRUE)} else {NA}} else {NA}}, mc.cores = ncores))

gff_ranges$mean.rho = mean.rho
gff_ranges$median.rho = median.rho
gff_ranges$mean.rho.control = mean.rho.control
gff_ranges$median.rho.control = median.rho.control
gff_ranges$weighted.mean.rho = unlist(wmean.rho)
gff_ranges$weighted.mean.rho.control = wmean.rho.control


# 10. Estimating mean Rho (strict mode) ----
# Keep only LD map windows with a start after the start of the feature (i>e. start position overlapping)
# To eliminate any artefactual gradient due to overlapping large hotspot windows
cat("Mean Rho (strict mode)\n")
# index = seq(1, length(sample_ranges))
# res = split(sample_ranges, as.factor(index))
hits = findOverlaps(sample_ranges, ldmap_ranges, type = "any")
query = queryHits(hits)
subj = subjectHits(hits)
# Remove hits which have a start before the query start
hits_startafter = hits[-which(start(ldmap_ranges)[subj] < start(gff_ranges)[query])]
query = queryHits(hits_startafter)
subj = subjectHits(hits_startafter)

mean.rho.startafter = unlist(pbmclapply(1:length(sample_ranges),
                                    function(x) {if (x %in% query) {
                                      if (length(subj[which(query == x)]) > 0) {
                                        mean(ldmap_ranges$Mean_rho[subj[which(query == x)]], na.rm = TRUE)}
                                      else {NA}} else {NA}}, mc.cores = ncores))
# sample_ranges[1]
# mean.rho.startafter[1]
# ldmap_ranges[subj[which(query == 1)]]


# cat("Weighted mean (strict mode)\n")
# # Restrict weigths to sample range
# # Weighted mean Rho overlap
wmeanrho.strict = function(x) {
  # x is a GRanges object of size 1
  hits = findOverlaps(x, ldmap_ranges, type = "any")
  subj = subjectHits(hits)
  # Remove hits which have a start before the query start
  hits_startafter = hits[-which(start(ldmap_ranges)[subj] < start(x))]
  subj = subjectHits(hits_startafter)
  # hits = subjectHits(findOverlaps(x, ldmap_ranges, type = "any"))
  rho = ldmap_ranges$Mean_rho[subj]
  w = width(restrict(ldmap_ranges[subj], start(x), end(x)))
  wmean = weighted.mean(rho, w)
  if (is.numeric(wmean) & length(wmean) == 1 & !is.nan(wmean)) {
    return(wmean)
  } else {
    return(NA)
  }
}

weighted.mean.rho.startafter = unlist(pbmclapply(1:length(gff_ranges), function(x) {wmeanrho.strict(gff_ranges[x])}))


gff_ranges$mean.rho.startafter = mean.rho.startafter
gff_ranges$weighted.mean.rho.startafter = weighted.mean.rho.startafter


# 11. Relative Rho values calculated for each gene ----
# Debug
cat("Relative Rho\n")

# TEST
# Empirical values
# saveRDS(gff_ranges, file = "Output/test.rds")
# gff_ranges = readRDS("Output/test.rds")
# Random values
# gff_ranges$weighted.mean.rho = runif(n = length(gff_ranges), min = 0, max = 10)


# Normalize by the Rho value of the gene.
# Divide each row by the value of its gene.
relative.meanrho = unlist(pbmclapply(1:length(gff_ranges), function(x) {ifelse(length(gff_ranges$mean.rho[x]) > 0 & length(which(gff_ranges$id == gff_ranges$gene_id[x] & gff_ranges$feature == "gene")) > 0,
                                                                               gff_ranges$mean.rho[x]/mean(gff_ranges$mean.rho[which(gff_ranges$id == gff_ranges$gene_id[x] & gff_ranges$feature == "gene")]),
                                                                               NA)
}))

# relative.medianrho = unlist(pbmclapply(1:length(gff_ranges), function(x) {ifelse(length(gff_ranges$mean.rho[x]) > 0 & length(which(gff_ranges$id == gff_ranges$gene_id[x] & gff_ranges$feature == "gene")) > 0,
#                                                                                  gff_ranges$median.rho[x]/mean(gff_ranges$median.rho[which(gff_ranges$id == gff_ranges$gene_id[x] & gff_ranges$feature == "gene")]),
#                                                                                  NA)
# }))

relative.wmeanrho = unlist(pbmclapply(1:length(gff_ranges), function(x) {ifelse(length(gff_ranges$mean.rho[x]) > 0 & length(which(gff_ranges$id == gff_ranges$gene_id[x] & gff_ranges$feature == "gene")) > 0,
                                                                                gff_ranges$weighted.mean.rho[x]/mean(gff_ranges$weighted.mean.rho[which(gff_ranges$id == gff_ranges$gene_id[x] & gff_ranges$feature == "gene")]),
                                                                                NA)
}))

gff_ranges$relative.meanrho = relative.meanrho
# gff_ranges$relative.medianrho = relative.medianrho
gff_ranges$relative.wmeanrho = relative.wmeanrho



# 12. Relative position from ATG ----
cat("Relative position from ATG\n")
dist_atg = function(x) {
  # Compute the distance to the ATG position
  # Take an index of a row in a gff_ranges object in argument
  atg = start(gff_ranges)[which(gff_ranges$gene_id == gff_ranges$gene_id[x] &
                                  gff_ranges$feature == "CDS" & gff_ranges$rank == 1)]
  if (length(atg) > 0) {
    atg = max(atg) # Take the longest distance (longest transcript)
    st = start(gff_ranges)[x]
    if (as.character(strand(gff_ranges)[x]) == "+") {
      d_atg = st - atg
    }
    if (as.character(strand(gff_ranges)[x]) == "-") {
      d_atg = atg - st
    }
  } else {
    d_atg = NA
  }
  return(d_atg)
}
# x = genes[[1]]
# dist_atg(1)

atg = pbmclapply(1:length(gff_ranges), function(x) {dist_atg(x)})

gff_ranges$dist_atg = unlist(atg)

# for (x in 255200:257900) {
#   dist_atg(x)
# }
# x = 256278
# dist_atg(x)



# Save ----
cat("Save ", paste("Output/gff_rho_", dataset, "_", chromosome, ".Rda", sep = ""),"...\n")
save(gff_ranges, file = paste("Output/gff_rho_", dataset, "_", chromosome, ".Rda", sep = ""))

df = as.data.frame(gff_ranges, row.names = NULL)
# df$weighted.mean.rho.control = unlist(df$weighted.mean.rho.control)

cat("Save .csv file\n")
write.table(df, file = gzfile(paste("Data/Genomic_landscapes/Rho/gff_rho_", dataset, "_", chromosome, ".csv.gz", sep = "")),
            row.names = F, col.names = T, sep = "\t", quote = F)

# load(file = paste("Output/gff_rho_", dataset, "_", chromosome, ".Rda", sep = ""))



#============================================================================#
# End of script ----
#============================================================================#
