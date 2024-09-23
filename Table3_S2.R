# Loading env
source("Source/init_knitr.R")
source("Source/set_ggplot_theme.R")

# BiocManager::install("nullranges")
# library(nullranges)


ldhot_trimmed_nofiltering = read.ldhot.all(max.length = 10^8,
                               peak.rate = 0,
                               intensity = 0,
                               max.intensity = 10^8)
ldhot_trimmed_soft = read.ldhot.all(max.length = 10^4,
                               peak.rate = 0,
                               intensity = 0,
                               max.intensity = 10^8)
ldhot_trimmed_hard = read.ldhot.all(max.length = 10^4,
                               peak.rate = 0,
                               intensity = 4,
                               max.intensity = 200)

gff = readRDS("Data/Genome/gff_all.rds")


ds = sort(as.character(unique(ldhot_trimmed_hard$dataset)))
sp = ds
for (i in 1:length(sp)) {
  sp[i] = dataset_metadata$species[which(dataset_metadata$dataset == sp[i])]
}

df = data.frame(dataset = ds,
                species = sp,
                n.hotspots.genic = NA,
                n.hotspots.intergenic = NA,
                expected.hotspots.genic = NA,
                expected.hotspots.genic.lowerCI = NA,
                expected.hotspots.genic.upperCI = NA,
                expected.hotspots.intergenic = NA,
                expected.hotspots.intergenic.lowerCI = NA,
                expected.hotspots.intergenic.upperCI = NA,
                n.hotspots.genic.filtered2 = NA,
                n.hotspots.intergenic.filtered2 = NA,
                expected.hotspots.genic.filtered2 = NA,
                expected.hotspots.genic.filtered2.lowerCI = NA,
                expected.hotspots.genic.filtered2.upperCI = NA,
                expected.hotspots.intergenic.filtered2 = NA,
                expected.hotspots.intergenic.filtered2.lowerCI = NA,
                expected.hotspots.intergenic.filtered2.upperCI = NA,
                n.hotspots.genic.filtered4 = NA,
                n.hotspots.intergenic.filtered4 = NA,
                expected.hotspots.genic.filtered4 = NA,
                expected.hotspots.genic.filtered4.lowerCI = NA,
                expected.hotspots.genic.filtered4.upperCI = NA,
                expected.hotspots.intergenic.filtered4 = NA,
                expected.hotspots.intergenic.filtered4.lowerCI = NA,
                expected.hotspots.intergenic.filtered4.upperCI = NA)

for (i in 1:nrow(df)) {
  species = df$species[i]
  dat = df$dataset[i]
  cat(species, "\n")
  
  genes = data_all[which(data_all$species == species & data_all$feature == "gene"),]
  hot = makeGRangesFromDataFrame(ldhot_trimmed_nofiltering[which(ldhot_trimmed_nofiltering$dataset == dat),],
                                       ignore.strand = T)
  # Remove discarded chromosomes in GFF
  genes = genes[(genes$chromosome %in% seqnames(hot)),]  
    
    
  # Remove D sub-genomes in Wheat
  genes = genes[!(genes$chromosome %in% c("chr1D", "chr2D", "chr3D", "chr4D", "chr5D", "chr6D", "chr7D",
  "chr1A", "chr2A", "chr3A", "chr4A", "chr5A", "chr6A", "chr7A")),]
  
  gffRanges = makeGRangesFromDataFrame(genes,
                                       ignore.strand = T)
  
 
  
  # UNFILTERED HOTSPOTS
  hot = makeGRangesFromDataFrame(ldhot_trimmed_nofiltering[which(ldhot_trimmed_nofiltering$dataset == dat),],
                                       ignore.strand = T)
  
  # Remove discarded chromosomes
  gffRanges = gffRanges[which(seqnames(gffRanges) %in% seqnames(hot))]
  
  seqlen = aggregate(end ~ seqnames, as.data.frame(gffRanges), max)
  colnames(seqlen)[1] = "chromosome"
  
  for (j in 1:nrow(seqlen)) {
    seqlen$length[j] = chromosome_metadata$chrsize.bp[which(chromosome_metadata$set == dat & chromosome_metadata$ldmapname == seqlen$chromosome[j])]
  }
  
  seqlengths(gffRanges) = seqlen$length
  
  nbgenes = aggregate(end ~ chromosome, genes, length)
  
  # Count overlaps
  hits = countOverlaps(hot, gffRanges)
  # table(hits)
  df$n.hotspots.intergenic[i] = sum(hits == 0)
  df$n.hotspots.genic[i] = sum(hits > 0)
  
  # Randomize genes
  nboot = 1000
  boot = data.frame(genic = numeric(nboot),
                    intergenic = numeric(nboot))
  
  for (j in 1:nboot) {
    # cat(j, "\n")
    # resampling genes
    # newpos = unlist(lapply(1:nrow(seqlen), function(x) {sample(1:seqlen$end[x], nbgenes$end[x], replace = F)}))
    # bootRanges = GRanges(seqnames = seqnames(gffRanges),
    #                        IRanges(start = newpos,
    #                                width = width(gffRanges)))
    
    # resampling hotspots
    nbhot = aggregate(end ~ seqnames, as.data.frame(hot), length)

    newpos = unlist(lapply(1:nrow(seqlen), function(x) {sample(1:seqlen$end[x], nbhot$end[x], replace = F)}))
    bootRanges = GRanges(seqnames = seqnames(hot),
                           IRanges(start = newpos,
                                   width = width(hot)))
    
    
    # Count expected overlap
    hits = countOverlaps(hot, bootRanges)
    # table(hits)
    boot$intergenic[j] = sum(hits == 0)
    boot$genic[j] = sum(hits > 0)
  }
  
  df$expected.hotspots.genic[i] = mean(boot$genic)
  df$expected.hotspots.genic.lowerCI[i] = quantile(boot$genic, c(0.025))
  df$expected.hotspots.genic.upperCI[i] = quantile(boot$genic, c(0.975))
  
  df$expected.hotspots.intergenic[i] = mean(boot$intergenic)
  df$expected.hotspots.intergenic.lowerCI[i] = quantile(boot$intergenic, c(0.025))
  df$expected.hotspots.intergenic.upperCI[i] = quantile(boot$intergenic, c(0.975))
  
  
  # SOFT FILTERED HOTSPOTS
  hot = makeGRangesFromDataFrame(ldhot_trimmed_soft[which(ldhot_trimmed_soft$dataset == dat),],
                                       ignore.strand = T)
  
  # Remove discarded chromosomes
  gffRanges = gffRanges[which(seqnames(gffRanges) %in% seqnames(hot))]
  
  seqlen = aggregate(end ~ seqnames, as.data.frame(gffRanges), max)
  colnames(seqlen)[1] = "chromosome"
  
  for (j in 1:nrow(seqlen)) {
    seqlen$length[j] = chromosome_metadata$chrsize.bp[which(chromosome_metadata$set == dat & chromosome_metadata$ldmapname == seqlen$chromosome[j])]
  }
  
  seqlengths(gffRanges) = seqlen$length
  
  nbgenes = aggregate(end ~ chromosome, genes, length)
  
  # Count overlaps
  hits = countOverlaps(hot, gffRanges)
  # table(hits)
  df$n.hotspots.intergenic.filtered2[i] = sum(hits == 0)
  df$n.hotspots.genic.filtered2[i] = sum(hits > 0)
  
  # Randomize genes
  nboot = 1000
  boot = data.frame(genic = numeric(nboot),
                    intergenic = numeric(nboot))
  
  for (j in 1:nboot) {
    # cat(j, "\n")
    # resampling genes
    # newpos = unlist(lapply(1:nrow(seqlen), function(x) {sample(1:seqlen$end[x], nbgenes$end[x], replace = F)}))
    # bootRanges = GRanges(seqnames = seqnames(gffRanges),
    #                        IRanges(start = newpos,
    #                                width = width(gffRanges)))
    
    # resampling hotspots
    nbhot = aggregate(end ~ seqnames, as.data.frame(hot), length)

    newpos = unlist(lapply(1:nrow(seqlen), function(x) {sample(1:seqlen$end[x], nbhot$end[x], replace = F)}))
    bootRanges = GRanges(seqnames = seqnames(hot),
                           IRanges(start = newpos,
                                   width = width(hot)))
    
    
    # Count expected overlap
    hits = countOverlaps(hot, bootRanges)
    # table(hits)
    boot$intergenic[j] = sum(hits == 0)
    boot$genic[j] = sum(hits > 0)
  }
  
  df$expected.hotspots.genic.filtered2[i] = mean(boot$genic)
  df$expected.hotspots.genic.filtered2.lowerCI[i] = quantile(boot$genic, c(0.025))
  df$expected.hotspots.genic.filtered2.upperCI[i] = quantile(boot$genic, c(0.975))
  
  df$expected.hotspots.intergenic.filtered2[i] = mean(boot$intergenic)
  df$expected.hotspots.intergenic.filtered2.lowerCI[i] = quantile(boot$intergenic, c(0.025))
  df$expected.hotspots.intergenic.filtered2.upperCI[i] = quantile(boot$intergenic, c(0.975))


  # HARD FILTERED HOTSPOTS
  hot = makeGRangesFromDataFrame(ldhot_trimmed_hard[which(ldhot_trimmed_hard$dataset == dat),],
                                       ignore.strand = T)
  
  # Remove discarded chromosomes
  gffRanges = gffRanges[which(seqnames(gffRanges) %in% seqnames(hot))]
  
  seqlen = aggregate(end ~ seqnames, as.data.frame(gffRanges), max)
  colnames(seqlen)[1] = "chromosome"
  
  for (j in 1:nrow(seqlen)) {
    seqlen$length[j] = chromosome_metadata$chrsize.bp[which(chromosome_metadata$set == dat & chromosome_metadata$ldmapname == seqlen$chromosome[j])]
  }
  
  seqlengths(gffRanges) = seqlen$length
  
  nbgenes = aggregate(end ~ chromosome, genes, length)
  
  # Count overlaps
  hits = countOverlaps(hot, gffRanges)
  # table(hits)
  df$n.hotspots.intergenic.filtered4[i] = sum(hits == 0)
  df$n.hotspots.genic.filtered4[i] = sum(hits > 0)
  
  # Randomize genes
  nboot = 1000
  boot = data.frame(genic = numeric(nboot),
                    intergenic = numeric(nboot))
  
  for (j in 1:nboot) {
    # cat(j, "\n")
    # resampling genes
    # newpos = unlist(lapply(1:nrow(seqlen), function(x) {sample(1:seqlen$end[x], nbgenes$end[x], replace = F)}))
    # bootRanges = GRanges(seqnames = seqnames(gffRanges),
    #                        IRanges(start = newpos,
    #                                width = width(gffRanges)))
    
    # resampling hotspots
    nbhot = aggregate(end ~ seqnames, as.data.frame(hot), length)

    newpos = unlist(lapply(1:nrow(seqlen), function(x) {sample(1:seqlen$end[x], nbhot$end[x], replace = F)}))
    bootRanges = GRanges(seqnames = seqnames(hot),
                           IRanges(start = newpos,
                                   width = width(hot)))
    
    
    # Count expected overlap
    hits = countOverlaps(hot, bootRanges)
    # table(hits)
    boot$intergenic[j] = sum(hits == 0)
    boot$genic[j] = sum(hits > 0)
  }
  
  df$expected.hotspots.genic.filtered4[i] = mean(boot$genic)
  df$expected.hotspots.genic.filtered4.lowerCI[i] = quantile(boot$genic, c(0.025))
  df$expected.hotspots.genic.filtered4.upperCI[i] = quantile(boot$genic, c(0.975))
  
  df$expected.hotspots.intergenic.filtered4[i] = mean(boot$intergenic)
  df$expected.hotspots.intergenic.filtered4.lowerCI[i] = quantile(boot$intergenic, c(0.025))
  df$expected.hotspots.intergenic.filtered4.upperCI[i] = quantile(boot$intergenic, c(0.975))
}

### Table 3
# SAVE TABLE 3
write.xlsx(x = df[,c(2,11:18)], file = paste("Table/Table3.xls", sep = ""), sheetName = "Hotspot_position", row.names = FALSE, append = FALSE)


write.xlsx(x = df, file = paste("Table/TableS2.xls", sep = ""), sheetName = "Hotspot_position", row.names = FALSE, append = FALSE)
