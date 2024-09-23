#!/usr/bin/env Rscript

source("Source/init.R")

# df_rho = readRDS(file = "Data/Recombination/Gradient/gff_rho_all.rds")

# Simplify, use a pooled GFF
gff_total = readRDS("Data/Genome/gff_all.rds")
gff_total = gff_total[which(gff_total$type == "gene"),]
gff_total$strand[which(!(gff_total$strand == "+") & !(gff_total$strand == "-"))] = "*"
gff_total = makeGRangesFromDataFrame(gff_total,
                                     seqnames.field = "seqid",
                                     keep.extra.columns = TRUE)

snps = read.table(gzfile("Data/Polymorphism/vcf_alleles.gz"),
                  header = TRUE)
snps$Start = snps$Position
snps$End = snps$Position + 1
snps = makeGRangesFromDataFrame(snps,
                                seqnames.field = "Chromosome",
                                keep.extra.columns = TRUE)


# The function to estimate Rho around hotspot centre ----
rho_hotspot_centre = function(range.size = 5000,
                              bin.size = 200,
                              max.length = 10^4,
                              peak.rate = 0,
                              intensity = 2,
                              max.intensity = 10^6,
                              hotspot.type = c("All", "Genic", "Intergenic", "TSS", "TTS")) {
  
  hotspots = read.ldhot.all(max.length = max.length,
                            peak.rate = peak.rate,
                            intensity = intensity,
                            max.intensity = max.intensity)
  
  idx = seq(-range.size, range.size, bin.size)
  
  dist2hotspot = data.frame(idx = idx,
                            dataset = rep(list_dataset, each = length(idx)))
  dist2hotspot$meanRho = NA
  dist2hotspot$meanRho.control = NA
  dist2hotspot$medianRho = NA
  dist2hotspot$medianRho.control = NA
  dist2hotspot$speciesRho = NA
  dist2hotspot$weighted.meanRho = NA
  dist2hotspot$weighted.meanRho.control = NA
  dist2hotspot$weighted.speciesRho = NA
  dist2hotspot$SNPCount = NA
  dist2hotspot$SNPCount.control = NA
  dist2hotspot$speciesSNPCount = NA
  # dist2hotspot$SNPCountAT = NA
  # dist2hotspot$SNPCountAT.control = NA
  # dist2hotspot$SNPCountGC = NA
  # dist2hotspot$SNPCountGC.control = NA
  dist2hotspot$Ts = NA
  dist2hotspot$Ts.control = NA
  dist2hotspot$Tv = NA
  dist2hotspot$Tv.control = NA
  dist2hotspot$geneCount = NA
  dist2hotspot$geneCount.control = NA
  dist2hotspot$geneCount.start = NA
  dist2hotspot$geneCount.start.control = NA
  dist2hotspot$geneCount.end = NA
  dist2hotspot$geneCount.end.control = NA
  dist2hotspot$speciesGeneCount = NA

  hotspots$midpoint = round(hotspots$start + hotspots$length/2, digits = 0)
  
  for (i in 1:nrow(dist2hotspot)) {
    cat(i, ":", round(i/nrow(dist2hotspot) * 100, digits = 0), "%\n")
    # Get hotspot midpoint
    s = dist2hotspot$dataset[i]
    hot = hotspots[which(hotspots$dataset == s),]

    # hot$strand = "*"
    hot$start = hot$midpoint + dist2hotspot$idx[i]
    hot$end = hot$midpoint + dist2hotspot$idx[i] + bin.size - 1
    hot = makeGRangesFromDataFrame(hot, keep.extra.columns=TRUE)
    hot
    
    gff = gff_total[which(gff_total$dataset == s)]
    gff
    
    # Translate chromosome names
    # seqnames(gff) = droplevels(seqnames(gff))
    lev = seqlevels(gff)
    
    for (seq in 1:length(lev)) {
      if (sum(chromosome_metadata$set == s & chromosome_metadata$annotname == lev[seq]) > 0) {
        chromosome_metadata$ldmapname[which(chromosome_metadata$set == s & chromosome_metadata$annotname == lev[seq])]
        lev[seq] = chromosome_metadata$ldmapname[which(chromosome_metadata$set == s & chromosome_metadata$annotname == lev[seq])]
      }
    }
    
    seqlevels(gff) = lev
    
    # Define TSS and TTS regions as TSS/TTS +- 500bp
    startGenes = gff
    forward = which(strand(startGenes) == "+")
    reverse = which(strand(startGenes) == "-")
    end(startGenes[forward]) = start(startGenes[forward]) + 500
    start(startGenes[forward]) = start(startGenes[forward]) - 500
    end(startGenes[reverse]) = end(startGenes[reverse]) + 500
    start(startGenes[reverse]) = end(startGenes[reverse]) - 500
    
    endGenes = gff
    forward = which(strand(endGenes) == "+")
    reverse = which(strand(endGenes) == "-")
    start(endGenes[forward]) = end(endGenes[forward]) - 500
    end(endGenes[forward]) = end(endGenes[forward]) + 500
    start(endGenes[reverse]) = start(endGenes[reverse]) - 500
    end(endGenes[reverse]) = start(endGenes[reverse]) + 500
    
    # Keep or remove only genic/intergenic hotspots
    if (hotspot.type == "Genic") {
      hits = findOverlaps(hot, gff)
      hot = hot[unique(queryHits(hits))]
    }
    if (hotspot.type == "Intergenic") {
      hits = findOverlaps(hot, gff)
      hot = hot[-unique(queryHits(hits))]
    }
    if (hotspot.type == "TSS") {
      hits = findOverlaps(hot, startGenes)
      hot = hot[unique(queryHits(hits))]
    }
    if (hotspot.type == "TTS") {
      hits = findOverlaps(hot, endGenes)
      hot = hot[unique(queryHits(hits))]
    }
    
    # Rho
    # Import Rho maps only if we changed of set
    if (i == 1) {
      cat(s, "\n")
      rho = read.table(file = gzfile(paste("Data/Recombination/LD/ldhat/", dist2hotspot$dataset[i], ".csv.gz", sep = "")), header = TRUE)
      # rho$strand = "*"
      rho = makeGRangesFromDataFrame(rho, keep.extra.columns=TRUE)
    } else {
      if (dist2hotspot$dataset[i] != dist2hotspot$dataset[i - 1]) {
        cat(s, "\n")
        rho = read.table(file = gzfile(paste("Data/Recombination/LD/ldhat/", dist2hotspot$dataset[i], ".csv.gz", sep = "")), header = TRUE)
        # rho$strand = "*"
        rho = makeGRangesFromDataFrame(rho, keep.extra.columns=TRUE)
      }
    }
    
    # Random control
    # Random resample in a nearby region (+- 50 kb)
    resample = sample(c(5000:50000, -50000:-5000), size = length(hot), replace = TRUE)
    hot.control = GRanges(seqnames = as.character(seqnames(hot)),
                          strand = "*",
                          ranges = IRanges(start = start(hot) + resample, width = width(hot)))
    
    # # Random resample anywhere on the chromosome
    # resample = sample(c(5000:50000, -50000:-5000), size = length(hot), replace = TRUE)
    # hot.control = GRanges(seqnames = as.character(seqnames(hot)),
    #                       strand = "*",
    #                       ranges = IRanges(start = start(hot) + resample, width = width(hot)))
    
    dist2hotspot$n.hotspots[i] = length(hot)
    
    # Mean Rho in the window and average Rho of the species
    hits = findOverlaps(hot, rho)
    dist2hotspot$meanRho[i] = mean(rho$Mean_rho[subjectHits(hits)], na.rm = TRUE)
    dist2hotspot$medianRho[i] = median(rho$Mean_rho[subjectHits(hits)], na.rm = TRUE)
    dist2hotspot$weighted.meanRho[i] = weighted.mean(rho$Mean_rho[subjectHits(hits)], width(rho)[subjectHits(hits)], na.rm = TRUE)
    dist2hotspot$SNPCount[i] = length(rho$Mean_rho[subjectHits(hits)]) - 1
    
    # hits = findOverlaps(hot.bg, rho)
    # dist2hotspot$meanRho.bg[i] = mean(rho$Mean_rho[subjectHits(hits)], na.rm = TRUE)
    # dist2hotspot$weighted.meanRho.bg[i] = weighted.mean(rho$Mean_rho[subjectHits(hits)], width(rho)[subjectHits(hits)], na.rm = TRUE)
    # dist2hotspot$SNPCount.bg[i] = length(rho$Mean_rho[subjectHits(hits)]) - 1
    
    hits = findOverlaps(hot.control, rho)
    dist2hotspot$meanRho.control[i] = mean(rho$Mean_rho[subjectHits(hits)], na.rm = TRUE)
    dist2hotspot$medianRho.control[i] = median(rho$Mean_rho[subjectHits(hits)], na.rm = TRUE)
    dist2hotspot$weighted.meanRho.control[i] = weighted.mean(rho$Mean_rho[subjectHits(hits)], width(rho)[subjectHits(hits)], na.rm = TRUE)
    dist2hotspot$SNPCount.control[i] = length(rho$Mean_rho[subjectHits(hits)]) - 1
    
    dist2hotspot$speciesRho[i] = mean(rho$Mean_rho, na.rm = TRUE)
    dist2hotspot$speciesMedianRho[i] = median(rho$Mean_rho, na.rm = TRUE)
    dist2hotspot$weighted.speciesRho[i] = weighted.mean(rho$Mean_rho, width(rho), na.rm = TRUE)
    dist2hotspot$speciesSNPCount[i] = length(rho$Mean_rho) - 1
    
    # Count AT vs GC SNPs
    # AT.snps = snps[which(snps$Dataset == s)]
    # AT.snps = AT.snps[which(AT.snps$Ref %in% c("A", "T") & AT.snps$Alt %in% c("A", "T"))]
    # GC.snps = snps[which(snps$Dataset == s)]
    # GC.snps = GC.snps[which(GC.snps$Ref %in% c("G", "C") & GC.snps$Alt %in% c("G", "C"))]
    # 
    # hits = findOverlaps(hot, AT.snps)
    # dist2hotspot$SNPCountAT[i] = length(unique(subjectHits(hits)))
    # 
    # hits = findOverlaps(hot.control, AT.snps)
    # dist2hotspot$SNPCountAT.control[i] = length(unique(subjectHits(hits)))
    # 
    # hits = findOverlaps(hot, GC.snps)
    # dist2hotspot$SNPCountGC[i] = length(unique(subjectHits(hits)))
    # 
    # hits = findOverlaps(hot.control, GC.snps)
    # dist2hotspot$SNPCountGC.control[i] = length(unique(subjectHits(hits)))
    
    # Count Ts/Tv
    Ts.snps = snps[which(snps$Dataset == s)]
    Ts.idx = which((Ts.snps$Ref == "A" & Ts.snps$Alt == "G") | (Ts.snps$Ref == "G" & Ts.snps$Alt == "A") | (Ts.snps$Ref == "C" & Ts.snps$Alt == "T") | (Ts.snps$Ref == "T" & Ts.snps$Alt == "C"))
    Ts.snps = Ts.snps[Ts.idx]
    
    hits = findOverlaps(hot, Ts.snps)
    dist2hotspot$Ts[i] = length(unique(subjectHits(hits)))
    dist2hotspot$Tv[i] = dist2hotspot$SNPCount[i] - dist2hotspot$Ts[i]
    
    hits = findOverlaps(hot.control, Ts.snps)
    dist2hotspot$Ts.control[i] = length(unique(subjectHits(hits)))
    dist2hotspot$Tv.control[i] = dist2hotspot$SNPCount.control[i] - dist2hotspot$Ts.control[i]
    
    
    # Count genes overlapping the window
    hits = findOverlaps(hot, gff)
    dist2hotspot$geneCount[i] = length(unique(subjectHits(hits)))
    hits = findOverlaps(hot.control, gff)
    dist2hotspot$geneCount.control[i] = length(unique(subjectHits(hits)))
    
    # startGenes = gff
    # forward = which(strand(startGenes) == "+")
    # reverse = which(strand(startGenes) == "-")
    # end(startGenes[forward]) = start(startGenes[forward]) + 500
    # start(startGenes[forward]) = start(startGenes[forward]) - 500
    # end(startGenes[reverse]) = end(startGenes[reverse]) + 500
    # start(startGenes[reverse]) = end(startGenes[reverse]) - 500
    hits = findOverlaps(hot, startGenes)
    dist2hotspot$geneCount.start[i] = length(unique(subjectHits(hits)))
    hits = findOverlaps(hot.control, startGenes)
    dist2hotspot$geneCount.start.control[i] = length(unique(subjectHits(hits)))
    
    # endGenes = gff
    # forward = which(strand(endGenes) == "+")
    # reverse = which(strand(endGenes) == "-")
    # start(endGenes[forward]) = end(endGenes[forward]) - 500
    # end(endGenes[forward]) = end(endGenes[forward]) + 500
    # start(endGenes[reverse]) = start(endGenes[reverse]) - 500
    # end(endGenes[reverse]) = start(endGenes[reverse]) + 500
    hits = findOverlaps(hot, endGenes)
    dist2hotspot$geneCount.end[i] = length(unique(subjectHits(hits)))
    hits = findOverlaps(hot.control, endGenes)
    dist2hotspot$geneCount.end.control[i] = length(unique(subjectHits(hits)))
    
    dist2hotspot$speciesGeneCount[i] = length(gff)
  }
  
  filename = paste0("Data/Recombination/rho_hotspot_center.size",
                     as.integer(max.length),".intensity", intensity,
                     ".maxintensity", as.integer(max.intensity), "_", hotspot.type, ".rds")
  saveRDS(dist2hotspot, file = filename)
  
  rm(dist2hotspot)
}


#========================================================================#
# Estimate ----

for (hot_type in c("All", "Genic", "Intergenic", "TSS", "TTS")) {
  cat("======================\nProcessing", hot_type, "...\n")
  # Raw
  rho_hotspot_centre(max.length = 10^8,
                     peak.rate = 0,
                     intensity = 0,
                     max.intensity = 10^8,
                     hotspot.type = hot_type)

  # Filter by size only
  rho_hotspot_centre(max.length = 50000,
                     peak.rate = 0,
                     intensity = 0,
                     max.intensity = 10^8,
                     hotspot.type = hot_type)
  
  rho_hotspot_centre(max.length = 10^4,
                     peak.rate = 0,
                     intensity = 0,
                     max.intensity = 10^8,
                     hotspot.type = hot_type)

  # Filter by intensity only
  rho_hotspot_centre(max.length = 10^8,
                     peak.rate = 0,
                     intensity = 4,
                     max.intensity = 10^8,
                     hotspot.type = hot_type)
  
  # Filtered
  rho_hotspot_centre(max.length = 10^4,
                     peak.rate = 0,
                     intensity = 2,
                     max.intensity = 10^8,
                     hotspot.type = hot_type)

  rho_hotspot_centre(max.length = 10^4,
                     peak.rate = 0,
                     intensity = 4,
                     max.intensity = 10^8,
                     hotspot.type = hot_type)

  # Filtered by max intensity too

  rho_hotspot_centre(max.length = 10^4,
                     peak.rate = 0,
                     intensity = 2,
                     max.intensity = 400,
                     hotspot.type = hot_type)

  rho_hotspot_centre(max.length = 10^4,
                     peak.rate = 0,
                     intensity = 4,
                     max.intensity = 200,
                     hotspot.type = hot_type)
}
