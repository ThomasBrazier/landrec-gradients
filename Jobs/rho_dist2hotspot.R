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

# The function to estimate Rho around hotspot centre ----
rho_hotspot_centre = function(range.size = 5000,
                              bin.size = 200,
                              max.length = 10^4,
                              peak.rate = 0,
                              intensity = 2,
                              max.intensity = 10^6) {
  hotspots = read.ldhot.all(max.length = max.length,
                            peak.rate = peak.rate,
                            intensity = intensity,
                            max.intensity = max.intensity)
  
  interval = bin.size
  bounds = range.size
  idx = seq(-bounds, bounds, interval)
  
  dist2hotspot = data.frame(idx = idx,
                            dataset = rep(list_dataset, each = length(idx)))
  dist2hotspot$meanRho = NA
  dist2hotspot$SNPCount = NA
  dist2hotspot$meanRho.control = NA
  dist2hotspot$SNPCount.control = NA
  dist2hotspot$speciesRho = NA
  hotspots$midpoint = round(hotspots$start + hotspots$length/2, digits = 0)
  dist2hotspot$geneCount = NA
  dist2hotspot$geneCount.control = NA
  dist2hotspot$geneCount.start = NA
  dist2hotspot$geneCount.start.control = NA
  dist2hotspot$geneCount.end = NA
  dist2hotspot$geneCount.end.control = NA
  dist2hotspot$speciesGeneCount = NA
  
  for (i in 1:nrow(dist2hotspot)) {
    cat(i, ":", round(i/nrow(dist2hotspot) * 100, digits = 0), "%\n")
    # Get hotspot midpoint
    s = dist2hotspot$dataset[i]
    hot = hotspots[which(hotspots$dataset == s),]
    
    # hot$strand = "*"
    hot$start = hot$midpoint + dist2hotspot$idx[i]
    hot$end = hot$midpoint + dist2hotspot$idx[i] + interval - 1
    hot = makeGRangesFromDataFrame(hot, keep.extra.columns=TRUE)
    hot
    
    gff = gff_total[which(gff_total$dataset == s)]
    # Translate chromosome names
    # seqnames(gff) = droplevels(seqnames(gff))
    lev = seqlevels(gff)
    
    for (seq in 1:length(lev)) {
      if (sum(chromosome_metadata$set == s & chromosome_metadata$annotname == lev[seq]) > 0) {
        lev[seq] = chromosome_metadata$ldmapname[which(chromosome_metadata$set == s & chromosome_metadata$annotname == lev[seq])]
      }
    }
    
    seqlevels(gff) = lev
    
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
    # Random resample in a nearby region (+- 10 kb)
    resample = sample(c(5000:15000, -15000:-5000), size = length(hot), replace = TRUE)
    hot.control = GRanges(seqnames = as.character(seqnames(hot)),
                          strand = "*",
                          ranges = IRanges(start = start(hot) + resample, width = width(hot)))
    
    dist2hotspot$n.hotspots[i] = length(hot)
    
    # Mean Rho in the window and average Rho of the species
    hits = findOverlaps(hot, rho)
    dist2hotspot$meanRho[i] = mean(rho$Mean_rho[subjectHits(hits)], na.rm = TRUE)
    dist2hotspot$SNPCount[i] = length(rho$Mean_rho[subjectHits(hits)]) - 1
    
    hits = findOverlaps(hot.control, rho)
    dist2hotspot$meanRho.control[i] = mean(rho$Mean_rho[subjectHits(hits)], na.rm = TRUE)
    dist2hotspot$SNPCount.control[i] = length(rho$Mean_rho[subjectHits(hits)]) - 1
    
    dist2hotspot$speciesRho[i] = mean(rho$Mean_rho, na.rm = TRUE)
    dist2hotspot$speciesSNPCount[i] = length(rho$Mean_rho, na.rm = TRUE)
    
    # Count genes overlapping the window
    hits = findOverlaps(hot, gff)
    dist2hotspot$geneCount[i] = length(unique(subjectHits(hits)))
    hits = findOverlaps(hot.control, gff)
    dist2hotspot$geneCount.control[i] = length(unique(subjectHits(hits)))
    
    startGenes = gff
    forward = which(strand(startGenes) == "+")
    reverse = which(strand(startGenes) == "-")
    end(startGenes[forward]) = start(startGenes[forward]) + 500
    start(startGenes[forward]) = start(startGenes[forward]) - 500
    end(startGenes[reverse]) = end(startGenes[reverse]) + 500
    start(startGenes[reverse]) = end(startGenes[reverse]) - 500
    hits = findOverlaps(hot, startGenes)
    dist2hotspot$geneCount.start[i] = length(unique(subjectHits(hits)))
    hits = findOverlaps(hot.control, startGenes)
    dist2hotspot$geneCount.start.control[i] = length(unique(subjectHits(hits)))
    
    endGenes = gff
    forward = which(strand(endGenes) == "+")
    reverse = which(strand(endGenes) == "-")
    start(endGenes[forward]) = end(endGenes[forward]) - 500
    end(endGenes[forward]) = end(endGenes[forward]) + 500
    start(endGenes[reverse]) = start(endGenes[reverse]) - 500
    end(endGenes[reverse]) = start(endGenes[reverse]) + 500
    hits = findOverlaps(hot, endGenes)
    dist2hotspot$geneCount.end[i] = length(unique(subjectHits(hits)))
    hits = findOverlaps(hot.control, endGenes)
    dist2hotspot$geneCount.end.control[i] = length(unique(subjectHits(hits)))
    
    dist2hotspot$speciesGeneCount[i] = length(gff)
  }
  
  filename = paste0("Data/Recombination/dist2hotspot_size",
                    max.length,"_intensity", intensity,
                    "_maxintensity", max.intensity, ".rds")
  saveRDS(dist2hotspot, file = filename)
  
  rm(dist2hotspot)
}


#========================================================================#
# Estimate ----

# Raw
rho_hotspot_centre(max.length = 10^6,
                   peak.rate = 0,
                   intensity = 0,
                   max.intensity = 10^8)

rho_hotspot_centre(max.length = 50000,
                   peak.rate = 0,
                   intensity = 0,
                   max.intensity = 10^8)

# Filter by size only
rho_hotspot_centre(max.length = 10^4,
                   peak.rate = 0,
                   intensity = 0,
                   max.intensity = 10^8)

# Filtered
rho_hotspot_centre(max.length = 10^4,
                   peak.rate = 0,
                   intensity = 2,
                   max.intensity = 10^8)

rho_hotspot_centre(max.length = 10^4,
                   peak.rate = 0,
                   intensity = 4,
                   max.intensity = 10^8)

# Filtered by max intensity too

rho_hotspot_centre(max.length = 10^4,
                   peak.rate = 0,
                   intensity = 2,
                   max.intensity = 400)

rho_hotspot_centre(max.length = 10^4,
                   peak.rate = 0,
                   intensity = 4,
                   max.intensity = 400)