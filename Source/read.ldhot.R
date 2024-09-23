read.ldhot = function(dataset,
                      chromosome,
                      bpen = 5,
                      max.length = 10000,
                      peak.rate = 0,
                      detailed.results = T,
                      intensity = 2,
                      max.intensity = 10^8) {
  require(GenomicRanges)
  require(pbmcapply)
  cat("Read LDhot for", dataset, "chromosome", chromosome, "\n")
  ldhot = data.frame(seqnames = character(0), start = numeric(0), end = numeric(0),	p = numeric(0),	rho_across_hotspot = numeric(0),	peak_rate = numeric(0))
  hotfile = paste("Data/Recombination/LD/ldhot/", dataset, ".", chromosome, ".bpen", bpen,".hot_summary.txt.gz", sep = "")
  mapfile = paste("Data/Recombination/LD/ldhat/", dataset, ".", chromosome, ".bpen", bpen,".res.txt.gz", sep = "")
  if (file.exists(mapfile)) {
    if (file.exists(hotfile)) {
      hot = read.table(gzfile(hotfile), header = F)
      hot = cbind(chromosome, hot)
      colnames(hot) = colnames(ldhot)
      hot = cbind(dataset, hot)
      # Distances are converted in bp
      hot$start = hot$start*10^3
      hot$end = hot$end*10^3
      
      if (detailed.results == TRUE) {
        
        map = read.ldmap(dataset, chromosome, bpen = 5)
        
        map = makeGRangesFromDataFrame(map, keep.extra.columns = TRUE)
        
        estimated_bg_rate = function(idx) {
          hot_summary = makeGRangesFromDataFrame(hot[idx,],
                                                 keep.extra.columns = TRUE)
          # Extend region to background +- 50kb
          hot_summary
          start(hot_summary) = start(hot_summary) - 50*10^3
          end(hot_summary) = end(hot_summary) + 50*10^3
          
          overlaps = findOverlaps(hot_summary, map)
          
          me = map$Mean_rho[subjectHits(overlaps)]
          we = width(map)[subjectHits(overlaps)]
          mean_estimated_bg_rate = weighted.mean(me,
                                                 we,
                                                 na.rm = TRUE)
          return(mean_estimated_bg_rate)
        }
        
        estimated_mean_rho = function(idx) {
          hot_summary = makeGRangesFromDataFrame(hot[idx,],
                                                 keep.extra.columns = TRUE)

          overlaps = findOverlaps(hot_summary, map)
          
          hot_overlap = map[subjectHits(overlaps)]
          
          hot_overlap = restrict(hot_overlap,
                                 start= start(hot_summary),
                                 end = end(hot_summary))
          
          me = hot_overlap$Mean_rho
          we = width(hot_overlap)
          mean_rho = weighted.mean(me,
                                   we,
                                   na.rm = TRUE)
          return(mean_rho)
        }
        
        hot$estimated_bg_rate = unlist(pbmclapply(1:nrow(hot), estimated_bg_rate))
        hot$estimated_mean_rho = unlist(pbmclapply(1:nrow(hot), estimated_mean_rho))
        hot$intensity = hot$peak_rate/hot$estimated_bg_rate
        hot$intensity.meanrho = hot$estimated_mean_rho/hot$estimated_bg_rate
        
        hotspotfile = paste("Data/Recombination/LD/ldhot/", dataset, ".", chromosome, ".bpen", bpen,".hotspots.txt.gz", sep = "")
        hotspots = read.table(gzfile(hotspotfile),
                              header = T,
                              comment.char = "",
                              sep = "\t",
                              fill = TRUE)
        hotspots = hotspots[which(!is.na(hotspots$HotStart)),]
        hotspots = hotspots[which(!is.na(hotspots$HotEnd)),]
        
        hotspots$start = hotspots$HotStart * 10^3
        hotspots$end = hotspots$HotEnd * 10^3
        hotspots$chromosome = chromosome
        hotspots = makeGRangesFromDataFrame(hotspots,
                                            keep.extra.columns = TRUE)
        
        MLE_bg_rate = function(idx) {
          hot_summary = makeGRangesFromDataFrame(hot[idx,],
                                                 keep.extra.columns = TRUE)
          overlaps = findOverlaps(hot_summary, hotspots)
          mean_MLE_bg_rate = mean(hotspots$MLE_bg_rate[subjectHits(overlaps)])
          return(mean_MLE_bg_rate)
        }
        
        MLE_hotspot_rate = function(idx) {
          hot_summary = makeGRangesFromDataFrame(hot[idx,],
                                                 keep.extra.columns = TRUE)
          overlaps = findOverlaps(hot_summary, hotspots)
          mean_MLE_hotspot_rate = mean(hotspots$MLE_hotspot_rate[subjectHits(overlaps)])
          return(mean_MLE_hotspot_rate)
        }
        
        hot$mean_MLE_bg_rate = unlist(pbmclapply(1:nrow(hot), MLE_bg_rate))
        hot$mean_MLE_hotspot_rate = unlist(pbmclapply(1:nrow(hot), MLE_hotspot_rate))
        
        hot$MLE_intensity = hot$mean_MLE_hotspot_rate/hot$mean_MLE_bg_rate
        
        hot = hotspot_filter(hot, max.length, peak.rate, intensity, max.intensity)
      } else {
        hot = hotspot_filter(hot, max.length, peak.rate)
      }
      
      return(hot)
    } else {
      cat("File", hotfile, "does not exists\n")
      return(invisible(NULL))
    }
  } else {
    cat("File", mapfile, "does not exists\n")
    return(invisible(NULL))
  }
}


read.ldhot.all = function(max.length = 10000,
                          peak.rate = 0,
                          intensity = 2,
                          max.intensity = 10^8) {
  ldhot = read.table(file = gzfile("Data/Recombination/LD/ldhotspots_raw.csv.gz"), header = T)
  
  ldhot$dataset = as.factor(ldhot$dataset)
  ldhot = hotspot_filter(ldhot, max.length, peak.rate, intensity, max.intensity)
  
  return(ldhot)
}


hotspot_filter = function(hot,
                          max.length,
                          peak.rate,
                          intensity,
                          max.intensity) {
  # Hotspot filtering
  # (1) by size
  hot$length = hot$end - hot$start + 1
  hot = subset(hot, hot$length <= max.length)
  # (2) by peak rate
  hot = hot[which(hot$peak_rate >= peak.rate),]
  
  # (3) by intensity
  if ("intensity" %in% colnames(hot)) {
    cat(paste("Filter by intensity >", intensity, "and max intensity <", max.intensity, "\n"))
    hot = hot[which(hot$intensity >= intensity & hot$intensity <= max.intensity),]
  } else {
    cat("Cannot filter by intensity\n")
  }
  
  return(hot)
}
