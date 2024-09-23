read.ldmap = function(set, chr, bpen = 5) {
  ldmap = data.frame(set = character(0), chromosome = character(0), Loci = character(0), Mean_rho = numeric(0), Median = numeric(0), L95 = numeric(0), U95 = numeric(0))
  for (i in 1:length(chr)) {
    chromosome = chr[i]
    mapfile = paste("Data/Recombination/LD/ldhat/", set, ".", chromosome, ".bpen", bpen,".res.txt.gz", sep = "")
    if (file.exists(mapfile)) {
      map = read.table(gzfile(mapfile), header = T)
      # Put positions in bp
      map$Loci = map$Loci * 1000
      # Remove negative position at start
      map = map[-which(map$Loci < 0),]
      # Map ranges
      map$start = map$Loci
      map$end = c((map$Loci[-1] - 1), 0)
      # map$start = c(0, head(map$Loci, -1))
      # map$end = (map$Loci - 1)
      map = head(map, -1)
      
      ######### Filter segments with constant rho over > 100kb
      # Save unfiltered Rho and keep compatibility with downstream analyses
      map$Mean_rho_unfiltered = map$Mean_rho
      
      # Get segments
      library(dostats)
      
      # Detect groups
      grp = seq_consecutive(map$Mean_rho)
      map$group = grp
      
      # Get the genomic size of each group
      grp_start = aggregate(start ~ grp, map, min)
      grp_end = aggregate(end ~ grp, map, max)
      grp_end$length = grp_end$end - grp_start$start
      
      # Size filter: remove segments/groups larger than this
      size_filter = 10^5
      idx = grp_end$grp[which(grp_end$length > size_filter)]
      idx
      # Set Rho to NA in these groupes
      map$Mean_rho[which(map$group %in% idx)] = NA
      
      # Combine
      map = cbind(set, chromosome, map)
      ldmap = rbind(ldmap, map)
      return(ldmap)
    } else {
      return(invisible(NULL))
    }
  }
}

# TESTING
# ggplot(data = map, aes(x = end/10^6, y = Mean_rho)) +
#   geom_step(aes(x = start/10^6, y = Mean_rho), alpha = 0.3) +
#   geom_segment(data = map[which(is.na(map$Mean_rho)),], aes(x = start/10^6, xend = end/10^6, y = Mean_rho_unfiltered, yend = Mean_rho_unfiltered), color = "Red") +
#   xlab("Genomic position (Mb)") + ylab("Mean Rho/kb") +
#   ggtitle(paste(set, "chromosome", chromosome))




plot.ldmap = function(dataset, chromosome, start, end, gff = data.frame()) {
  ldmap = read.ldmap(dataset, chromosome)
  ldmap = ldmap[which(ldmap$Loci >= start & ldmap$Loci <= end),]
  library(ggplot2)
  p = ggplot(data = ldmap, aes(x = end/10^6, y = Mean_rho)) +
    geom_step(aes(x = start/10^6, y = Mean_rho)) +
    xlab("Genomic position (Mb)") + ylab("Mean Rho/kb") +
    ggtitle(paste(dataset, "chromosome", chromosome))
  
  # Add hotspots in red
  ldhot = read.ldhot(dataset, chromosome, max.length = 10000, peak.rate = 3)
  ldhot = ldhot[which((ldhot$start >= start & ldhot$start <= end) | (ldhot$end >= start & ldhot$end <= end)),]
  if (nrow(ldhot) > 0) {
    ldhot$centre = (ldhot$start + ldhot$end)/2
    ldhot$zero = 0
    
    p = p + geom_vline(xintercept = ldhot$centre/10^6, color = "Red") +
      geom_segment(data = ldhot, aes(x = start/10^6, xend = end/10^6, y = zero, yend = zero),
                   lineend = "round",
                   color = "red")
                   # arrow = arrow(ends = "both", length = unit(0.1, "inches")))
  }
  
  
  # Add genes in blue
  if (nrow(gff) > 0) {
    gff = subset(gff, gff$feature == "gene")
    gff = subset(gff, gff$set == dataset)
    annot_name = chromosome_metadata$annotname[which(chromosome_metadata$set == dataset & chromosome_metadata$ldmapname == chromosome)]
    gff = subset(gff, gff$seqnames == annot_name)
    gff = gff[which((gff$start >= start & gff$start <= end) | (gff$end >= start & gff$end <= end)),]
    gff$y = max(ldmap$Mean_rho) + 0.3
    p = p + geom_segment(data = gff, aes(x = start/10^6, xend = end/10^6, y = y, yend = y),
                   lineend = "butt",
                   color = "blue")
  }
  
  
  return(p)
}

