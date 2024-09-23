########################################################################## #
#                    ECOBIO - PhD
#
# 
#         Generate LD maps figures ----
#
########################################################################## #
rm(list=ls(all=TRUE))
# args = commandArgs(trailingOnly=TRUE)
# wd = args[1]
# wd = dirname(rstudioapi::getSourceEditorContext()$path)
# wd = gsub("/Jobs", "", wd)
# setwd(wd)

source("Source/init.R")


#==================================================================#
# - [x] LDmap landscape (Rho/kb) ----
# s = "Arabidopsis_thaliana_1001genomes"
# chr = "1"
cat("LDmap landscape (Rho/kb)...\n")

source("Source/read.ldmap.R")

for (s in list_dataset) {
  cat(s, "\n")
  list_chr = chromosome_metadata$ldmapname[which(chromosome_metadata$set == s)]
  for (chr in list_chr) {
    cat(chr, "\n")
    if (file.exists(paste0("Data/Recombination/LD/ldhat/", s, ".", chr, ".bpen", bpen, ".res.txt.gz"))) {
      ldmap = read.ldmap(s, chr, 5)
      
      # ggplot(ldmap, aes(x = end/10^6, y = Mean_rho)) +
      #   geom_line() +
      #   xlab("Genomic position (Mb)") + ylab("Mean Rho/kb") +
      #   ggtitle(paste(s, "chromosome", chr))
      
      ggplot(data = ldmap, aes(x = end/10^6, y = Mean_rho_unfiltered)) +
        geom_line(alpha = 0.4) +
        # geom_step(aes(x = start/10^6, y = Mean_rho), alpha = 0.3) +
        geom_segment(data = ldmap[which(is.na(ldmap$Mean_rho)),], aes(x = start/10^6, xend = end/10^6, y = Mean_rho_unfiltered, yend = Mean_rho_unfiltered), color = "Red") +
        xlab("Genomic position (Mb)") + ylab("Mean Rho/kb") +
        ggtitle(paste(s, "chromosome", chr))

      ggsave(paste0("Figure/Landscapes/Rho/", s, ".", chr, ".jpeg"),
             width = 22, height = 15, units = "cm")
      
      rm(ldmap)
      gc()
    }
  }
}


#==================================================================#
# - [x] Rho in 100k windows ----
# s = "Arabidopsis_thaliana_1001genomes"
# chr = "1"

source("Source/read.ldmap.R")

for (s in list_dataset) {
  cat(s, "\n")
  list_chr = chromosome_metadata$ldmapname[which(chromosome_metadata$set == s)]
  for (chr in list_chr) {
    cat(chr, "\n")
    if (file.exists(paste0("Data/Recombination/LD/ldhat/", s, ".", chr, ".bpen", bpen, ".res.txt.gz"))) {
      ldmap = read.ldmap(s, chr, 5)
      ldmap$width = ldmap$end - ldmap$start
      
      intervals = data.frame(start = seq(1, max(ldmap$end), 100*10^3),
                             end = seq(1, max(ldmap$end), 100*10^3) + 100*10^3 -1,
                             rho = NA)
      intervals$midpoint = round((intervals$start + intervals$end)/2, digits = 2)
      
      for (i in 1:nrow(intervals)) {
        tmp = subset(ldmap, ldmap$Loci >= intervals$start[i] & ldmap$Loci <= intervals$end[i])
        intervals$rho[i] = weighted.mean(tmp$Mean_rho, tmp$width)
      }
      
      ggplot(intervals, aes(x = midpoint/10^6, y = rho)) +
        geom_line() +
        xlab("Genomic position (Mb)") + ylab("Mean Rho/kb (100kb windows)") +
        ggtitle(paste(s, "chromosome", chr))
      
      ggsave(paste0("Figure/Landscapes/Rho_sliding100kb/", s, ".", chr, ".jpeg"),
             width = 22, height = 15, units = "cm")
      rm(tmp)
      rm(intervals)
      rm(ldmap)
      gc()
    }
  }
}


#==================================================================#
# - [x] LDhot - hotspots peak rates ----
# s = "Arabidopsis_thaliana_1001genomes"
# chr = "1"
cat("LDhot hotspots peak rate...\n")

# if (file.exists(paste0("Data/Recombination/LD/ldhotspots_filtered.csv.gz"))) {
# ldhot = read.table(file = gzfile(paste0("Data/Recombination/LD/ldhotspots_filtered.csv.gz")), header = T)
ldhot = read.ldhot.all(max.length = 10^4,
                               peak.rate = 0,
                               intensity = 4,
                               max.intensity = 200)

# ldhot = subset(ldhot, ldhot$dataset %in% list_dataset)

ldhot$dataset = as.factor(ldhot$dataset)
if ("width" %in% colnames(ldhot) & !("length" %in% colnames(ldhot))) {
  ldhot$length = ldhot$width
}

ldhot$midpoint = (ldhot$end + ldhot$start)/2

for (s in list_dataset) {
  cat(s, "\n")
  hotspots = subset(ldhot, ldhot$dataset == s)
  list_chr = chromosome_metadata$ldmapname[which(chromosome_metadata$set == s)]
  for (chr in list_chr) {
    cat(chr, "\n")
    hot_chrom = subset(hotspots, hotspots$seqnames == chr)
    
    if (nrow(hot_chrom) > 0) {
      ggplot(hot_chrom, aes(x = midpoint/10^6, y = peak_rate)) +
        geom_line() +
        xlab("Genomic position (Mb)") + ylab("Hotspot peak rate") +
        ggtitle(paste(s, "chromosome", chr))
      
      ggsave(paste0("Figure/Landscapes/Hotspot_peak_rate/", s, ".", chr, ".jpeg"),
             width = 22, height = 15, units = "cm")
      
      ggplot(hot_chrom, aes(x = midpoint/10^6, y = rho_across_hotspot)) +
        geom_line() +
        xlab("Genomic position (Mb)") + ylab("Rho across hotspot") +
        ggtitle(paste(s, "chromosome", chr))
      
      ggsave(paste0("Figure/Landscapes/Hotspot_rho/", s, ".", chr, ".jpeg"),
             width = 22, height = 15, units = "cm")
      
      
      ggplot(hot_chrom, aes(x = midpoint/10^6, y = length/10^3)) +
        geom_line() +
        xlab("Genomic position (Mb)") + ylab("Hotspot size (kb)") +
        ggtitle(paste(s, "chromosome", chr))
      
      ggsave(paste0("Figure/Landscapes/Hotspot_size/", s, ".", chr, ".jpeg"),
             width = 22, height = 15, units = "cm")
      
      rm(hot_chrom)
    }
  }
  rm(hotspots)
}

gc()
# }

#==================================================================#
# - [x] Number of hotspots in 100kb windows ----
# s = "Arabidopsis_thaliana_1001genomes"
# chr = "1"
cat("Number of hotspots in 100kb windows...\n")

# if (file.exists(paste0("Data/Recombination/LD/ldhotspots_filtered.csv.gz"))) {
# ldhot = read.table(file = gzfile(paste0("Data/Recombination/LD/ldhotspots_filtered.csv.gz")), header = T)
# ldhot = subset(ldhot, ldhot$dataset %in% list_dataset)
# ldhot$dataset = as.factor(ldhot$dataset)
# ldhot$length = ldhot$end - ldhot$start

for (s in list_dataset) {
  cat(s, "\n")
  hotspots = subset(ldhot, ldhot$dataset == s)
  list_chr = chromosome_metadata$ldmapname[which(chromosome_metadata$set == s)]
  for (chr in list_chr) {
    cat(chr, "\n")
    hot_chrom = subset(hotspots, hotspots$seqnames == chr)
    
    if (nrow(hot_chrom) > 0) {
      hot = GRanges(seqnames = chr,
                    ranges = IRanges(start = hot_chrom$start, end = hot_chrom$end),
                    strand = "*")
      
      intervals = GRanges(seqnames = chr,
                          ranges = IRanges(start = seq(1, max(end(hot)), 10*10^4), width = 10*10^4),
                          strand = "*")
      
      hits = countOverlaps(intervals, hot)
      
      intervals$hotspot_count = hits
      
      intervals = as.data.frame(intervals)
      intervals$midpoint = (intervals$start + intervals$end)/2
      ggplot(intervals, aes(x = midpoint/10^6, y = hotspot_count)) +
        geom_line() +
        xlab("Genomic position (Mb)") + ylab("Number of hotspots overlapping 100 kb windows") +
        ggtitle(paste(s, "chromosome", chr))
      
      ggsave(paste0("Figure/Landscapes/Hotspot_count/", s, ".", chr, ".jpeg"),
             width = 22, height = 15, units = "cm")
      
      rm(hot)
      rm(intervals)
      rm(hits)
    }
    rm(hot_chrom)
  }
  rm(hotspots)
}

rm(ldhot)
gc()
# }


cat("End of script")
