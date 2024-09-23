#============================================================================#
# Loading variables & objects ----
#============================================================================#
# Get the directory of the file & set working directory
# wd=dirname(rstudioapi::getSourceEditorContext()$path)
# wd=gsub("/Source", "", wd)
# setwd(wd)
# setwd("~/Academic/landrec-gradients")
source("Source/init.R")
# ncpus = 8
# ncores = min(ncpus, detectCores())


source("Source/read.ldhot.R")


#============================================================================#
# Hotspot overlap ----
#============================================================================#
cat("Import GFF\n")
for (s in sets) {
  cat(s, "\n")
  gff_recombination = readRDS(paste("Data/Recombination/Gradient/gff_rho_", s, ".rds", sep = ""))
  cat("Import LD hotspots\n")
  # Import LD hotspots
  # Raw data
  ldhot = read.table(file = gzfile("Data/Recombination/LD/ldhotspots_raw.csv.gz"), header = T)
  ldhot = ldhot[which(ldhot$set == s),]
  ldhot = makeGRangesFromDataFrame(ldhot,
                                   keep.extra.columns=TRUE)
  # Translate chromosome names
  chrnames = seqlevels(ldhot)
  for (j in 1:length(chrnames)) {
    # cat(chr, '\n')
    chrnames[j] = chromosome_metadata$annotname[which(chromosome_metadata$ldmapname == chrnames[j] & chromosome_metadata$set == s)]
  }
  seqlevels(ldhot) = chrnames

  # Trimmed
  ldhot_6kb = read.table(file = gzfile("Data/Recombination/LD/ldhotspots_6kb.csv.gz"), header = T)
  # ldhot_3kb = read.table(file = gzfile("Data/Recombination/LD/ldhotspots_3kb.csv.gz"), header = T)
  ldhot_6kb = ldhot_6kb[which(ldhot_6kb$set == s),]
  # ldhot_3kb = ldhot_3kb[which(ldhot_3kb$set == s),]
  ldhot_6kb = makeGRangesFromDataFrame(ldhot_6kb,
                                       keep.extra.columns=TRUE)
  # ldhot_3kb = makeGRangesFromDataFrame(ldhot_3kb,
  #                                      keep.extra.columns=TRUE)
  seqlevels(ldhot_6kb) = chrnames
  # seqlevels(ldhot_3kb) = chrnames
  
  # GFF ranges ----
  gff_ranges = makeGRangesFromDataFrame(gff_recombination,
                                        keep.extra.columns=TRUE)
  
  # Assess if feature overlaps hotspot/hotspot position ----
  cat("Hotspot overlap\n")
  
  hotoverlap = countOverlaps(gff_ranges, ldhot, type = "any")
  
  if (length(hotoverlap) != length(gff_ranges)) {
    warning("Parsing hotspot overlap introduced missing rows")
  }
  gff_ranges$hotspot_overlap = hotoverlap
  
  # Trimmed hotspots
  # 6kb
  ldhot_6kb = makeGRangesFromDataFrame(ldhot_6kb,
                                              keep.extra.columns=TRUE)
  hotoverlap = countOverlaps(gff_ranges, ldhot_6kb, type = "any")
  
  if (length(hotoverlap) != length(gff_ranges)) {
    warning("Parsing hotspot overlap introduced missing rows")
  }
  gff_ranges$hotspot_overlap_6kb = hotoverlap
  
  
  # 3kb
  # ldhot_3kb = makeGRangesFromDataFrame(ldhot_3kb,
  #                                             keep.extra.columns=TRUE)
  # hotoverlap = countOverlaps(gff_ranges, ldhot_3kb, type = "any")
  # 
  # if (length(hotoverlap) != length(gff_ranges)) {
  #   warning("Parsing hotspot overlap introduced missing rows")
  # }
  # gff_ranges$hotspot_overlap_3kb = hotoverlap

  
  # 4b. Hotspot coverage - How much of the sequence is covered by a hostpot
  # Size of intersect
  hits = findOverlaps(gff_ranges, ldhot, ignore.strand = TRUE)
  grl = extractList(ldhot, as(hits, "List"))
  hotcoverage = width(pintersect(gff_ranges, grl, ignore.strand = TRUE))
  hotcoverage[which(unlist(lapply(hotcoverage, length)) == 0)] = 0
  hotcoverage = lapply(hotcoverage, function(x){mean(unlist(x))})
  hotcoverage = unlist(hotcoverage)
  
  
  if (length(hotcoverage) != length(gff_ranges)) {
    warning("Parsing hotspot coverage introduced missing rows")
  }
  gff_ranges$hotspot_coverage_bp = hotcoverage

  # hotoverlap = countOverlaps(gff_ranges, ldhot, type = "any")
  # if (length(hotoverlap) != length(gff_ranges)) {
  #   warning("Parsing hotspot overlap introduced missing rows")
  # }
  # gff_ranges$hotspot_overlap = hotoverlap
  # table(gff_ranges$hotspot_overlap)
  # 
  # # Trimmed hotspots
  # # 6kb
  # hotoverlap = countOverlaps(gff_ranges, ldhot_6kb, type = "any")
  # if (length(hotoverlap) != length(gff_ranges)) {
  #   warning("Parsing hotspot overlap introduced missing rows")
  # }
  # gff_ranges$hotspot_overlap_6kb = hotoverlap
  # table(gff_ranges$hotspot_overlap_6kb)
  # 
  # # # 3kb
  # # hotoverlap = countOverlaps(gff_ranges, ldhot_3kb, type = "any")
  # # if (length(hotoverlap) != length(gff_ranges)) {
  # #   warning("Parsing hotspot overlap introduced missing rows")
  # # }
  # # gff_ranges$hotspot_overlap_3kb = hotoverlap
  df_rho = as.data.frame(gff_ranges)
  saveRDS(df_rho, file = paste("Data/Recombination/Gradient/gff_rho_", s, ".rds", sep = ""))
  rm(df_rho)
}

