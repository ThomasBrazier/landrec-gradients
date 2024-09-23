#!/usr/bin/env Rscript
source("Source/init.R")


# The function to smooth maps ======================================================

ldmap_win = function(win.size) {
  ld_allmaps = data.table()
  
  for (ds in list_dataset) {
    maps = chromosome_metadata[which(chromosome_metadata$set == ds),]
    gfffile = paste("Data/Genomic_landscapes/GFF_parsed/", ds, ".csv.gz", sep = "")
    
    if (file.exists(gfffile)) {
      gff = read.table(gzfile(gfffile), header = T, sep = "\t", quote = "")
      gff = gff[which(gff$seqname %in% unique(maps$ldmapname)),]
      gff = gff[which(gff$feature == "gene"),]
      gff = GRanges(gff)
      # gff_startpos = gff
      # end(gff_startpos) = start(gff_startpos) + 1
      cat(ds, "\n")
      for (i in 1:nrow(maps)) {
        chromosome = maps$ldmapname[i]
        cat(chromosome, "\n")
        mapfile = paste("Data/Recombination/LD/ldhat/", ds, ".", chromosome, ".bpen", bpen,".res.txt.gz", sep = "")
        if (file.exists(mapfile)) {
          # Import LD map
          ldmap = read.ldmap(ds, chromosome)
          ldmap = makeGRangesFromDataFrame(ldmap, keep.extra.columns = TRUE)
          
          # Create a grid of 10kb windows
          pos = seq(1, max(end(ldmap), na.rm = TRUE), by = win.size)
          grid = GRanges(seqnames = chromosome,
                         IRanges(start = pos,
                                 width =  win.size),
                         strand = "*")
          
          
          # Estimate weighted mean rho in windows
          wmean = function(x) {
            hits = findOverlaps(grid[x], ldmap)
            me = ldmap$Mean_rho[subjectHits(hits)]
            wi = width(ldmap[subjectHits(hits)])
            wme = weighted.mean(me, wi, na.rm = TRUE)
            return(wme)
          }
          
          grid$weighted.mean.rho.kb = unlist(pbmclapply(X = 1:length(grid), wmean))
          grid = grid[!is.na(grid$weighted.mean.rho.kb)]
          
          # Number of intervals
          nintervals = function(x) {
            hits = findOverlaps(grid[x], ldmap)
            ninter = length(hits)
            return(ninter)
          }
          
          grid$n.intervals = unlist(pbmclapply(X = 1:length(grid), nintervals))
          
          # SNP density
          grid$snp.density.bp = grid$n.intervals/win.size
          
          # Gene density
          # idx = which(seqnames(gff) == chromosome_metadata$annotname[which(chromosome_metadata$set == ds & chromosome_metadata$ldmapname == chromosome)])
          idx = which(seqnames(gff) == chromosome)
          gff_chromosome = as.data.frame(gff[idx])
          # gff_chromosome$seqnames = chromosome
          gff_chromosome = GRanges(gff_chromosome)
          if (sum(is.na(gff_chromosome$gene_biotype)) == 0) {
            gff_chromosome = gff_chromosome[(gff_chromosome$gene_biotype == "protein_coding") | (gff_chromosome$gene_biotype == "")]
          }
          
          genedens = function(x) {
            hits = findOverlaps(grid[x], gff_chromosome)
            gdens = length(hits)/win.size
            return(gdens)
          }
          
          grid$gene.density.bp = unlist(pbmclapply(X = 1:length(grid), genedens))
          
          # Proportion of genic sequence
          propgenic = function(x) {
            hits = findOverlaps(grid[x], gff_chromosome)
            genic = pintersect(grid[x][queryHits(hits)],
                               gff_chromosome[subjectHits(hits)],
                               ignore.strand = TRUE)
            genic = reduce(genic, ignore.strand = TRUE)
            lengenic = sum(width(genic))
            gprop = lengenic/win.size
            return(gprop)
          }
          
          grid$proportion.genic = unlist(pbmclapply(X = 1:length(grid), propgenic))
          
          # Order by descending recombination rate
          # Index and cumulative proportion of genome
          grid = data.table(as.data.frame(grid))
          
          grid = arrange(grid, plyr::desc(weighted.mean.rho.kb))
          
          grid$index = seq(1, nrow(grid))
          grid$cumulative.genome.proportion = grid$index/nrow(grid)
          grid$cumulative.recombination.proportion = cumsum(grid$weighted.mean.rho.kb)/sum(grid$weighted.mean.rho.kb, na.rm = TRUE)
          
          grid$species = paste(unlist(strsplit(ds, "_"))[1:2], collapse = " ")
          grid$chromosome = chromosome
          
          ld_allmaps = rbind(ld_allmaps, grid)
          rm(grid)
          gc()
        }
      }
    }
  }
  
  return(ld_allmaps)
}


cat("======================================\n")
cat("LD maps smoothed at a 10 kb scale\n")

win.size = 10*10^3

ld_maps_10 = ldmap_win(win.size)

saveRDS(ld_maps_10, "Data/Recombination/ldmap_windows_10kb.rds")



cat("======================================\n")
cat("LD maps smoothed at a 100 kb scale\n")


win.size = 100*10^3

ld_maps_100 = ldmap_win(win.size)

saveRDS(ld_maps_100, "Data/Recombination/ldmap_windows_100kb.rds")


