########################################################################## #
#                    ECOBIO - PhD
#
# 
#         Generate LD maps figures ----
#
########################################################################## #
rm(list=ls(all=TRUE))

# source("Source/init.R")


#==================================================================#
# - [x] LDmap landscape (Rho/kb) ----
# s = "Arabidopsis_thaliana_1001genomes"
# chr = "1"
cat("LDmap landscape (Rho/kb)")

source("Source/init.R")

source("Source/read.ldmap.R")
source("Source/set_ggplot_theme.R")


pdf("Figure/Paper/FigS1.pdf", width = 5, height = 5, pointsize = 9)

for (s in sort(list_dataset)) {
  cat(s, "\n")
  list_chr = chromosome_metadata$ldmapname[which(chromosome_metadata$set == s)]
  for (chr in list_chr) {
    cat(chr, "\n")
    if (file.exists(paste0("Data/Recombination/LD/ldhat/", s, ".", chr, ".bpen", bpen, ".res.txt.gz"))) {
      ldmap = read.ldmap(s, chr, 5)
      
      p1 = ggplot(data = ldmap, aes(x = end/10^6, y = Mean_rho_unfiltered)) +
        geom_line(alpha = 0.4) +
        # geom_step(aes(x = start/10^6, y = Mean_rho), alpha = 0.3) +
        geom_segment(data = ldmap[which(is.na(ldmap$Mean_rho)),], aes(x = start/10^6, xend = end/10^6, y = Mean_rho_unfiltered, yend = Mean_rho_unfiltered), color = "Red") +
        xlab("Genomic position (Mb)") + ylab("Mean Rho/kb") +
        ggtitle(paste(s, "chromosome", chr))
      
      # Which Marey map for the LD map?
      mset = marey_data$set[which(marey_data$species == gsub(" ", "_",dataset_metadata$species[which(dataset_metadata$dataset == s)]))]
      mchr = chromosome_metadata$litteralname[which(chromosome_metadata$set == mset & chromosome_metadata$ldmapname == chr)]
      
      mapfile = paste0("Data/Recombination/Marey/", mset, "_chromosome", mchr, ".txt")
      # mapfile = "Data/Recombination/Marey/Arabidopsis_thaliana_Serin2017_chromosome1.txt"
        
      if (file.exists(mapfile)) {
        
        mareymap = read.table(mapfile, header = T, sep = "\t")
        
        p2 = ggplot(mareymap, aes(x = phys, y = rec.rate)) +
          geom_line() +
          geom_ribbon(aes(x = phys, ymin = lower, ymax= upper), alpha = 0.2) +
          xlab("Genomic position (Mb)") + ylab("cM/Mb")
        
        rm(mareymap)
      } else {
        p2 = ggplot()
      }
      
      p = ggarrange(p1, p2, nrow = 2, labels = "AUTO")
      
      print(p)
      
      rm(p1)
      rm(p2)
      rm(p)
      rm(ldmap)
      gc()
    }
  }
}

dev.off()
