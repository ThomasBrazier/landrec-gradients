########################################################################## #
#     JOB - Figures for LD recombination maps with LDhat
########################################################################## #

#============================================================================#
# LOADING ENVIRONMENT ----
#============================================================================#
# clear global environment: remove all variables
rm(list=ls(all=TRUE))
library(readODS)
library(ggplot2)
library(hrbrthemes)
library(viridis)
library(ggpubr)
library(tidyr)
library(tidyverse)
library(gridExtra)
library(purrr)

#============================================================================#
# Loading variables & objects ----
#============================================================================#
# Get the directory of the file & set working directory
wd=dirname(rstudioapi::getSourceEditorContext()$path)
wd=gsub("/Source", "", wd)
setwd(wd)

#============================================================================#
# Get the list of species and chromosomes to process ----
#============================================================================#
meetadata_sets = read_ods("Data/LDRecMaps_TODO.ods", sheet = 1, na = "")
# complete_sets = complete_sets[which(complete_sets$done == "yes"),c("set", "chromosome")]

#============================================================================#
# Import gff files ----
#============================================================================#
gff_set = list("Oryza_sativa_DeLeon2016")
# gff_set = unique(complete_sets$set)
# Read each csv file into a list
read_func <- function(x){
  file = paste("Data/Genomic_landscapes/GC_genes/", x, ".csv.gz", sep = "")
  gff = read.table(gzfile(file), header = T, sep = "\t")
  return(gff)
}
df = purrr::map(gff_set, read_func)
df = bind_rows(df, .id = 'set')



#============================================================================#
# Get the list of data produced and make a metadata file ----
#============================================================================#
# - set
# - chromosome
# - bpen
list_res = list.files(paste(wd, "/Data/Recombination/LD/ldhat", sep = ""))
res_bpen5 = list_res[which(grepl("bpen5", list_res))]
res_bpen15 = list_res[which(grepl("bpen15", list_res))]
res_bpen25 = list_res[which(grepl("bpen25", list_res))]

complete_sets = data.frame(file = c(res_bpen5, res_bpen15, res_bpen25),
                           set = NA,
                           chromosome = NA,
                           bpen = c(rep(5, length(res_bpen5)),
                                    rep(15, length(res_bpen15)),
                                    rep(25, length(res_bpen25))))
complete_sets$set = gsub(".[A-Za-z0-9]*.bpen[0-9]*.res.txt.gz", "", complete_sets$file)
complete_sets$chromosome =  gsub("^[A-Za-z0-9_]*.", "", gsub(".bpen[0-9]*.res.txt.gz", "", complete_sets$file))
  
#============================================================================#
# Print LD recombination maps - LDHAT ----
#============================================================================#
# Save all Marey maps and recombination maps in supplementary
# Create multiple plots using lapply()
p = lapply(1:nrow(complete_sets), function(x) {
  # Import LD recombination map
  file = paste("Data/Recombination/LD/ldhat/", complete_sets$set[x], ".",
               complete_sets$chromosome[x], ".bpen", 
               complete_sets$bpen[x], ".res.txt.gz", sep = "")
  ldmap = read.table(gzfile(file), header = T)
  # Put positions in bp
  ldmap$Loci = ldmap$Loci * 1000
  # Remove negative positions
  ldmap = ldmap[-which(ldmap$Loci < 0),]
  ldmap$start = c(0, head(ldmap$Loci, -1))
  ldmap$end = (ldmap$Loci - 1)
  # Make the figure as a ggplot object
  ggplot(data = ldmap, aes(x = Loci/10^6, y = Mean_rho)) +
    ggtitle(paste(complete_sets$set[x],
                  "chromosome", complete_sets$chromosome[x],
                  "bpen", complete_sets$bpen[x], sep = " ")) +
    geom_line() + # Add markers
    xlab("Genomic position (Mb)") + ylab("Rho") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(color="black", size=14,hjust = 0.5),
          axis.title.x = element_text(color="black", size=14),
          axis.title.y = element_text(color="black", size=14),
          axis.text=element_text(size=14, colour="black"),
          legend.key = element_rect(fill = "white", size = 1),
          legend.key.height = unit(2,"line"),
          legend.key.width = unit(5,"line"),
          legend.text=element_text(size=14, face = "italic"),
          legend.title=element_text(size=14),
          legend.position='right')
})

# Save list of plots
ggsave(filename = paste(wd, "/Figure/LD_maps.pdf", sep = ""), 
       plot = marrangeGrob(p, nrow=1, ncol=1), 
       width = 11, height = 8
)



#============================================================================#
# Print LD recombination maps - LDHAT ----
#============================================================================#
# 1 - Marey map on top
# 2 - LD maps for bpen = 5,15,25
# 3 - SNP density on bottom





#============================================================================#
# Print LD recombination maps - LDHAT 100kb averaged ----
#============================================================================#
p = lapply(1:nrow(complete_sets), function(x) {
  # Import LD recombination map
  file = paste("Data/Recombination/LD/ldhat/", complete_sets$set[x], ".",
               complete_sets$chromosome[x], ".bpen", 
               complete_sets$bpen[x], ".res.txt.gz", sep = "")
  ldmap = read.table(gzfile(file), header = T)
  # Put positions in bp
  ldmap$Loci = ldmap$Loci * 1000
  # Remove negative positions
  ldmap = ldmap[-which(ldmap$Loci < 0),]
  ldmap$start = c(0, head(ldmap$Loci, -1))
  ldmap$end = (ldmap$Loci - 1)
  # Recombination rates averaged in 100kb windows
  broadmap = data.frame(start = seq(1, max(ldmap$Loci), by = 100000))
  broadmap$end = c(broadmap$start[-1] - 1, max(ldmap$Loci))
  recaveraged = function(start, end) {
    mean(ldmap$Mean_rho[which(ldmap$Loci > start & ldmap$Loci < end)])
  }
  broadmap$recaveraged = mapply(recaveraged, broadmap$start, broadmap$end)
  broadmap$pos = ((broadmap$end + broadmap$start)/2)/10^6

  # Make the figure as a ggplot object
  ggplot(data = broadmap, aes(x = pos, y = recaveraged)) +
    ggtitle(paste(complete_sets$set[x],
                  "chromosome", complete_sets$chromosome[x],
                  "bpen", complete_sets$bpen[x], sep = " ")) +
    geom_line() + # Add markers
    xlab("Genomic position (Mb)") + ylab("Average Rho (100kb)") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(color="black", size=14,hjust = 0.5),
          axis.title.x = element_text(color="black", size=14),
          axis.title.y = element_text(color="black", size=14),
          axis.text=element_text(size=14, colour="black"),
          legend.key = element_rect(fill = "white", size = 1),
          legend.key.height = unit(2,"line"),
          legend.key.width = unit(5,"line"),
          legend.text=element_text(size=14, face = "italic"),
          legend.title=element_text(size=14),
          legend.position='right')
})

# Save list of plots
ggsave(filename = paste(wd, "/Figure/LD_maps_100kb.pdf", sep = ""), 
       plot = marrangeGrob(p, nrow=1, ncol=1), 
       width = 11, height = 8
)

#============================================================================#
# Print recombination gradients and GC gradients in genes - LDHAT ----
#============================================================================#
# Import GFF
# file = paste("Data/Genomic_landscapes/GC_genes/", complete_sets$set[x], ".csv.gz", sep = "")
# gff = read.table(gzfile(file), header = T, sep = "\t")





#============================================================================#
# End ----
#============================================================================#