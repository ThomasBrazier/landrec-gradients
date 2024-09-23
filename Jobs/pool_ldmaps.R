#!/usr/bin/env Rscript
source("Source/init.R")

source("Source/read.ldmap.R")

cat("==========================\n")
cat("Pool LD maps per species\n")
for (s in list_dataset) {
  cat(s, "\n")
  df = chromosome_metadata[which(chromosome_metadata$set == s),]
  maps = lapply(1:nrow(df),
                function(x) {read.ldmap(s, df$ldmapname[x])})
  maps = rbindlist(maps)
  
  write.table(maps, file = gzfile(paste("Data/Recombination/LD/ldhat/", s, ".csv.gz", sep = "")), col.names = T,
              row.names = F, sep = "\t", quote = F)
  rm(df)
}


cat("==========================\n")
cat("Pool LD maps\n")
rm(maps)
for (s in list_dataset) {
  cat(s, "\n")
  tmp = read.table(file = gzfile(paste("Data/Recombination/LD/ldhat/", s, ".csv.gz", sep = "")), header = TRUE)
  tmp$chromosome = as.character(tmp$chromosome)

  if (exists("maps")) {
    maps = dplyr::bind_rows(maps, tmp)
    rm(tmp)
  } else {
    maps = tmp
    rm(tmp)
  }
}


maps$length = maps$end - maps$start + 1
write.table(maps, file = gzfile(paste("Data/Recombination/LD/ldhat/LD_maps.csv.gz", sep = "")),
            quote = F, col.names = T, row.names = F, sep = "\t")
maps = maps[, c(1,12)]
write.table(maps, file = gzfile(paste("Data/Recombination/LD/ldhat/LD_widths.csv.gz", sep = "")),
           quote = F, col.names = T, row.names = F, sep = "\t")
rm(maps)
gc()

