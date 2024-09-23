#!/usr/bin/env Rscript
source("Source/init.R")

#============================================================================#
# Pool GFF files in a single file ----
#============================================================================#
cat("======================================\n")
cat("Pool GFF all in one\n")

for (s in list_dataset) {
  cat(s, "\n")
  accession = dataset_metadata$accession[which(dataset_metadata$dataset == s)]
  species = gsub(" ", "_", dataset_metadata$species[which(dataset_metadata$dataset == s)])
  cat(accession, "\n")
  gfffile = paste("Data/Genome/", species, "_", accession, "/", species, "_", accession, ".gff.gz", sep = "")
  df_tmp = read.table(gfffile, header = F, sep = "\t", quote = "")
  colnames(df_tmp) = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  
  # TODO Keep only protein-coding genes
  # Or all if not annotation of gene type
  
  df_tmp$dataset = s
  df_tmp$species = species
  df_tmp$accession = accession
  
  # Combine datasets
  if (exists("df_total")) {
    df_total = dplyr::bind_rows(df_total, df_tmp)
    rm(df_tmp)
  } else {
    df_total = df_tmp
    rm(df_tmp)
  }
}

saveRDS(df_total, file = "Data/Genome/gff_all.rds")
rm(df_total)
gc()


