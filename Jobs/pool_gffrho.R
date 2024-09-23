#!/usr/bin/env Rscript
source("Source/init.R")


#============================================================================#
# Debug
#============================================================================#
cat("===================\n")
cat("Recode GFF rho csv file not encoded well\n")

# for (s in list_dataset) {
#   cat(s, "\n")

#   for (c in chromosome_metadata$ldmapname[which(chromosome_metadata$set == s)]) {
#     if (file.exists(paste("Output/gff_rho_", s, "_", c, ".Rda", sep = ""))) {
#       cat(c, "\n")
#       rm(df)
#       rm(gff_ranges)
#       load(file = paste("Output/gff_rho_", s, "_", c, ".Rda", sep = ""))

#       df = as.data.frame(gff_ranges, row.names = NULL)
#       df$weighted.mean.rho.control = unlist(df$weighted.mean.rho.control)

#       cat("Save .csv file\n")
#       write.table(df, file = gzfile(paste("Data/Genomic_landscapes/Rho/gff_rho_", s, "_", c, ".csv.gz", sep = "")),
#                   row.names = F, col.names = T, sep = "\t", quote = F)
# }}}



#============================================================================#
# Pool per species
#============================================================================#
cat("===================\n")
cat("Pool GFF per species\n")

for (s in list_dataset) {
  cat(s, "\n")
  rm(df_rho)
  for (c in chromosome_metadata$ldmapname[which(chromosome_metadata$set == s)]) {
    if (file.exists(paste("Data/Genomic_landscapes/Rho/gff_rho_", s, "_", c, ".csv.gz", sep = ""))) {
      cat(c, "\n")
      # load(paste("Data/Genomic_landscapes/Rho/gff_rho_", s, "_", c, ".Rda", sep = ""))
      gff_ranges = readr::read_tsv(gzfile(paste("Data/Genomic_landscapes/Rho/gff_rho_", s, "_", c, ".csv.gz", sep = "")))
      tmp = as.data.frame(gff_ranges, row.names = NULL)
      # tmp = tmp[,-which(colnames(tmp) == "X0")]
      # Filter p  rotein coding genes only
      # if ("protein_coding" %in% tmp$gene_biotype) {
      #   pcgenes = tmp$id[which(tmp$gene_biotype == "protein_coding")]
      #   tmp = tmp[which(tmp$id %in% pcgenes | tmp$parent %in% pcgenes | tmp$parent %in% tmp$id[which(tmp$parent %in% pcgenes)]),]
      # }

      # Keep only a single transcript per gene
      if (!(s %in% c("Brassica_napus_Wu2019", "Citrullus_lanatus_Guo2019"))) {
        # List of transcripts (mRNA) per gene
        transcripts = data.frame(mRNA = tmp$id[which(tmp$feature == "mRNA")],
                                 gene = tmp$gene_id[which(tmp$feature == "mRNA")])

        # Which transcripts to remove?
        # Keep only first transcript per gene
        transcr = split(transcripts, transcripts$gene)
        transcr_filtered = rbindlist(pbmclapply(transcr, function(x) {x[order(x$mRNA),][-1,]}))
        # Drop duplicated transcripts
        if (nrow(transcr_filtered) > 0) {
          tmp = tmp[-which(tmp$id %in% transcr_filtered$mRNA | tmp$parent %in% transcr_filtered$mRNA),]
        }
      }


      tmp$gene_length_bp = abs(tmp$end - tmp$start + 1)
      gff_gene = tmp[which(tmp$feature %in% c("gene", "mRNA")),]
      # gene_id_trim = unique(gff_gene$gene_id[which(gff_gene$nb_exons > max.exons | gff_gene$snp_count < 2 | gff_gene$gene_length_bp > 10000)])
      gene_id_trim = unique(gff_gene$gene_id[which(gff_gene$nb_exons > max.exons)])

      if (!(s %in% c("Brassica_napus_Wu2019", "Citrullus_lanatus_Guo2019"))) {
        tmp = tmp[-which(tmp$gene_id %in% gene_id_trim),]
      }

      rm(gff_ranges)
      
      colnames(tmp)[which(colnames(tmp) == "seqnames")] = "chromosome"
      tmp$chromosome = as.character(tmp$chromosome)
      
      # Count the number of SNPs
      if (file.exists(paste0("Data/Recombination/LD/ldhat/", s, ".", c, ".bpen5.res.txt.gz"))) {
        gffrange = tmp
        gffrange$chromosome = c
        gffrange = makeGRangesFromDataFrame(gffrange, keep.extra.columns = TRUE)
        # Load the LD map
        ldmap = read.ldmap(s, c)
        ldmap = makeGRangesFromDataFrame(ldmap, keep.extra.columns = TRUE)
        end(ldmap) = start(ldmap) + 1
        snp_count = countOverlaps(gffrange, ldmap, type = "any")
        
        tmp$nb_snp = snp_count
      } else {
        tmp$nb_snp = NA
      }
      

      
      if (exists("df_rho")) {
        df_rho = dplyr::bind_rows(df_rho, tmp)
        rm(tmp)
      } else {
        df_rho = tmp
        rm(tmp)
      }
    }
  }
  saveRDS(df_rho, file = paste("Data/Genomic_landscapes/Rho/gff_rho_", s, ".rds", sep = ""))
  rm(df_rho)
}

gc()


#============================================================================#
# Hotspot overlap ----
#============================================================================#
# Add a new column for Hotspot overlap based on trimmed data
# cat("======================================\n")
# cat("Hotspot overlap")
# 
# for (s in list_dataset) {
#   cat(s, "\n")
#   gff_recombination = readRDS(paste("Data/Recombination/Gradient/gff_rho_", s, ".rds", sep = ""))
#   # Import LD hotspots
#   # Raw data
#   ldhot = read.table(file = gzfile("Data/Recombination/LD/ldhotspots_raw.csv.gz"), header = T)
#   ldhot = ldhot[which(ldhot$set == s),]
#   ldhot = makeGRangesFromDataFrame(ldhot,
#                                    keep.extra.columns=TRUE)
#   # Translate chromosome names
#   chrnames = seqlevels(ldhot)
#   for (j in 1:length(chrnames)) {
#     # cat(chr, '\n')
#     chrnames[j] = chromosome_metadata$annotname[which(chromosome_metadata$ldmapname == chrnames[j] & chromosome_metadata$set == s)]
#   }
#   seqlevels(ldhot) = chrnames
#   
#   # Trimmed
#   ldhot_filtered = read.table(file = gzfile("Data/Recombination/LD/ldhotspots_filtered.csv.gz"), header = T)
#   # ldhot_3kb = read.table(file = gzfile("Data/Recombination/LD/ldhotspots_3kb.csv.gz"), header = T)
#   ldhot_filtered = ldhot_filtered[which(ldhot_filtered$set == s),]
#   # ldhot_3kb = ldhot_3kb[which(ldhot_3kb$set == s),]
#   ldhot_filtered = makeGRangesFromDataFrame(ldhot_filtered,
#                                        keep.extra.columns=TRUE)
#   # ldhot_3kb = makeGRangesFromDataFrame(ldhot_3kb,
#   #                                      keep.extra.columns=TRUE)
#   seqlevels(ldhot_filtered) = chrnames
#   # seqlevels(ldhot_3kb) = chrnames
#   
#   # GFF ranges
#   gff_ranges = makeGRangesFromDataFrame(gff_recombination,
#                                         keep.extra.columns=TRUE)
#   
#   # Assess if feature overlaps hotspot/hotspot position
#   hotoverlap = countOverlaps(gff_ranges, ldhot, type = "any")
#   if (length(hotoverlap) != length(gff_ranges)) {
#     warning("Parsing hotspot overlap introduced missing rows")
#   }
#   gff_ranges$hotspot_overlap = hotoverlap
#   table(gff_ranges$hotspot_overlap)
#   
#   # Trimmed hotspots
#   # 6kb
#   hotoverlap = countOverlaps(gff_ranges, ldhot_filtered, type = "any")
#   if (length(hotoverlap) != length(gff_ranges)) {
#     warning("Parsing hotspot overlap introduced missing rows")
#   }
#   gff_ranges$hotspot_overlap_filtered = hotoverlap
#   table(gff_ranges$hotspot_overlap_filtered)
#   
#   # # 3kb
#   # hotoverlap = countOverlaps(gff_ranges, ldhot_3kb, type = "any")
#   # if (length(hotoverlap) != length(gff_ranges)) {
#   #   warning("Parsing hotspot overlap introduced missing rows")
#   # }
#   # gff_ranges$hotspot_overlap_3kb = hotoverlap
#   gff_ranges = as.data.frame(gff_ranges)
#   saveRDS(gff_ranges, file = paste("Data/Recombination/Gradient/gff_rho_", s, ".rds", sep = ""))
#   rm(gff_ranges)
# }
# gc()





#============================================================================#
# GFF Rho all in one dataset ----
#============================================================================#
cat("======================================\n")
cat("GFF rho all in one")

# max_number_exons = 60

# Pool GFF Rho in a single dataset
rm(df_rho)
for (s in list_dataset) {
  cat(s, "\n")
  df_tmp = readRDS(file = paste("Data/Genomic_landscapes/Rho/gff_rho_", s, ".rds", sep = ""))
  df_tmp$set = s
  df_tmp$species = paste0(strsplit(s, "_")[[1]][1], " ", strsplit(s, "_")[[1]][2])

  df_tmp$weighted.mean.rho.control = as.numeric(df_tmp$weighted.mean.rho.control)

  # Add a simulated gradient ----
  cat("Add a simulated gradient.\n")
  df_tmp$simulated.gradient = NA
  
  # Simulation. Power to detect a true gradient if there is one.
  # Re-ordering mean Rho per ascending order (within genes)
  reorder_ranks = function(id) {
    # A function that takes a gene and reorder Rho values
    # of introns/exons by descending order (simulated gradient)
    # x$simulated.gradient = NA
    # cat(id, "\n")
    x = which(df_tmp$gene_id == id & df_tmp$rank > 0 & df_tmp$rank <= max.exons & df_tmp$feature == "CDS")
    ordered_rho = sort(df_tmp$weighted.mean.rho[x], decreasing = TRUE, na.last = TRUE)
    # order_rank = df_tmp$rank[x]
    # ordered_rho = ordered_rho[order_rank]
    
    # if (length(df_tmp$gene_id == id & df_tmp$rank > 0 & df_tmp$rank <= max.exons & df_tmp$feature == "exon") > 0) {
    #   x = which(df_tmp$gene_id == id & df_tmp$rank > 0 & df_tmp$rank <= max.exons & df_tmp$feature == "exon")
    #   ordered_rho = sort(df_tmp$weighted.mean.rho[x], decreasing = TRUE, na.last = TRUE)
    #   order_rank = df_tmp$rank[x]
    #   ordered_rho = ordered_rho[order_rank]
    #   # df_tmp$simulated.gradient[x] = ordered_rho
    # }
    # if (length(df_tmp$gene_id == id & df_tmp$rank > 0 & df_tmp$rank <= max.exons & df_tmp$feature == "CDS") > 0) {
    #   x = which(df_tmp$gene_id == id & df_tmp$rank > 0 & df_tmp$rank <= max.exons & df_tmp$feature == "CDS")
    #   ordered_rho = sort(df_tmp$weighted.mean.rho[x], decreasing = TRUE, na.last = TRUE)
    #   order_rank = df_tmp$rank[x]
    #   ordered_rho = ordered_rho[order_rank]
    #   # df_tmp$simulated.gradient[x] = ordered_rho
    #   }
    # if (length(df_tmp$gene_id == id & df_tmp$rank > 0 & df_tmp$rank <= max.exons & df_tmp$feature == "intron") > 0) {
    #   x = which(df_tmp$gene_id == id & df_tmp$rank > 0 & df_tmp$rank <= max.exons & df_tmp$feature == "intron")
    #   ordered_rho = sort(df_tmp$weighted.mean.rho[x], decreasing = TRUE, na.last = TRUE)
    #   order_rank = df_tmp$rank[x]
    #   ordered_rho = ordered_rho[order_rank]
      # df_tmp$simulated.gradient[x] = ordered_rho
    # }
    return(ordered_rho)
  }
  
  reorder_ranks_idx = function(id) {
    # A function that takes a gene and reorder Rho values
    # of introns/exons by descending order (simulated gradient)
    # x$simulated.gradient = NA
    # cat(id, "\n")
    x = which(df_tmp$gene_id == id & df_tmp$rank > 0 & df_tmp$rank <= max.exons & df_tmp$feature == "CDS")
    
    # if (length(df_tmp$gene_id == id & df_tmp$rank > 0 & df_tmp$rank <= max.exons & df_tmp$feature == "exon") > 0) {
    #   x = which(df_tmp$gene_id == id & df_tmp$rank > 0 & df_tmp$rank <= max.exons & df_tmp$feature == "exon")
    #       }
    # if (length(df_tmp$gene_id == id & df_tmp$rank > 0 & df_tmp$rank <= max.exons & df_tmp$feature == "CDS") > 0) {
    #   x = which(df_tmp$gene_id == id & df_tmp$rank > 0 & df_tmp$rank <= max.exons & df_tmp$feature == "CDS")
    #   }
    # if (length(df_tmp$gene_id == id & df_tmp$rank > 0 & df_tmp$rank <= max.exons & df_tmp$feature == "intron") > 0) {
    #   x = which(df_tmp$gene_id == id & df_tmp$rank > 0 & df_tmp$rank <= max.exons & df_tmp$feature == "intron")
    #   }
    return(x)
  }

  list_genes = unique(df_tmp$gene_id)
  
  new_rank = unlist(pbmclapply(list_genes, function(x) {reorder_ranks(x)}))
  idx = unlist(pbmclapply(list_genes, function(x) {reorder_ranks_idx(x)}))
  
  df_tmp$simulated.gradient[idx] = new_rank


  
  # Combine datasets
  if (exists("df_rho")) {
    df_rho = dplyr::bind_rows(df_rho, df_tmp)
    # rm(new_order)
    rm(df_tmp)
  } else {
    df_rho = df_tmp
    # rm(new_order)
    rm(df_tmp)
  }
}
df_rho


saveRDS(df_rho, file = "Data/Recombination/Gradient/gff_rho_all.rds")
rm(df_rho)
gc()


