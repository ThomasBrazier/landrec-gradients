#!/usr/bin/env Rscript

source("Source/init.R")

interval = 200
upstream = 5000 # Size of the upstream region
downstream = 5000 # Max size of the downstream region

#============================================================================#
# Pool Rho gradient 5kb ATG ----
#============================================================================#rm(df_pos)
rm(df_set)
# s = "Arabidopsis_thaliana_1001genomes"
for (s in list_dataset) {
  cat(s, "\n")
  rm(df_set)
  for (c in chromosome_metadata$ldmapname[which(chromosome_metadata$set == s)]) {
    cat(c, "\n")
    if (file.exists(paste("Output/RhoGradient_ATG_", s, "_", c, ".Rda", sep = ""))) {
      load(paste("Output/RhoGradient_ATG_", s, "_", c, ".Rda", sep = ""))
      
      if (exists("df_set")) {
        df_set = rbind.fill(df_pos, df_set)
        rm(df_pos)
      } else {
        df_set = df_pos
        rm(df_pos)
      }
    }
  }
  saveRDS(df_set, file = paste("Data/Recombination/Gradient/RhoGradient_ATG_", s, ".rds", sep = ""))
  rm(df_set)
}
# df = readRDS(file = paste("Data/Recombination/Gradient/RhoGradient_ATG_", s, ".rds", sep = ""))



#============================================================================#
# Pool Rho gradient 5kb TSS + TTS ----
#============================================================================#
cat("======================================\n")
cat("Rho gradient TSS + TTS")
# Pool Rho gradient per species
for (s in list_dataset) {
  cat(s, "\n")
  rm(df_all)
  for (c in chromosome_metadata$ldmapname[which(chromosome_metadata$set == s)]) {
    if (file.exists(paste("Output/RhoGradient_5kbTTS_", s, "_", c, ".Rda", sep = ""))) {
      cat(c, "\n")
      load(paste("Output/RhoGradient_5kbTTS_", s, "_", c, ".Rda", sep = ""))
      
      if (exists("df_all")) {
        df_all = rbind.fill(df_all, df_pos)
        rm(df_pos)
      } else {
        df_all = df_pos
        rm(df_pos)
      }
    }
  }
  saveRDS(df_all, file = paste("Data/Recombination/Gradient/RhoGradient_5kbTTS_", s, ".rds", sep = ""))
  rm(df_all)
}

for (s in list_dataset) {
  cat(s, "\n")
  rm(df_all)
  for (c in chromosome_metadata$ldmapname[which(chromosome_metadata$set == s)]) {
    if (file.exists(paste("Output/RhoGradient_5kbTSS_", s, "_", c, ".Rda", sep = ""))) {
      cat(c, "\n")
      load(paste("Output/RhoGradient_5kbTSS_", s, "_", c, ".Rda", sep = ""))
      
      if (exists("df_all")) {
        df_all = rbind.fill(df_all, df_pos)
        rm(df_pos)
      } else {
        df_all = df_pos
        rm(df_pos)
      }
    }
  }
  saveRDS(df_all, file = paste("Data/Recombination/Gradient/RhoGradient_5kbTSS_", s, ".rds", sep = ""))
  rm(df_all)
}


#============================================================================#
# Compute the average gradient ----
#============================================================================#
cat("======================================\n")
cat("Compute the average gradient ")

rm(df_gradient)
for (s in list_dataset) {
  cat(s, "\n")
  rm(df_pos)
  rm(df_pos_exons)
  rm(df_pos_introns)
  cat("ATG position\n")
  df_pos = readRDS(paste("Data/Recombination/Gradient/RhoGradient_ATG_", s, ".rds", sep = ""))
  
  # Keep only genes with less than 15 exons
  list_genes = data_all$gene_id[which(data_all$nb_exons <= max.exons)]
  df_pos = df_pos[which(df_pos$gene_id %in% list_genes),]
  
  # Remove non genic parts downstream
  # df_pos = df_pos[-which(df_pos$nb_gene == 0 & df_pos$idx >= 0),]
  # # In fact, keep only 70% of the gene; do not count rho at the end or after gene's end
  # p.keep = 0.7
  # df_pos$relative_pos = df_pos$idx/df_pos$gene_size
  # df_pos = df_pos[which(df_pos$relative_pos <= p.keep),]
  # ATG ----
  
  # Mean Rho over all intervals
  meanRho.ATG = aggregate(wmean.rho ~ idx, data = df_pos, mean)
  meanRho.ATG$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos, length)$wmean.rho
  meanRho.ATG$condition = "observed"
  meanRho.ATG$hotoverlap = "both"
  m.obs = mean(meanRho.ATG$wmean.rho, na.rm = TRUE)
  sd.obs = sd(meanRho.ATG$wmean.rho, na.rm = TRUE)
  
  # meanRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, mean)
  meanRho.control = aggregate(mean.rho.control ~ idx, data = df_pos, mean)
  meanRho.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos, length)$mean.rho.control
  meanRho.control$condition = "control"
  meanRho.control$hotoverlap = "both"
  m.control = mean(meanRho.control$mean.rho.control, na.rm = TRUE)
  sd.control = sd(meanRho.control$mean.rho.control, na.rm = TRUE)
  
  colnames(meanRho.ATG) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  colnames(meanRho.control) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  
  meanRho.ATG = rbind(meanRho.ATG, meanRho.control)
  meanRho.ATG$m.obs = m.obs
  meanRho.ATG$sd.obs = sd.obs
  meanRho.ATG$m.control = m.control
  meanRho.ATG$sd.control = sd.control
  
  meanRho.ATG$meanRho.rescaled = meanRho.ATG$meanRho
    
  meanRho.ATG$meanRho.rescaled[which(meanRho.ATG$condition == "control")] = (meanRho.control$meanRho - m.control) + m.obs
  
  # HOT vs COLD ----
  cat("Discriminate between hot and cold genes\n")
  
  df_pos$hotoverlap = ifelse(df_pos$hotspot_overlap.filtered.4 > 0, "hot", "cold")
  
  # Mean Rho over all intervals
  meanRho.ATG.hot = aggregate(wmean.rho ~ idx + hotoverlap, data = df_pos, mean)
  meanRho.ATG.hot$n.intervals = aggregate(wmean.rho ~ idx + hotoverlap, data = df_pos, length)$wmean.rho
  meanRho.ATG.hot$condition = "observed"
  m.obs = mean(meanRho.ATG.hot$wmean.rho, na.rm = TRUE)
  sd.obs = sd(meanRho.ATG.hot$wmean.rho, na.rm = TRUE)
  
  # meanRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, mean)
  meanRho.control.hot = aggregate(mean.rho.control ~ idx + hotoverlap, data = df_pos, mean)
  meanRho.control.hot$n.intervals = aggregate(mean.rho.control ~ idx + hotoverlap, data = df_pos, length)$mean.rho.control
  meanRho.control.hot$condition = "control"
  m.control = mean(meanRho.control.hot$mean.rho.control, na.rm = TRUE)
  sd.control = sd(meanRho.control.hot$mean.rho.control, na.rm = TRUE)
  
  colnames(meanRho.ATG.hot) = c("idx", "hotoverlap",  "meanRho", "n.intervals", "condition")
  colnames(meanRho.control.hot) = c("idx", "hotoverlap", "meanRho", "n.intervals", "condition")
  
  meanRho.ATG.hot = rbind(meanRho.ATG.hot, meanRho.control.hot)
  meanRho.ATG.hot$m.obs = m.obs
  meanRho.ATG.hot$sd.obs = sd.obs
  meanRho.ATG.hot$m.control = m.control
  meanRho.ATG.hot$sd.control = sd.control
  
  meanRho.ATG.hot$meanRho.rescaled = meanRho.ATG.hot$meanRho
  
  meanRho.ATG.hot$meanRho.rescaled[which(meanRho.ATG.hot$condition == "control")] = (meanRho.control.hot$meanRho - m.control) + m.obs
  
  meanRho.ATG = rbind(meanRho.ATG, meanRho.ATG.hot)
  
  meanRho.ATG$level = "Genic"
  meanRho.ATG$position = "ATG"
  
  # Add exons ----
  cat("Add exons\n")
  df_pos_exons = df_pos[which(df_pos$exon_overlap == TRUE & df_pos$idx > 0),]
  # df_pos_introns = df_pos[which((df_pos$nb_gene == 1) & df_pos$intron_overlap == TRUE),]
  meanRho_exons = aggregate(wmean.rho ~ idx, data = df_pos_exons, mean)
  meanRho_exons$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_exons, length)$wmean.rho
  meanRho_exons$condition = "observed"
  meanRho_exons$hotoverlap = "both"
  m.obs = mean(meanRho_exons$wmean.rho, na.rm = TRUE)
  sd.obs = sd(meanRho_exons$wmean.rho, na.rm = TRUE)
  
  # meanRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, mean)
  meanRho.control_exons = aggregate(mean.rho.control ~ idx, data = df_pos_exons, mean)
  meanRho.control_exons$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_exons, length)$mean.rho.control
  meanRho.control_exons$condition = "control"
  meanRho.control_exons$hotoverlap = "both"
  m.control = mean(meanRho.control_exons$mean.rho.control, na.rm = TRUE)
  sd.control = sd(meanRho.control_exons$mean.rho.control, na.rm = TRUE)
  
  colnames(meanRho_exons) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  colnames(meanRho.control_exons) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  
  meanRho_exons = rbind(meanRho_exons, meanRho.control_exons)
  meanRho_exons$m.obs = m.obs
  meanRho_exons$sd.obs = sd.obs
  meanRho_exons$m.control = m.control
  meanRho_exons$sd.control = sd.control
  
  meanRho_exons$meanRho.rescaled = meanRho_exons$meanRho
  
  meanRho_exons$meanRho.rescaled[which(meanRho_exons$condition == "control")] = (meanRho.control_exons$meanRho - m.control) + m.obs
  
  meanRho_exons$level = "Exons"
  meanRho_exons$position = "ATG"
  
  meanRho.ATG = rbind(meanRho.ATG, meanRho_exons)
  
  # Add introns ----
  cat("Add introns\n")
  df_pos_introns = df_pos[which(df_pos$intron_overlap == TRUE & df_pos$idx > 0),]

  meanRho_introns = aggregate(wmean.rho ~ idx, data = df_pos_introns, mean)
  meanRho_introns$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_introns, length)$wmean.rho
  meanRho_introns$condition = "observed"
  meanRho_introns$hotoverlap = "both"
  m.obs = mean(meanRho_introns$wmean.rho, na.rm = TRUE)
  sd.obs = sd(meanRho_introns$wmean.rho, na.rm = TRUE)
  
  # meanRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, mean)
  meanRho.control_introns = aggregate(mean.rho.control ~ idx, data = df_pos_introns, mean)
  meanRho.control_introns$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_introns, length)$mean.rho.control
  meanRho.control_introns$condition = "control"
  meanRho.control_introns$hotoverlap = "both"
  m.control = mean(meanRho.control_introns$mean.rho.control, na.rm = TRUE)
  sd.control = sd(meanRho.control_introns$mean.rho.control, na.rm = TRUE)
  
  colnames(meanRho_introns) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  colnames(meanRho.control_introns) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  
  meanRho_introns = rbind(meanRho_introns, meanRho.control_introns)
  meanRho_introns$m.obs = m.obs
  meanRho_introns$sd.obs = sd.obs
  meanRho_introns$m.control = m.control
  meanRho_introns$sd.control = sd.control
  
  meanRho_introns$meanRho.rescaled = meanRho_introns$meanRho
  
  meanRho_introns$meanRho.rescaled[which(meanRho_introns$condition == "control")] = (meanRho.control_introns$meanRho - m.control) + m.obs
  
  meanRho_introns$level = "Introns"
  meanRho_introns$position = "ATG"
  
  meanRho.ATG = rbind(meanRho.ATG, meanRho_introns)
  
  # df_pos_exons = df_pos[which((df_pos$nb_gene == 1) & df_pos$exon_overlap == TRUE),]
  # df_pos_introns = df_pos[which((df_pos$nb_gene == 1) & df_pos$intron_overlap == TRUE),]
  # 
  # # mean as a function of position
  # meanRho_exons = aggregate(wmean.rho ~ idx, data = df_pos_exons, mean)
  # meanRho_exons$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_exons, length)$wmean.rho
  # meanRho_exons$condition = "observed"
  # meanRho_exons$hotoverlap = "both"
  # # meanRho_exons.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos_exons, mean)
  # meanRho_exons.control = aggregate(mean.rho.control ~ idx, data = df_pos_exons, mean)
  # meanRho_exons.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_exons, length)$mean.rho.control
  # meanRho_exons.control$condition = "control"
  # meanRho_exons.control$hotoverlap = "both"
  # colnames(meanRho_exons) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  # colnames(meanRho_exons.control) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  # meanRho_exons = rbind(meanRho_exons, meanRho_exons.control)
  # 
  # 
  # meanRho_introns = aggregate(wmean.rho ~ idx, data = df_pos_introns, mean)
  # meanRho_introns$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_introns, length)$wmean.rho
  # meanRho_introns$condition = "observed"
  # meanRho_introns$hotoverlap = "both"
  # # meanRho_introns.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos_introns, mean)
  # meanRho_introns.control = aggregate(mean.rho.control ~ idx, data = df_pos_introns, mean)
  # meanRho_introns.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_introns, length)$mean.rho.control
  # meanRho_introns.control$condition = "control"
  # meanRho_introns.control$hotoverlap = "both"
  # colnames(meanRho_introns) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  # colnames(meanRho_introns.control) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  # meanRho_introns = rbind(meanRho_introns, meanRho_introns.control)
  # 
  # 
  # meanRho.ATG$level = "all"
  # meanRho_exons$level = "exon" 
  # meanRho_introns$level = "intron"
  # 
  # meanRho_exons = meanRho_exons[which(meanRho_exons$idx >= 0),]
  # meanRho_introns = meanRho_introns[which(meanRho_introns$idx >= 0),]
  # 
  # meanRho.ATG = rbind(meanRho.ATG, meanRho_exons, meanRho_introns)
  # 
  # meanRho.ATG$position = "ATG"
  # 
  # # Mean ATG-TSS distance of the species
  atg_tss_dist = mean(df_pos$dist_tss, na.rm = TRUE)
  meanRho.ATG$atg_tss_dist = atg_tss_dist
  # 
  # 
  # # ------------------------------ #
  # TSS ----
  cat("TSS position\n")
  df_pos = readRDS(paste("Data/Recombination/Gradient/RhoGradient_5kbTSS_", s, ".rds", sep = ""))
  
  # Keep only genes with less than 15 exons
  list_genes = data_all$gene_id[which(data_all$nb_exons <= max.exons)]
  df_pos = df_pos[which(df_pos$gene_id %in% list_genes),]
  
  
  meanRho.TSS = aggregate(wmean.rho ~ idx, data = df_pos, mean)
  meanRho.TSS$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos, length)$wmean.rho
  meanRho.TSS$condition = "observed"
  meanRho.TSS$hotoverlap = "both"
  m.obs = mean(meanRho.TSS$wmean.rho, na.rm = TRUE)
  sd.obs = sd(meanRho.TSS$wmean.rho, na.rm = TRUE)
  
  # meanRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, mean)
  meanRho.control = aggregate(mean.rho.control ~ idx, data = df_pos, mean)
  meanRho.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos, length)$mean.rho.control
  meanRho.control$condition = "control"
  meanRho.control$hotoverlap = "both"
  m.control = mean(meanRho.control$mean.rho.control, na.rm = TRUE)
  sd.control = sd(meanRho.control$mean.rho.control, na.rm = TRUE)
  
  colnames(meanRho.TSS) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  colnames(meanRho.control) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  
  meanRho.TSS = rbind(meanRho.TSS, meanRho.control)
  meanRho.TSS$m.obs = m.obs
  meanRho.TSS$sd.obs = sd.obs
  meanRho.TSS$m.control = m.control
  meanRho.TSS$sd.control = sd.control
  
  meanRho.TSS$meanRho.rescaled = meanRho.TSS$meanRho
  
  meanRho.TSS$meanRho.rescaled[which(meanRho.TSS$condition == "control")] = (meanRho.control$meanRho - m.control) + m.obs
  
  meanRho.TSS$position = "TSS"
  meanRho.TSS$level = "Genic"
  
  # df_pos$hotoverlap = ifelse(df_pos$hotspot_overlap.filtered.4 > 0, "hot", "cold")
  # 
  # meanRho.TSS = aggregate(wmean.rho ~ idx, data = df_pos, mean)
  # meanRho.TSS$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos, length)$wmean.rho
  # meanRho.TSS$condition = "observed"
  # meanRho.TSS$hotoverlap = "both"
  # # meanRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, mean)
  # meanRho.control = aggregate(mean.rho.control ~ idx, data = df_pos, mean)
  # meanRho.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos, length)$mean.rho.control
  # meanRho.control$condition = "control"
  # meanRho.control$hotoverlap = "both"
  # 
  # colnames(meanRho.TSS) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  # colnames(meanRho.control) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  # meanRho.TSS = rbind(meanRho.TSS, meanRho.control)
  # 
  # 
  # meanRho.TSS.hot = aggregate(wmean.rho ~ idx + hotoverlap, data = df_pos, mean)
  # meanRho.TSS.hot$n.intervals = aggregate(wmean.rho ~ idx + hotoverlap, data = df_pos, length)$wmean.rho
  # meanRho.TSS.hot$condition = "observed"
  # # meanRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, mean)
  # meanRho.control.hot = aggregate(mean.rho.control ~ idx + hotoverlap, data = df_pos, mean)
  # meanRho.control.hot$n.intervals = aggregate(mean.rho.control ~ idx + hotoverlap, data = df_pos, length)$mean.rho.control
  # meanRho.control.hot$condition = "control"
  # 
  # colnames(meanRho.TSS.hot) = c("idx", "hotoverlap",  "meanRho", "n.intervals", "condition")
  # colnames(meanRho.control.hot) = c("idx", "hotoverlap", "meanRho", "n.intervals", "condition")
  # meanRho.TSS = rbind(meanRho.TSS, meanRho.TSS.hot, meanRho.control.hot)
  # 
  # # Remove intervals overlapping two gene or more than the gene size limit
  # df_pos_exons = df_pos[which((df_pos$nb_gene == 1) & df_pos$exon_overlap == TRUE),]
  # df_pos_introns = df_pos[which((df_pos$nb_gene == 1) & df_pos$intron_overlap == TRUE),]
  # 
  # # mean as a function of position
  # meanRho_exons = aggregate(wmean.rho ~ idx, data = df_pos_exons, mean)
  # meanRho_exons$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_exons, length)$wmean.rho
  # meanRho_exons$condition = "observed"
  # meanRho_exons$hotoverlap = "both"
  # # meanRho_exons.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos_exons, mean)
  # meanRho_exons.control = aggregate(mean.rho.control ~ idx, data = df_pos_exons, mean)
  # meanRho_exons.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_exons, length)$mean.rho.control
  # meanRho_exons.control$condition = "control"
  # meanRho_exons.control$hotoverlap = "both"
  # colnames(meanRho_exons) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  # colnames(meanRho_exons.control) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  # meanRho_exons = rbind(meanRho_exons, meanRho_exons.control)
  # 
  # 
  # meanRho_introns = aggregate(wmean.rho ~ idx, data = df_pos_introns, mean)
  # meanRho_introns$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_introns, length)$wmean.rho
  # meanRho_introns$condition = "observed"
  # meanRho_introns$hotoverlap = "both"
  # # meanRho_introns.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos_introns, mean)
  # meanRho_introns.control = aggregate(mean.rho.control ~ idx, data = df_pos_introns, mean)
  # meanRho_introns.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_introns, length)$mean.rho.control
  # meanRho_introns.control$condition = "control"
  # meanRho_introns.control$hotoverlap = "both"
  # colnames(meanRho_introns) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  # colnames(meanRho_introns.control) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  # meanRho_introns = rbind(meanRho_introns, meanRho_introns.control)
  # 
  # 
  # meanRho.TSS$level = "all"
  # meanRho_exons$level = "exon" 
  # meanRho_introns$level = "intron"
  # 
  # meanRho_exons = meanRho_exons[which(meanRho_exons$idx >= 0),]
  # meanRho_introns = meanRho_introns[which(meanRho_introns$idx >= 0),]
  # 
  # meanRho.TSS = rbind(meanRho.TSS, meanRho_exons, meanRho_introns)
  # 

  # 
  # # ------------------------------ #
  # TTS ----
  cat("TTS position\n")
  df_pos = readRDS(paste("Data/Recombination/Gradient/RhoGradient_5kbTTS_", s, ".rds", sep = ""))
  
  # Keep only genes with less than 15 exons
  list_genes = data_all$gene_id[which(data_all$nb_exons <= max.exons)]
  df_pos = df_pos[which(df_pos$gene_id %in% list_genes),]
  
  meanRho.TTS = aggregate(wmean.rho ~ idx, data = df_pos, mean)
  meanRho.TTS$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos, length)$wmean.rho
  meanRho.TTS$condition = "observed"
  meanRho.TTS$hotoverlap = "both"
  m.obs = mean(meanRho.TTS$wmean.rho, na.rm = TRUE)
  sd.obs = sd(meanRho.TTS$wmean.rho, na.rm = TRUE)
  
  # meanRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, mean)
  meanRho.control = aggregate(mean.rho.control ~ idx, data = df_pos, mean)
  meanRho.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos, length)$mean.rho.control
  meanRho.control$condition = "control"
  meanRho.control$hotoverlap = "both"
  m.control = mean(meanRho.control$mean.rho.control, na.rm = TRUE)
  sd.control = sd(meanRho.control$mean.rho.control, na.rm = TRUE)
  
  colnames(meanRho.TTS) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  colnames(meanRho.control) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  
  meanRho.TTS = rbind(meanRho.TTS, meanRho.control)
  meanRho.TTS$m.obs = m.obs
  meanRho.TTS$sd.obs = sd.obs
  meanRho.TTS$m.control = m.control
  meanRho.TTS$sd.control = sd.control
  
  meanRho.TTS$meanRho.rescaled = meanRho.TTS$meanRho
  
  meanRho.TTS$meanRho.rescaled[which(meanRho.TTS$condition == "control")] = (meanRho.control$meanRho - m.control) + m.obs
  
  meanRho.TTS$position = "TTS"
  meanRho.TTS$level = "Genic"
  
  # df_pos$hotoverlap = ifelse(df_pos$hotspot_overlap.filtered.4 > 0, "hot", "cold")
  # 
  # meanRho.TTS = aggregate(wmean.rho ~ idx, data = df_pos, mean)
  # meanRho.TTS$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos, length)$wmean.rho
  # meanRho.TTS$condition = "observed"
  # meanRho.TTS$hotoverlap = "both"
  # # meanRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, mean)
  # meanRho.control = aggregate(mean.rho.control ~ idx, data = df_pos, mean)
  # meanRho.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos, length)$mean.rho.control
  # meanRho.control$condition = "control"
  # meanRho.control$hotoverlap = "both"
  # 
  # colnames(meanRho.TTS) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  # colnames(meanRho.control) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  # meanRho.TTS = rbind(meanRho.TTS, meanRho.control)
  # 
  # 
  # meanRho.TTS.hot = aggregate(wmean.rho ~ idx + hotoverlap, data = df_pos, mean)
  # meanRho.TTS.hot$n.intervals = aggregate(wmean.rho ~ idx + hotoverlap, data = df_pos, length)$wmean.rho
  # meanRho.TTS.hot$condition = "observed"
  # # meanRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, mean)
  # meanRho.control.hot = aggregate(mean.rho.control ~ idx + hotoverlap, data = df_pos, mean)
  # meanRho.control.hot$n.intervals = aggregate(mean.rho.control ~ idx + hotoverlap, data = df_pos, length)$mean.rho.control
  # meanRho.control.hot$condition = "control"
  # 
  # colnames(meanRho.TTS.hot) = c("idx", "hotoverlap",  "meanRho", "n.intervals", "condition")
  # colnames(meanRho.control.hot) = c("idx", "hotoverlap", "meanRho", "n.intervals", "condition")
  # meanRho.TTS = rbind(meanRho.TTS, meanRho.TTS.hot, meanRho.control.hot)
  # 
  # # Remove intervals overlapping two gene or more than the gene size limit
  # # df_pos_exons = df_pos[which((df_pos$nb_gene == 1) & df_pos$exon_overlap == TRUE),]
  # # df_pos_introns = df_pos[which((df_pos$nb_gene == 1) & df_pos$intron_overlap == TRUE),]
  # df_pos_exons = df_pos[which(df_pos$exon_overlap == TRUE),]
  # df_pos_introns = df_pos[which(df_pos$intron_overlap == TRUE),]
  # 
  # # mean as a function of position
  # meanRho_exons = aggregate(wmean.rho ~ idx, data = df_pos_exons, mean)
  # meanRho_exons$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_exons, length)$wmean.rho
  # meanRho_exons$condition = "observed"
  # meanRho_exons$hotoverlap = "both"
  # # meanRho_exons.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos_exons, mean)
  # meanRho_exons.control = aggregate(mean.rho.control ~ idx, data = df_pos_exons, mean)
  # meanRho_exons.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_exons, length)$mean.rho.control
  # meanRho_exons.control$condition = "control"
  # meanRho_exons.control$hotoverlap = "both"
  # colnames(meanRho_exons) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  # colnames(meanRho_exons.control) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  # meanRho_exons = rbind(meanRho_exons, meanRho_exons.control)
  # 
  # 
  # meanRho_introns = aggregate(wmean.rho ~ idx, data = df_pos_introns, mean)
  # meanRho_introns$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_introns, length)$wmean.rho
  # meanRho_introns$condition = "observed"
  # meanRho_introns$hotoverlap = "both"
  # # meanRho_introns.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos_introns, mean)
  # meanRho_introns.control = aggregate(mean.rho.control ~ idx, data = df_pos_introns, mean)
  # meanRho_introns.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_introns, length)$mean.rho.control
  # meanRho_introns.control$condition = "control"
  # meanRho_introns.control$hotoverlap = "both"
  # colnames(meanRho_introns) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  # colnames(meanRho_introns.control) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
  # meanRho_introns = rbind(meanRho_introns, meanRho_introns.control)
  # 
  # 
  # meanRho.TTS$level = "all"
  # meanRho_exons$level = "exon" 
  # meanRho_introns$level = "intron"
  # 
  # meanRho_exons = meanRho_exons[which(meanRho_exons$idx >= 0),]
  # meanRho_introns = meanRho_introns[which(meanRho_introns$idx >= 0),]
  # 
  # meanRho.TTS = rbind(meanRho.TTS, meanRho_exons, meanRho_introns)
  # 
  # meanRho.TTS$position = "TSS"
  # 
  # 
  meanRho.TSS$atg_tss_dist = NA
  meanRho.TTS$atg_tss_dist = NA

  meanRho = rbind(meanRho.ATG, meanRho.TSS, meanRho.TTS)
  
  # Dataset
  meanRho$set = s
  
  if (exists("df_gradient")) {
    df_gradient = rbind(df_gradient, meanRho)
    rm(meanRho)
  } else {
    df_gradient = meanRho
    rm(meanRho)
  }
}
# Save the file ----
saveRDS(df_gradient, file = paste("Data/Recombination/Gradient/meanRho.rds", sep = ""))
rm(df_gradient)






#============================================================================#
# Compute the average gradient - CONTROL No Filtering ----
#============================================================================#


# rm(df_gradient)
# for (s in list_dataset[-7]) {
#   cat(s, "\n")
#   rm(df_pos)
#   rm(df_pos_exons)
#   rm(df_pos_introns)
#   cat("ATG position\n")
#   df_pos = readRDS(paste("Data/Recombination/Gradient/RhoGradient_ATG_", s, ".rds", sep = ""))
#   df_pos$hotoverlap = ifelse(df_pos$hotspot_overlap.filtered.4 > 0, "hot", "cold")
#   
#   # # Remove non genic parts downstream
#   # df_pos = df_pos[-which(df_pos$nb_gene == 0 & df_pos$idx >= 0),]
#   # # In fact, keep only 70% of the gene; do not count rho at the end or after gene's end
#   # p.keep = 0.7
#   # df_pos$relative_pos = df_pos$idx/df_pos$gene_size
#   # df_pos = df_pos[which(df_pos$relative_pos <= p.keep),]
#   
#   meanRho.ATG = aggregate(wmean.rho ~ idx, data = df_pos, mean)
#   meanRho.ATG$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos, length)$wmean.rho
#   meanRho.ATG$condition = "observed"
#   meanRho.ATG$hotoverlap = "both"
#   # meanRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, mean)
#   meanRho.control = aggregate(mean.rho.control ~ idx, data = df_pos, mean)
#   meanRho.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos, length)$mean.rho.control
#   meanRho.control$condition = "control"
#   meanRho.control$hotoverlap = "both"
#   
#   colnames(meanRho.ATG) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
#   colnames(meanRho.control) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
#   meanRho.ATG = rbind(meanRho.ATG, meanRho.control)
#   
#   
#   meanRho.ATG.hot = aggregate(wmean.rho ~ idx + hotoverlap, data = df_pos, mean)
#   meanRho.ATG.hot$n.intervals = aggregate(wmean.rho ~ idx + hotoverlap, data = df_pos, length)$wmean.rho
#   meanRho.ATG.hot$condition = "observed"
#   # meanRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, mean)
#   meanRho.control.hot = aggregate(mean.rho.control ~ idx + hotoverlap, data = df_pos, mean)
#   meanRho.control.hot$n.intervals = aggregate(mean.rho.control ~ idx + hotoverlap, data = df_pos, length)$mean.rho.control
#   meanRho.control.hot$condition = "control"
#   
#   colnames(meanRho.ATG.hot) = c("idx", "hotoverlap",  "meanRho", "n.intervals", "condition")
#   colnames(meanRho.control.hot) = c("idx", "hotoverlap", "meanRho", "n.intervals", "condition")
#   meanRho.ATG = rbind(meanRho.ATG, meanRho.ATG.hot, meanRho.control.hot)
#   
#   # Remove intervals overlapping two gene or more than the gene size limit
#   df_pos_exons = df_pos[which((df_pos$nb_gene == 1) & df_pos$exon_overlap == TRUE),]
#   df_pos_introns = df_pos[which((df_pos$nb_gene == 1) & df_pos$intron_overlap == TRUE),]
#   
#   # mean as a function of position
#   meanRho_exons = aggregate(wmean.rho ~ idx, data = df_pos_exons, mean)
#   meanRho_exons$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_exons, length)$wmean.rho
#   meanRho_exons$condition = "observed"
#   meanRho_exons$hotoverlap = "both"
#   # meanRho_exons.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos_exons, mean)
#   meanRho_exons.control = aggregate(mean.rho.control ~ idx, data = df_pos_exons, mean)
#   meanRho_exons.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_exons, length)$mean.rho.control
#   meanRho_exons.control$condition = "control"
#   meanRho_exons.control$hotoverlap = "both"
#   colnames(meanRho_exons) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
#   colnames(meanRho_exons.control) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
#   meanRho_exons = rbind(meanRho_exons, meanRho_exons.control)
#   
#   
#   meanRho_introns = aggregate(wmean.rho ~ idx, data = df_pos_introns, mean)
#   meanRho_introns$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_introns, length)$wmean.rho
#   meanRho_introns$condition = "observed"
#   meanRho_introns$hotoverlap = "both"
#   # meanRho_introns.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos_introns, mean)
#   meanRho_introns.control = aggregate(mean.rho.control ~ idx, data = df_pos_introns, mean)
#   meanRho_introns.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_introns, length)$mean.rho.control
#   meanRho_introns.control$condition = "control"
#   meanRho_introns.control$hotoverlap = "both"
#   colnames(meanRho_introns) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
#   colnames(meanRho_introns.control) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
#   meanRho_introns = rbind(meanRho_introns, meanRho_introns.control)
#   
#   
#   meanRho.ATG$level = "all"
#   meanRho_exons$level = "exon" 
#   meanRho_introns$level = "intron"
#   
#   meanRho_exons = meanRho_exons[which(meanRho_exons$idx >= 0),]
#   meanRho_introns = meanRho_introns[which(meanRho_introns$idx >= 0),]
#   
#   meanRho.ATG = rbind(meanRho.ATG, meanRho_exons, meanRho_introns)
#   
#   meanRho.ATG$position = "ATG"
#   
#   # Mean ATG-TSS distance of the species
#   atg_tss_dist = mean(df_pos$dist_tss, na.rm = TRUE)
#   meanRho.ATG$atg_tss_dist = atg_tss_dist
#   
#   
#   # ------------------------------ #
#   cat("TSS position\n")
#   df_pos = readRDS(paste("Data/Recombination/Gradient/RhoGradient_5kbTSS_", s, ".rds", sep = ""))
#   df_pos$hotoverlap = ifelse(df_pos$hotspot_overlap.filtered.4 > 0, "hot", "cold")
#   
#   meanRho.TSS = aggregate(wmean.rho ~ idx, data = df_pos, mean)
#   meanRho.TSS$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos, length)$wmean.rho
#   meanRho.TSS$condition = "observed"
#   meanRho.TSS$hotoverlap = "both"
#   # meanRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, mean)
#   meanRho.control = aggregate(mean.rho.control ~ idx, data = df_pos, mean)
#   meanRho.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos, length)$mean.rho.control
#   meanRho.control$condition = "control"
#   meanRho.control$hotoverlap = "both"
#   
#   colnames(meanRho.TSS) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
#   colnames(meanRho.control) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
#   meanRho.TSS = rbind(meanRho.TSS, meanRho.control)
#   
#   
#   meanRho.TSS.hot = aggregate(wmean.rho ~ idx + hotoverlap, data = df_pos, mean)
#   meanRho.TSS.hot$n.intervals = aggregate(wmean.rho ~ idx + hotoverlap, data = df_pos, length)$wmean.rho
#   meanRho.TSS.hot$condition = "observed"
#   # meanRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, mean)
#   meanRho.control.hot = aggregate(mean.rho.control ~ idx + hotoverlap, data = df_pos, mean)
#   meanRho.control.hot$n.intervals = aggregate(mean.rho.control ~ idx + hotoverlap, data = df_pos, length)$mean.rho.control
#   meanRho.control.hot$condition = "control"
#   
#   colnames(meanRho.TSS.hot) = c("idx", "hotoverlap",  "meanRho", "n.intervals", "condition")
#   colnames(meanRho.control.hot) = c("idx", "hotoverlap", "meanRho", "n.intervals", "condition")
#   meanRho.TSS = rbind(meanRho.TSS, meanRho.TSS.hot, meanRho.control.hot)
#   
#   # Remove intervals overlapping two gene or more than the gene size limit
#   df_pos_exons = df_pos[which((df_pos$nb_gene == 1) & df_pos$exon_overlap == TRUE),]
#   df_pos_introns = df_pos[which((df_pos$nb_gene == 1) & df_pos$intron_overlap == TRUE),]
#   
#   # mean as a function of position
#   meanRho_exons = aggregate(wmean.rho ~ idx, data = df_pos_exons, mean)
#   meanRho_exons$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_exons, length)$wmean.rho
#   meanRho_exons$condition = "observed"
#   meanRho_exons$hotoverlap = "both"
#   # meanRho_exons.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos_exons, mean)
#   meanRho_exons.control = aggregate(mean.rho.control ~ idx, data = df_pos_exons, mean)
#   meanRho_exons.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_exons, length)$mean.rho.control
#   meanRho_exons.control$condition = "control"
#   meanRho_exons.control$hotoverlap = "both"
#   colnames(meanRho_exons) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
#   colnames(meanRho_exons.control) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
#   meanRho_exons = rbind(meanRho_exons, meanRho_exons.control)
#   
#   
#   meanRho_introns = aggregate(wmean.rho ~ idx, data = df_pos_introns, mean)
#   meanRho_introns$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_introns, length)$wmean.rho
#   meanRho_introns$condition = "observed"
#   meanRho_introns$hotoverlap = "both"
#   # meanRho_introns.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos_introns, mean)
#   meanRho_introns.control = aggregate(mean.rho.control ~ idx, data = df_pos_introns, mean)
#   meanRho_introns.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_introns, length)$mean.rho.control
#   meanRho_introns.control$condition = "control"
#   meanRho_introns.control$hotoverlap = "both"
#   colnames(meanRho_introns) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
#   colnames(meanRho_introns.control) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
#   meanRho_introns = rbind(meanRho_introns, meanRho_introns.control)
#   
#   
#   meanRho.TSS$level = "all"
#   meanRho_exons$level = "exon" 
#   meanRho_introns$level = "intron"
#   
#   meanRho_exons = meanRho_exons[which(meanRho_exons$idx >= 0),]
#   meanRho_introns = meanRho_introns[which(meanRho_introns$idx >= 0),]
#   
#   meanRho.TSS = rbind(meanRho.TSS, meanRho_exons, meanRho_introns)
#   
#   meanRho.TSS$position = "TSS"
#   
#   # ------------------------------ #
#   cat("TTS position\n")
#   df_pos = readRDS(paste("Data/Recombination/Gradient/RhoGradient_5kbTTS_", s, ".rds", sep = ""))
#   
#   df_pos$hotoverlap = ifelse(df_pos$hotspot_overlap.filtered.4 > 0, "hot", "cold")
#   
#   meanRho.TTS = aggregate(wmean.rho ~ idx, data = df_pos, mean)
#   meanRho.TTS$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos, length)$wmean.rho
#   meanRho.TTS$condition = "observed"
#   meanRho.TTS$hotoverlap = "both"
#   # meanRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, mean)
#   meanRho.control = aggregate(mean.rho.control ~ idx, data = df_pos, mean)
#   meanRho.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos, length)$mean.rho.control
#   meanRho.control$condition = "control"
#   meanRho.control$hotoverlap = "both"
#   
#   colnames(meanRho.TTS) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
#   colnames(meanRho.control) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
#   meanRho.TTS = rbind(meanRho.TTS, meanRho.control)
#   
#   
#   meanRho.TTS.hot = aggregate(wmean.rho ~ idx + hotoverlap, data = df_pos, mean)
#   meanRho.TTS.hot$n.intervals = aggregate(wmean.rho ~ idx + hotoverlap, data = df_pos, length)$wmean.rho
#   meanRho.TTS.hot$condition = "observed"
#   # meanRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, mean)
#   meanRho.control.hot = aggregate(mean.rho.control ~ idx + hotoverlap, data = df_pos, mean)
#   meanRho.control.hot$n.intervals = aggregate(mean.rho.control ~ idx + hotoverlap, data = df_pos, length)$mean.rho.control
#   meanRho.control.hot$condition = "control"
#   
#   colnames(meanRho.TTS.hot) = c("idx", "hotoverlap",  "meanRho", "n.intervals", "condition")
#   colnames(meanRho.control.hot) = c("idx", "hotoverlap", "meanRho", "n.intervals", "condition")
#   meanRho.TTS = rbind(meanRho.TTS, meanRho.TTS.hot, meanRho.control.hot)
#   
#   # Remove intervals overlapping two gene or more than the gene size limit
#   # df_pos_exons = df_pos[which((df_pos$nb_gene == 1) & df_pos$exon_overlap == TRUE),]
#   # df_pos_introns = df_pos[which((df_pos$nb_gene == 1) & df_pos$intron_overlap == TRUE),]
#   df_pos_exons = df_pos[which(df_pos$exon_overlap == TRUE),]
#   df_pos_introns = df_pos[which(df_pos$intron_overlap == TRUE),]
#   
#   # mean as a function of position
#   meanRho_exons = aggregate(wmean.rho ~ idx, data = df_pos_exons, mean)
#   meanRho_exons$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_exons, length)$wmean.rho
#   meanRho_exons$condition = "observed"
#   meanRho_exons$hotoverlap = "both"
#   # meanRho_exons.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos_exons, mean)
#   meanRho_exons.control = aggregate(mean.rho.control ~ idx, data = df_pos_exons, mean)
#   meanRho_exons.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_exons, length)$mean.rho.control
#   meanRho_exons.control$condition = "control"
#   meanRho_exons.control$hotoverlap = "both"
#   colnames(meanRho_exons) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
#   colnames(meanRho_exons.control) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
#   meanRho_exons = rbind(meanRho_exons, meanRho_exons.control)
#   
#   
#   meanRho_introns = aggregate(wmean.rho ~ idx, data = df_pos_introns, mean)
#   meanRho_introns$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_introns, length)$wmean.rho
#   meanRho_introns$condition = "observed"
#   meanRho_introns$hotoverlap = "both"
#   # meanRho_introns.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos_introns, mean)
#   meanRho_introns.control = aggregate(mean.rho.control ~ idx, data = df_pos_introns, mean)
#   meanRho_introns.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_introns, length)$mean.rho.control
#   meanRho_introns.control$condition = "control"
#   meanRho_introns.control$hotoverlap = "both"
#   colnames(meanRho_introns) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
#   colnames(meanRho_introns.control) = c("idx", "meanRho", "n.intervals", "condition", "hotoverlap")
#   meanRho_introns = rbind(meanRho_introns, meanRho_introns.control)
#   
#   
#   meanRho.TTS$level = "all"
#   meanRho_exons$level = "exon" 
#   meanRho_introns$level = "intron"
#   
#   meanRho_exons = meanRho_exons[which(meanRho_exons$idx >= 0),]
#   meanRho_introns = meanRho_introns[which(meanRho_introns$idx >= 0),]
#   
#   meanRho.TTS = rbind(meanRho.TTS, meanRho_exons, meanRho_introns)
#   
#   meanRho.TTS$position = "TSS"
#   
#   
#   
#   meanRho.TSS$atg_tss_dist = NA
#   meanRho.TTS$atg_tss_dist = NA
#   
#   meanRho = rbind(meanRho.ATG, meanRho.TSS, meanRho.TTS)
#   
#   # Dataset
#   meanRho$set = s
#   
#   if (exists("df_gradient")) {
#     df_gradient = rbind(df_gradient, meanRho)
#     rm(meanRho)
#   } else {
#     df_gradient = meanRho
#     rm(meanRho)
#   }
# }
# saveRDS(df_gradient, file = paste("Data/Recombination/Gradient/meanRho_controlNoFiltering.rds", sep = ""))
# rm(df_gradient)

# # Aggregate mean Rho
# # Recombination gradient in exons (bp)
# # interval = 200
# # upstream = 5000 # Size of the upstream region
# # downstream = 5000 # Max size of the downstream region
# rm(df_gradient)
# for (s in list_dataset) {
#   cat(s, "\n")
#   rm(df_pos)
#   rm(df_pos_exons)
#   rm(df_pos_introns)
#   
#   if (file.exists(paste("Data/Recombination/Gradient/RhoGradient_5kbTSS_", s, ".rds", sep = ""))) {
#     df_pos_TSS = readRDS(paste("Data/Recombination/Gradient/RhoGradient_5kbTSS_", s, ".rds", sep = ""))
#     df_pos_TTS = readRDS(paste("Data/Recombination/Gradient/RhoGradient_5kbTTS_", s, ".rds", sep = ""))
#     
#     # Merge both ends with 2kb inside genes
#     df_pos_TSS = df_pos_TSS[which(df_pos_TSS$idx <= 1000),]
#     df_pos_TTS = df_pos_TTS[which(df_pos_TTS$idx >= -1000),]
#     # Change scale
#     df_pos_TSS$idx = df_pos_TSS$idx - 1000
#     df_pos_TTS$idx = df_pos_TTS$idx + 1000
#     
#     df_pos = rbind.fill(df_pos_TSS, df_pos_TTS)
#     # idx 0 count twice; remove it
#     df_pos = df_pos[-which(df_pos$idx == 0),]
#     
#     # To check if strand is well considered
#     # df_pos = df_pos[which(df_pos$strand == "+"),]
#     
#     overlap = 0.7 * interval
#     
#     meanRho = aggregate(wmean.rho ~ idx, data = df_pos, mean)
#     meanRho$condition = "observed"
#     # meanRho.control = aggregate(mean.rho ~ idx_resample, data = df_pos, mean)
#     meanRho.control = aggregate(mean.rho.control ~ idx, data = df_pos, mean)
#     meanRho.control$condition = "control"
#     colnames(meanRho) = c("idx", "meanRho", "condition")
#     colnames(meanRho.control) = c("idx", "meanRho", "condition")
#     meanRho = rbind(meanRho, meanRho.control)
#     
#     
#     df_pos_exons = df_pos[which((df_pos$nb_gene == 1) & df_pos$exon_overlap > overlap),]
#     df_pos_introns = df_pos[which((df_pos$nb_gene == 1) & df_pos$intron_overlap > overlap),]
#     
#     meanRho_exons = aggregate(wmean_exon ~ idx, data = df_pos_exons, mean)
#     meanRho_exons$condition = "observed"
#     # meanRho_exons.control = aggregate(mean.rho ~ idx_resample, data = df_pos_exons, mean)
#     meanRho_exons.control = aggregate(mean.rho.control ~ idx, data = df_pos_exons, mean)
#     meanRho_exons.control$condition = "control"
#     colnames(meanRho_exons) = c("idx", "meanRho", "condition")
#     colnames(meanRho_exons.control) = c("idx", "meanRho", "condition")
#     meanRho_exons = rbind(meanRho_exons, meanRho_exons.control)
#     
#     
#     meanRho_introns = aggregate(wmean_intron ~ idx, data = df_pos_introns, mean)
#     meanRho_introns$condition = "observed"
#     # meanRho_introns.control = aggregate(mean.rho ~ idx_resample, data = df_pos_introns, mean)
#     meanRho_introns.control = aggregate(mean.rho.control ~ idx, data = df_pos_introns, mean)
#     meanRho_introns.control$condition = "control"
#     colnames(meanRho_introns) = c("idx", "meanRho", "condition")
#     colnames(meanRho_introns.control) = c("idx", "meanRho", "condition")
#     meanRho_introns = rbind(meanRho_introns, meanRho_introns.control)
#     
#     
#     meanRho$level = "all" 
#     meanRho_exons$level = "exon" 
#     meanRho_introns$level = "intron"
#     
#     meanRho_exons = meanRho_exons[which(meanRho_exons$idx >= 0),]
#     meanRho_introns = meanRho_introns[which(meanRho_introns$idx >= 0),]
#     
#     meanRho = rbind(meanRho, meanRho_exons, meanRho_introns)
#     
#     # Hotspot count
#     
#     # hotoverlap = function(x) {
#     #   # A function that count the number of hotspot midpoint overlapping a given genomic window
#     #   c = NA
#     #   c = sum((hotspots$start >= x$start & hotspots$start <= x$end & hotspots$seqnames == df_pos$seqnames) | (hotspots$end >= x$start & hotspots$end <= x$end & hotspots$seqnames == df_pos$seqnames), na.rm = TRUE)
#     #   ifelse(length(c) > 0, return(c), return(NA))
#     # }
#     # hotcount = pbmclapply(1:100, function(x) hotoverlap(df_pos[x,]))
#     df_pos_ranges = makeGRangesFromDataFrame(df_pos, keep.extra.columns = TRUE)
#     hotCount = data.frame(idx = unique(df_pos_ranges$idx), hotcount_overlap = NA, hotcount_midpoint = NA)
#     
#     # Overlap hotspot range
#     hotspots = read.table(file = gzfile("Data/Recombination/LD/ldhotspots_filtered.csv.gz"), header = TRUE)
#     hotspots = hotspots[which(hotspots$set == s),]
#     hotspots$midpoint = hotspots$start + (hotspots$length/2)
#     hotspots = GRanges(seqnames = hotspots$seqnames, range = IRanges(start = hotspots$start, end = hotspots$end), strand = "+")
#     # Translate chromosome names
#     chrnames = seqlevels(hotspots)
#     for (j in 1:length(chrnames)) {
#       # cat(chr, '\n')
#       chrnames[j] = chromosome_metadata$annotname[which(chromosome_metadata$ldmapname == chrnames[j] & chromosome_metadata$set == s)]
#     }
#     seqlevels(hotspots) = chrnames
#     for (i in 1:nrow(hotCount)) {
#       hotCount$hotcount_overlap[i] = sum(countOverlaps(df_pos_ranges[df_pos_ranges$idx == hotCount$idx[i]], hotspots), na.rm = TRUE)
#     }
#     
#     # Overlap only midpoint
#     hotspots = read.table(file = gzfile("Data/Recombination/LD/ldhotspots_filtered.csv.gz"), header = TRUE)
#     hotspots = hotspots[which(hotspots$set == s),]
#     hotspots$midpoint = hotspots$start + (hotspots$length/2)
#     hotspots = GRanges(seqnames = hotspots$seqnames, range =IRanges(start = hotspots$midpoint, width = 1), strand = "+")
#     # Translate chromosome names
#     chrnames = seqlevels(hotspots)
#     for (j in 1:length(chrnames)) {
#       # cat(chr, '\n')
#       chrnames[j] = chromosome_metadata$annotname[which(chromosome_metadata$ldmapname == chrnames[j] & chromosome_metadata$set == s)]
#     }
#     seqlevels(hotspots) = chrnames
#     for (i in 1:nrow(hotCount)) {
#       hotCount$hotcount_midpoint[i] = sum(countOverlaps(df_pos_ranges[df_pos_ranges$idx == hotCount$idx[i]], hotspots), na.rm = TRUE)
#     }
#     
#     # # Exons
#     # hotcount_exons = unlist(pbmclapply(1:nrow(df_pos_exons), function(x) hotoverlap(df_pos_exons[x,])))
#     # df_pos_exons$hotcount = hotcount_exons
#     # hotCount_exons = aggregate(hotcount ~ idx, data = df_pos_exons, sum)
#     # # Introns
#     # hotcount_introns = unlist(pbmclapply(1:nrow(df_pos_introns), function(x) hotoverlap(df_pos_introns[x,])))
#     # df_pos_introns$hotcount = hotcount_introns
#     # hotCount_introns = aggregate(hotcount ~ idx, data = df_pos_introns, sum)
#     
#     meanRho$hotspot_overlap = NA
#     meanRho$hotspot_overlap[which(meanRho$condition == "observed" & meanRho$level == "all")] = hotCount$hotcount_overlap
#     meanRho$hotspot_count = NA
#     meanRho$hotspot_count[which(meanRho$condition == "observed" & meanRho$level == "all")] = hotCount$hotcount_midpoint
#     meanRho$hotspot_relative = NA
#     meanRho$hotspot_relative[which(meanRho$condition == "observed" & meanRho$level == "all")] = hotCount$hotcount_midpoint/sum(hotCount$hotcount, na.rm = TRUE)
#     # meanRho$hotspot_count[which(meanRho$condition == "observed" & meanRho$level == "exon")] = hotCount_exons$hotcount[which(hotCount_exons$idx >= 0)]
#     # meanRho$hotspot_count[which(meanRho$condition == "observed" & meanRho$level == "intron")] = hotCount_introns$hotcount[which(hotCount_exons$idx >= 0)]
#     
#     # Set
#     meanRho$set = s
#     
#     # TODO Number of data
#     
#     
#     # Mean ATG-TSS distance of the species
#     atg_tss_dist = mean(df_pos$dist_atg, na.rm = TRUE)
#     meanRho$atg_tss_dist = atg_tss_dist
#    
#     if (exists("df_gradient")) {
#       df_gradient = rbind(df_gradient, meanRho)
#       rm(meanRho)
#     } else {
#       df_gradient = meanRho
#       rm(meanRho)
#     }
#   }
# }
# # df_gradient
# 
# saveRDS(df_gradient, file = paste("Data/Recombination/Gradient/gradient_5kbTSS_TTS.rds", sep = ""))
# rm(df_gradient)
# 
# 
# 
# ### Gradient TTS ----
# # Aggregate mean Rho
# # Recombination gradient in exons (bp)
# # interval = 200
# # upstream = 2000 # Size of the upstream region
# # downstream = 2000 # Max size of the downstream region
# rm(df_gradient)
# for (s in list_dataset) {
#   cat(s, "\n")
#   rm(df_pos)
#   rm(df_pos_exons)
#   rm(df_pos_introns)
#   df_pos_TTS = readRDS(paste("Data/Recombination/Gradient/RhoGradient_5kbTTS_", s, ".rds", sep = ""))
#   
#   # Rescale
#   df_pos_TTS = df_pos_TTS[which(df_pos_TTS$idx >= -upstream & df_pos_TTS$idx <= downstream),]
#   
#   df_pos = df_pos_TTS
# 
#   # To check if strand is well considered
#   # df_pos = df_pos[which(df_pos$strand == "+"),]
#   
#   overlap = 0.7 * interval
#   
#   meanRho = aggregate(wmean.rho ~ idx, data = df_pos, mean)
#   meanRho$condition = "observed"
#   # meanRho.control = aggregate(mean.rho ~ idx_resample, data = df_pos, mean)
#   meanRho.control = aggregate(mean.rho.control ~ idx, data = df_pos, mean)
#   meanRho.control$condition = "control"
#   colnames(meanRho) = c("idx", "meanRho", "condition")
#   colnames(meanRho.control) = c("idx", "meanRho", "condition")
#   meanRho = rbind(meanRho, meanRho.control)
#   
#   meanRho$level = "all" 
#   df_pos_ranges = makeGRangesFromDataFrame(df_pos, keep.extra.columns = TRUE)
#   hotCount = data.frame(idx = unique(df_pos_ranges$idx), hotcount_overlap = NA, hotcount_midpoint = NA)
#   
#   # Overlap hotspot range
#   hotspots = read.table(file = gzfile("Data/Recombination/LD/ldhotspots_filtered.csv.gz"), header = TRUE)
#   hotspots = hotspots[which(hotspots$set == s),]
#   hotspots$midpoint = hotspots$start + (hotspots$length/2)
#   hotspots = GRanges(seqnames = hotspots$seqnames, range = IRanges(start = hotspots$start, end = hotspots$end), strand = "+")
#   # Translate chromosome names
#   chrnames = seqlevels(hotspots)
#   for (j in 1:length(chrnames)) {
#     # cat(chr, '\n')
#     chrnames[j] = chromosome_metadata$annotname[which(chromosome_metadata$ldmapname == chrnames[j] & chromosome_metadata$set == s)]
#   }
#   seqlevels(hotspots) = chrnames
#   for (i in 1:nrow(hotCount)) {
#     hotCount$hotcount_overlap[i] = sum(countOverlaps(df_pos_ranges[df_pos_ranges$idx == hotCount$idx[i]], hotspots), na.rm = TRUE)
#   }
#   
#   # Overlap only midpoint
#   hotspots = read.table(file = gzfile("Data/Recombination/LD/ldhotspots_filtered.csv.gz"), header = TRUE)
#   hotspots = hotspots[which(hotspots$set == s),]
#   hotspots$midpoint = hotspots$start + (hotspots$length/2)
#   hotspots = GRanges(seqnames = hotspots$seqnames, range =IRanges(start = hotspots$midpoint, width = 1), strand = "+")
#   # Translate chromosome names
#   chrnames = seqlevels(hotspots)
#   for (j in 1:length(chrnames)) {
#     # cat(chr, '\n')
#     chrnames[j] = chromosome_metadata$annotname[which(chromosome_metadata$ldmapname == chrnames[j] & chromosome_metadata$set == s)]
#   }
#   seqlevels(hotspots) = chrnames
#   for (i in 1:nrow(hotCount)) {
#     hotCount$hotcount_midpoint[i] = sum(countOverlaps(df_pos_ranges[df_pos_ranges$idx == hotCount$idx[i]], hotspots), na.rm = TRUE)
#   }
# 
#   
#   meanRho$hotspot_overlap = NA
#   meanRho$hotspot_overlap[which(meanRho$condition == "observed" & meanRho$level == "all")] = hotCount$hotcount_overlap
#   meanRho$hotspot_count = NA
#   meanRho$hotspot_count[which(meanRho$condition == "observed" & meanRho$level == "all")] = hotCount$hotcount_midpoint
#   meanRho$hotspot_relative = NA
#   meanRho$hotspot_relative[which(meanRho$condition == "observed" & meanRho$level == "all")] = hotCount$hotcount_midpoint/sum(hotCount$hotcount, na.rm = TRUE)
#   # meanRho$hotspot_count[which(meanRho$condition == "observed" & meanRho$level == "exon")] = hotCount_exons$hotcount[which(hotCount_exons$idx >= 0)]
#   # meanRho$hotspot_count[which(meanRho$condition == "observed" & meanRho$level == "intron")] = hotCount_introns$hotcount[which(hotCount_exons$idx >= 0)]
#   
#   # Set
#   meanRho$set = s
#   
#   # TODO Number of data
# 
#   
#   if (exists("df_gradient")) {
#     df_gradient = rbind(df_gradient, meanRho)
#     rm(meanRho)
#   } else {
#     df_gradient = meanRho
#     rm(meanRho)
#   }
# }
# # df_gradient
# 
# saveRDS(df_gradient, file = paste("Data/Recombination/Gradient/gradient_TTS.rds", sep = ""))
# rm(df_gradient)


#============================================================================#
# Rho gradient 5kb TSS + TTS
#============================================================================#
# Aggregate mean Rho
# Recombination gradient in exons (bp)
# interval = 100
# upstream = 2000 # Size of the upstream region
# genic = 2000 # Max size of the genic region
# rm(df_gradient)
# for (s in list_dataset) {
#   cat(s, "\n")
#   rm(df_pos)
#   rm(df_pos_exons)
#   rm(df_pos_introns)
#   df_pos = readRDS(paste("Data/Recombination/Gradient/RhoGradient_ATG_", s, ".rds", sep = ""))
#   # To check if strand is well considered
#   # df_pos = df_pos[which(df_pos$strand == "+"),]
#   
#   overlap = 0.7 * interval
#   
#   # TODO Random control by resampling idx and gene
#   # df_pos$idx_resample = sample(df_pos$idx, replace = FALSE)
#   # df_pos$gene_id_resample = sample(df_pos$gene_id, replace = TRUE)
#   
#   
#   meanRho = aggregate(wmean.rho ~ idx, data = df_pos, mean)
#   meanRho$condition = "observed"
#   # meanRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, mean)
#   meanRho.control = aggregate(mean.rho.control ~ idx, data = df_pos, mean)
#   meanRho.control$condition = "control"
#   colnames(meanRho) = c("idx", "meanRho", "condition")
#   colnames(meanRho.control) = c("idx", "meanRho", "condition")
#   meanRho = rbind(meanRho, meanRho.control)
#   
#   # Remove intervals overlapping two gene or more than the gene size limit
#   df_pos_exons = df_pos[which((df_pos$nb_gene == 1) & df_pos$exon_overlap > overlap),]
#   df_pos_introns = df_pos[which((df_pos$nb_gene == 1) & df_pos$intron_overlap > overlap),]
#   
#   # mean as a function of position
#   meanRho_exons = aggregate(wmean_exon ~ idx, data = df_pos_exons, mean)
#   meanRho_exons$condition = "observed"
#   # meanRho_exons.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos_exons, mean)
#   meanRho_exons.control = aggregate(mean.rho.control ~ idx, data = df_pos_exons, mean)
#   meanRho_exons.control$condition = "control"
#   colnames(meanRho_exons) = c("idx", "meanRho", "condition")
#   colnames(meanRho_exons.control) = c("idx", "meanRho", "condition")
#   meanRho_exons = rbind(meanRho_exons, meanRho_exons.control)
#       
#   
#   meanRho_introns = aggregate(wmean_intron ~ idx, data = df_pos_introns, mean)
#   meanRho_introns$condition = "observed"
#   # meanRho_introns.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos_introns, mean)
#   meanRho_introns.control = aggregate(mean.rho.control ~ idx, data = df_pos_introns, mean)
#   meanRho_introns.control$condition = "control"
#   colnames(meanRho_introns) = c("idx", "meanRho", "condition")
#   colnames(meanRho_introns.control) = c("idx", "meanRho", "condition")
#   meanRho_introns = rbind(meanRho_introns, meanRho_introns.control)
#   
#   
#   meanRho$level = "all" 
#   meanRho_exons$level = "exon" 
#   meanRho_introns$level = "intron"
#   
#   meanRho_exons = meanRho_exons[which(meanRho_exons$idx >= 0),]
#   meanRho_introns = meanRho_introns[which(meanRho_introns$idx >= 0),]
#   
#   meanRho = rbind(meanRho, meanRho_exons, meanRho_introns)
#   
#   # Hotspot count
#   hotspots = read.table(file = gzfile("Data/Recombination/LD/ldhotspots_filtered.csv.gz"), header = TRUE)
#   # Filter hotspots
#   hotspots$length = (hotspots$end - hotspots$start + 1) * 1000
#   # maxhotspotlength = 6000
#   hotspots = hotspots[which(hotspots$set == s),]
#   hotspots$midpoint = hotspots$start + (hotspots$length/2)
#   
#   # Range overlap
#   df_ranges = makeGRangesFromDataFrame(df_pos)
#   hotspots_ranges = GRanges(seqnames = hotspots$seqnames,
#                             IRanges(start = hotspots$start,
#                                     end = hotspots$end),
#                             strand = "*")
#   # Translate chromosome names
#   chrnames = seqlevels(hotspots_ranges)
#   for (j in 1:length(chrnames)) {
#     # cat(chr, '\n')
#     chrnames[j] = chromosome_metadata$annotname[which(chromosome_metadata$ldmapname == chrnames[j] & chromosome_metadata$set == s)]
#   }
#   seqlevels(hotspots_ranges) = chrnames
# 
#   hotoverlap = countOverlaps(df_ranges, hotspots_ranges)
#   df_pos$hotoverlap = hotoverlap
#   
#   # Midpoint overlap
#   hotoverlap = function(x) {
#     # A function that count the number of hotspot midpoint overlapping a given genomic window
#     c = NA
#     c = sum((hotspots$start >= x$start & hotspots$start <= x$end) | (hotspots$end >= x$start & hotspots$end <= x$end), na.rm = TRUE)
#     ifelse(length(c) > 0, return(c), return(NA))
#   }
#   # hotcount = pbmclapply(1:100, function(x) hotoverlap(df_pos[x,]))
#   hotcount = unlist(pbmclapply(1:nrow(df_pos), function(x) hotoverlap(df_pos[x,])))
#   df_pos$hotcount = hotcount
#   hotCount = aggregate(hotcount ~ idx, data = df_pos, sum)
#   # Exons
#   hotcount_exons = unlist(pbmclapply(1:nrow(df_pos_exons), function(x) hotoverlap(df_pos_exons[x,])))
#   df_pos_exons$hotcount = hotcount_exons
#   hotCount_exons = aggregate(hotcount ~ idx, data = df_pos_exons, sum)
#   # Introns
#   hotcount_introns = unlist(pbmclapply(1:nrow(df_pos_introns), function(x) hotoverlap(df_pos_introns[x,])))
#   df_pos_introns$hotcount = hotcount_introns
#   hotCount_introns = aggregate(hotcount ~ idx, data = df_pos_introns, sum)
#   
#   meanRho$hotspot_count = NA
#   meanRho$hotspot_count[which(meanRho$condition == "observed" & meanRho$level == "all")] = hotCount$hotcount
#   meanRho$hotspot_count[which(meanRho$condition == "observed" & meanRho$level == "exon")] = hotCount_exons$hotcount[which(hotCount_exons$idx >= 0)]
#   meanRho$hotspot_count[which(meanRho$condition == "observed" & meanRho$level == "intron")] = hotCount_introns$hotcount[which(hotCount_exons$idx >= 0)]
#   
#   # Set
#   meanRho$set = s
#   
#   # TODO Number of data
#   
#   
#   # Mean ATG-TSS distance of the species
#   atg_tss_dist = mean(df_pos$dist_tss, na.rm = TRUE)
#   meanRho$atg_tss_dist = atg_tss_dist
#   
#   
#   if (exists("df_gradient")) {
#     df_gradient = rbind(df_gradient, meanRho)
#     rm(meanRho)
#   } else {
#     df_gradient = meanRho
#     rm(meanRho)
#   }
# }
# saveRDS(df_gradient, file = paste("Data/Recombination/Gradient/meanRho.rds", sep = ""))
# rm(df_gradient)
# gc()

