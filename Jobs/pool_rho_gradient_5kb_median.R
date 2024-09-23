#!/usr/bin/env Rscript

source("Source/init.R")

interval = 200
upstream = 5000 # Size of the upstream region
downstream = 5000 # Max size of the downstream region



#============================================================================#
# Compute the average gradient MEDIAN instead of mean ----
#============================================================================#
cat("======================================\n")
cat("Compute the average gradient: MEDIAN ")

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
  
  # median Rho over all intervals
  medianRho.ATG = aggregate(wmean.rho ~ idx, data = df_pos, median)
  medianRho.ATG$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos, length)$wmean.rho
  medianRho.ATG$condition = "observed"
  medianRho.ATG$hotoverlap = "both"
  m.obs = median(medianRho.ATG$wmean.rho, na.rm = TRUE)
  sd.obs = sd(medianRho.ATG$wmean.rho, na.rm = TRUE)
  
  # medianRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, median)
  medianRho.control = aggregate(mean.rho.control ~ idx, data = df_pos, median)
  medianRho.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos, length)$mean.rho.control
  medianRho.control$condition = "control"
  medianRho.control$hotoverlap = "both"
  m.control = median(medianRho.control$mean.rho.control, na.rm = TRUE)
  sd.control = sd(medianRho.control$mean.rho.control, na.rm = TRUE)
  
  colnames(medianRho.ATG) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  colnames(medianRho.control) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  
  medianRho.ATG = rbind(medianRho.ATG, medianRho.control)
  medianRho.ATG$m.obs = m.obs
  medianRho.ATG$sd.obs = sd.obs
  medianRho.ATG$m.control = m.control
  medianRho.ATG$sd.control = sd.control
  
  medianRho.ATG$medianRho.rescaled = medianRho.ATG$medianRho
    
  medianRho.ATG$medianRho.rescaled[which(medianRho.ATG$condition == "control")] = (medianRho.control$medianRho - m.control) + m.obs
  
  # HOT vs COLD ----
  cat("Discriminate between hot and cold genes\n")
  
  df_pos$hotoverlap = ifelse(df_pos$hotspot_overlap.filtered.4 > 0, "hot", "cold")
  
  # median Rho over all intervals
  medianRho.ATG.hot = aggregate(wmean.rho ~ idx + hotoverlap, data = df_pos, median)
  medianRho.ATG.hot$n.intervals = aggregate(wmean.rho ~ idx + hotoverlap, data = df_pos, length)$wmean.rho
  medianRho.ATG.hot$condition = "observed"
  m.obs = median(medianRho.ATG.hot$wmean.rho, na.rm = TRUE)
  sd.obs = sd(medianRho.ATG.hot$wmean.rho, na.rm = TRUE)
  
  # medianRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, median)
  medianRho.control.hot = aggregate(mean.rho.control ~ idx + hotoverlap, data = df_pos, median)
  medianRho.control.hot$n.intervals = aggregate(mean.rho.control ~ idx + hotoverlap, data = df_pos, length)$mean.rho.control
  medianRho.control.hot$condition = "control"
  m.control = median(medianRho.control.hot$mean.rho.control, na.rm = TRUE)
  sd.control = sd(medianRho.control.hot$mean.rho.control, na.rm = TRUE)
  
  colnames(medianRho.ATG.hot) = c("idx", "hotoverlap",  "medianRho", "n.intervals", "condition")
  colnames(medianRho.control.hot) = c("idx", "hotoverlap", "medianRho", "n.intervals", "condition")
  
  medianRho.ATG.hot = rbind(medianRho.ATG.hot, medianRho.control.hot)
  medianRho.ATG.hot$m.obs = m.obs
  medianRho.ATG.hot$sd.obs = sd.obs
  medianRho.ATG.hot$m.control = m.control
  medianRho.ATG.hot$sd.control = sd.control
  
  medianRho.ATG.hot$medianRho.rescaled = medianRho.ATG.hot$medianRho
  
  medianRho.ATG.hot$medianRho.rescaled[which(medianRho.ATG.hot$condition == "control")] = (medianRho.control.hot$medianRho - m.control) + m.obs
  
  medianRho.ATG = rbind(medianRho.ATG, medianRho.ATG.hot)
  
  medianRho.ATG$level = "Genic"
  medianRho.ATG$position = "ATG"
  
  # Add exons ----
  cat("Add exons\n")
  df_pos_exons = df_pos[which(df_pos$exon_overlap == TRUE & df_pos$idx > 0),]
  # df_pos_introns = df_pos[which((df_pos$nb_gene == 1) & df_pos$intron_overlap == TRUE),]
  medianRho_exons = aggregate(wmean.rho ~ idx, data = df_pos_exons, median)
  medianRho_exons$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_exons, length)$wmean.rho
  medianRho_exons$condition = "observed"
  medianRho_exons$hotoverlap = "both"
  m.obs = median(medianRho_exons$wmean.rho, na.rm = TRUE)
  sd.obs = sd(medianRho_exons$wmean.rho, na.rm = TRUE)
  
  # medianRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, median)
  medianRho.control_exons = aggregate(mean.rho.control ~ idx, data = df_pos_exons, median)
  medianRho.control_exons$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_exons, length)$mean.rho.control
  medianRho.control_exons$condition = "control"
  medianRho.control_exons$hotoverlap = "both"
  m.control = median(medianRho.control_exons$mean.rho.control, na.rm = TRUE)
  sd.control = sd(medianRho.control_exons$mean.rho.control, na.rm = TRUE)
  
  colnames(medianRho_exons) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  colnames(medianRho.control_exons) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  
  medianRho_exons = rbind(medianRho_exons, medianRho.control_exons)
  medianRho_exons$m.obs = m.obs
  medianRho_exons$sd.obs = sd.obs
  medianRho_exons$m.control = m.control
  medianRho_exons$sd.control = sd.control
  
  medianRho_exons$medianRho.rescaled = medianRho_exons$medianRho
  
  medianRho_exons$medianRho.rescaled[which(medianRho_exons$condition == "control")] = (medianRho.control_exons$medianRho - m.control) + m.obs
  
  medianRho_exons$level = "Exons"
  medianRho_exons$position = "ATG"
  
  medianRho.ATG = rbind(medianRho.ATG, medianRho_exons)
  
  # Add introns ----
  cat("Add introns\n")
  df_pos_introns = df_pos[which(df_pos$intron_overlap == TRUE & df_pos$idx > 0),]

  medianRho_introns = aggregate(wmean.rho ~ idx, data = df_pos_introns, median)
  medianRho_introns$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_introns, length)$wmean.rho
  medianRho_introns$condition = "observed"
  medianRho_introns$hotoverlap = "both"
  m.obs = median(medianRho_introns$wmean.rho, na.rm = TRUE)
  sd.obs = sd(medianRho_introns$wmean.rho, na.rm = TRUE)
  
  # medianRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, median)
  medianRho.control_introns = aggregate(mean.rho.control ~ idx, data = df_pos_introns, median)
  medianRho.control_introns$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_introns, length)$mean.rho.control
  medianRho.control_introns$condition = "control"
  medianRho.control_introns$hotoverlap = "both"
  m.control = median(medianRho.control_introns$mean.rho.control, na.rm = TRUE)
  sd.control = sd(medianRho.control_introns$mean.rho.control, na.rm = TRUE)
  
  colnames(medianRho_introns) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  colnames(medianRho.control_introns) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  
  medianRho_introns = rbind(medianRho_introns, medianRho.control_introns)
  medianRho_introns$m.obs = m.obs
  medianRho_introns$sd.obs = sd.obs
  medianRho_introns$m.control = m.control
  medianRho_introns$sd.control = sd.control
  
  medianRho_introns$medianRho.rescaled = medianRho_introns$medianRho
  
  medianRho_introns$medianRho.rescaled[which(medianRho_introns$condition == "control")] = (medianRho.control_introns$medianRho - m.control) + m.obs
  
  medianRho_introns$level = "Introns"
  medianRho_introns$position = "ATG"
  
  medianRho.ATG = rbind(medianRho.ATG, medianRho_introns)
  
  # df_pos_exons = df_pos[which((df_pos$nb_gene == 1) & df_pos$exon_overlap == TRUE),]
  # df_pos_introns = df_pos[which((df_pos$nb_gene == 1) & df_pos$intron_overlap == TRUE),]
  # 
  # # median as a function of position
  # medianRho_exons = aggregate(wmean.rho ~ idx, data = df_pos_exons, median)
  # medianRho_exons$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_exons, length)$wmean.rho
  # medianRho_exons$condition = "observed"
  # medianRho_exons$hotoverlap = "both"
  # # medianRho_exons.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos_exons, median)
  # medianRho_exons.control = aggregate(mean.rho.control ~ idx, data = df_pos_exons, median)
  # medianRho_exons.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_exons, length)$mean.rho.control
  # medianRho_exons.control$condition = "control"
  # medianRho_exons.control$hotoverlap = "both"
  # colnames(medianRho_exons) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  # colnames(medianRho_exons.control) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  # medianRho_exons = rbind(medianRho_exons, medianRho_exons.control)
  # 
  # 
  # medianRho_introns = aggregate(wmean.rho ~ idx, data = df_pos_introns, median)
  # medianRho_introns$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_introns, length)$wmean.rho
  # medianRho_introns$condition = "observed"
  # medianRho_introns$hotoverlap = "both"
  # # medianRho_introns.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos_introns, median)
  # medianRho_introns.control = aggregate(mean.rho.control ~ idx, data = df_pos_introns, median)
  # medianRho_introns.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_introns, length)$mean.rho.control
  # medianRho_introns.control$condition = "control"
  # medianRho_introns.control$hotoverlap = "both"
  # colnames(medianRho_introns) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  # colnames(medianRho_introns.control) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  # medianRho_introns = rbind(medianRho_introns, medianRho_introns.control)
  # 
  # 
  # medianRho.ATG$level = "all"
  # medianRho_exons$level = "exon" 
  # medianRho_introns$level = "intron"
  # 
  # medianRho_exons = medianRho_exons[which(medianRho_exons$idx >= 0),]
  # medianRho_introns = medianRho_introns[which(medianRho_introns$idx >= 0),]
  # 
  # medianRho.ATG = rbind(medianRho.ATG, medianRho_exons, medianRho_introns)
  # 
  # medianRho.ATG$position = "ATG"
  # 
  # # Mean ATG-TSS distance of the species
  atg_tss_dist = mean(df_pos$dist_tss, na.rm = TRUE)
  medianRho.ATG$atg_tss_dist = atg_tss_dist
  # 
  # 
  # # ------------------------------ #
  # TSS ----
  cat("TSS position\n")
  df_pos = readRDS(paste("Data/Recombination/Gradient/RhoGradient_5kbTSS_", s, ".rds", sep = ""))
  
  # Keep only genes with less than 15 exons
  list_genes = data_all$gene_id[which(data_all$nb_exons <= max.exons)]
  df_pos = df_pos[which(df_pos$gene_id %in% list_genes),]
  
  
  medianRho.TSS = aggregate(wmean.rho ~ idx, data = df_pos,median)
  medianRho.TSS$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos, length)$wmean.rho
  medianRho.TSS$condition = "observed"
  medianRho.TSS$hotoverlap = "both"
  m.obs =median(medianRho.TSS$wmean.rho, na.rm = TRUE)
  sd.obs = sd(medianRho.TSS$wmean.rho, na.rm = TRUE)
  
  # medianRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos,median)
  medianRho.control = aggregate(mean.rho.control ~ idx, data = df_pos,median)
  medianRho.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos, length)$mean.rho.control
  medianRho.control$condition = "control"
  medianRho.control$hotoverlap = "both"
  m.control =median(medianRho.control$mean.rho.control, na.rm = TRUE)
  sd.control = sd(medianRho.control$mean.rho.control, na.rm = TRUE)
  
  colnames(medianRho.TSS) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  colnames(medianRho.control) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  
  medianRho.TSS = rbind(medianRho.TSS, medianRho.control)
  medianRho.TSS$m.obs = m.obs
  medianRho.TSS$sd.obs = sd.obs
  medianRho.TSS$m.control = m.control
  medianRho.TSS$sd.control = sd.control
  
  medianRho.TSS$medianRho.rescaled = medianRho.TSS$medianRho
  
  medianRho.TSS$medianRho.rescaled[which(medianRho.TSS$condition == "control")] = (medianRho.control$medianRho - m.control) + m.obs
  
  medianRho.TSS$position = "TSS"
  medianRho.TSS$level = "Genic"
  
  # df_pos$hotoverlap = ifelse(df_pos$hotspot_overlap.filtered.4 > 0, "hot", "cold")
  # 
  # medianRho.TSS = aggregate(wmean.rho ~ idx, data = df_pos,median)
  # medianRho.TSS$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos, length)$wmean.rho
  # medianRho.TSS$condition = "observed"
  # medianRho.TSS$hotoverlap = "both"
  # # medianRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos,median)
  # medianRho.control = aggregate(mean.rho.control ~ idx, data = df_pos,median)
  # medianRho.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos, length)$mean.rho.control
  # medianRho.control$condition = "control"
  # medianRho.control$hotoverlap = "both"
  # 
  # colnames(medianRho.TSS) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  # colnames(medianRho.control) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  # medianRho.TSS = rbind(medianRho.TSS, medianRho.control)
  # 
  # 
  # medianRho.TSS.hot = aggregate(wmean.rho ~ idx + hotoverlap, data = df_pos,median)
  # medianRho.TSS.hot$n.intervals = aggregate(wmean.rho ~ idx + hotoverlap, data = df_pos, length)$wmean.rho
  # medianRho.TSS.hot$condition = "observed"
  # # medianRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos,median)
  # medianRho.control.hot = aggregate(mean.rho.control ~ idx + hotoverlap, data = df_pos,median)
  # medianRho.control.hot$n.intervals = aggregate(mean.rho.control ~ idx + hotoverlap, data = df_pos, length)$mean.rho.control
  # medianRho.control.hot$condition = "control"
  # 
  # colnames(medianRho.TSS.hot) = c("idx", "hotoverlap",  "medianRho", "n.intervals", "condition")
  # colnames(medianRho.control.hot) = c("idx", "hotoverlap", "medianRho", "n.intervals", "condition")
  # medianRho.TSS = rbind(medianRho.TSS, medianRho.TSS.hot, medianRho.control.hot)
  # 
  # # Remove intervals overlapping two gene or more than the gene size limit
  # df_pos_exons = df_pos[which((df_pos$nb_gene == 1) & df_pos$exon_overlap == TRUE),]
  # df_pos_introns = df_pos[which((df_pos$nb_gene == 1) & df_pos$intron_overlap == TRUE),]
  # 
  # #median as a function of position
  # medianRho_exons = aggregate(wmean.rho ~ idx, data = df_pos_exons,median)
  # medianRho_exons$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_exons, length)$wmean.rho
  # medianRho_exons$condition = "observed"
  # medianRho_exons$hotoverlap = "both"
  # # medianRho_exons.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos_exons,median)
  # medianRho_exons.control = aggregate(mean.rho.control ~ idx, data = df_pos_exons,median)
  # medianRho_exons.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_exons, length)$mean.rho.control
  # medianRho_exons.control$condition = "control"
  # medianRho_exons.control$hotoverlap = "both"
  # colnames(medianRho_exons) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  # colnames(medianRho_exons.control) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  # medianRho_exons = rbind(medianRho_exons, medianRho_exons.control)
  # 
  # 
  # medianRho_introns = aggregate(wmean.rho ~ idx, data = df_pos_introns,median)
  # medianRho_introns$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_introns, length)$wmean.rho
  # medianRho_introns$condition = "observed"
  # medianRho_introns$hotoverlap = "both"
  # # medianRho_introns.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos_introns,median)
  # medianRho_introns.control = aggregate(mean.rho.control ~ idx, data = df_pos_introns,median)
  # medianRho_introns.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_introns, length)$mean.rho.control
  # medianRho_introns.control$condition = "control"
  # medianRho_introns.control$hotoverlap = "both"
  # colnames(medianRho_introns) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  # colnames(medianRho_introns.control) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  # medianRho_introns = rbind(medianRho_introns, medianRho_introns.control)
  # 
  # 
  # medianRho.TSS$level = "all"
  # medianRho_exons$level = "exon" 
  # medianRho_introns$level = "intron"
  # 
  # medianRho_exons = medianRho_exons[which(medianRho_exons$idx >= 0),]
  # medianRho_introns = medianRho_introns[which(medianRho_introns$idx >= 0),]
  # 
  # medianRho.TSS = rbind(medianRho.TSS, medianRho_exons, medianRho_introns)
  # 

  # 
  # # ------------------------------ #
  # TTS ----
  cat("TTS position\n")
  df_pos = readRDS(paste("Data/Recombination/Gradient/RhoGradient_5kbTTS_", s, ".rds", sep = ""))
  
  # Keep only genes with less than 15 exons
  list_genes = data_all$gene_id[which(data_all$nb_exons <= max.exons)]
  df_pos = df_pos[which(df_pos$gene_id %in% list_genes),]
  
  medianRho.TTS = aggregate(wmean.rho ~ idx, data = df_pos,median)
  medianRho.TTS$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos, length)$wmean.rho
  medianRho.TTS$condition = "observed"
  medianRho.TTS$hotoverlap = "both"
  m.obs =median(medianRho.TTS$wmean.rho, na.rm = TRUE)
  sd.obs = sd(medianRho.TTS$wmean.rho, na.rm = TRUE)
  
  # medianRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos,median)
  medianRho.control = aggregate(mean.rho.control ~ idx, data = df_pos, median)
  medianRho.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos, length)$mean.rho.control
  medianRho.control$condition = "control"
  medianRho.control$hotoverlap = "both"
  m.control = median(medianRho.control$mean.rho.control, na.rm = TRUE)
  sd.control = sd(medianRho.control$mean.rho.control, na.rm = TRUE)
  
  colnames(medianRho.TTS) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  colnames(medianRho.control) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  
  medianRho.TTS = rbind(medianRho.TTS, medianRho.control)
  medianRho.TTS$m.obs = m.obs
  medianRho.TTS$sd.obs = sd.obs
  medianRho.TTS$m.control = m.control
  medianRho.TTS$sd.control = sd.control
  
  medianRho.TTS$medianRho.rescaled = medianRho.TTS$medianRho
  
  medianRho.TTS$medianRho.rescaled[which(medianRho.TTS$condition == "control")] = (medianRho.control$medianRho - m.control) + m.obs
  
  medianRho.TTS$position = "TTS"
  medianRho.TTS$level = "Genic"
  
  # df_pos$hotoverlap = ifelse(df_pos$hotspot_overlap.filtered.4 > 0, "hot", "cold")
  # 
  # medianRho.TTS = aggregate(wmean.rho ~ idx, data = df_pos, median)
  # medianRho.TTS$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos, length)$wmean.rho
  # medianRho.TTS$condition = "observed"
  # medianRho.TTS$hotoverlap = "both"
  # # medianRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, median)
  # medianRho.control = aggregate(mean.rho.control ~ idx, data = df_pos, median)
  # medianRho.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos, length)$mean.rho.control
  # medianRho.control$condition = "control"
  # medianRho.control$hotoverlap = "both"
  # 
  # colnames(medianRho.TTS) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  # colnames(medianRho.control) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  # medianRho.TTS = rbind(medianRho.TTS, medianRho.control)
  # 
  # 
  # medianRho.TTS.hot = aggregate(wmean.rho ~ idx + hotoverlap, data = df_pos, median)
  # medianRho.TTS.hot$n.intervals = aggregate(wmean.rho ~ idx + hotoverlap, data = df_pos, length)$wmean.rho
  # medianRho.TTS.hot$condition = "observed"
  # # medianRho.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos, median)
  # medianRho.control.hot = aggregate(mean.rho.control ~ idx + hotoverlap, data = df_pos, median)
  # medianRho.control.hot$n.intervals = aggregate(mean.rho.control ~ idx + hotoverlap, data = df_pos, length)$mean.rho.control
  # medianRho.control.hot$condition = "control"
  # 
  # colnames(medianRho.TTS.hot) = c("idx", "hotoverlap",  "medianRho", "n.intervals", "condition")
  # colnames(medianRho.control.hot) = c("idx", "hotoverlap", "medianRho", "n.intervals", "condition")
  # medianRho.TTS = rbind(medianRho.TTS, medianRho.TTS.hot, medianRho.control.hot)
  # 
  # # Remove intervals overlapping two gene or more than the gene size limit
  # # df_pos_exons = df_pos[which((df_pos$nb_gene == 1) & df_pos$exon_overlap == TRUE),]
  # # df_pos_introns = df_pos[which((df_pos$nb_gene == 1) & df_pos$intron_overlap == TRUE),]
  # df_pos_exons = df_pos[which(df_pos$exon_overlap == TRUE),]
  # df_pos_introns = df_pos[which(df_pos$intron_overlap == TRUE),]
  # 
  # # median as a function of position
  # medianRho_exons = aggregate(wmean.rho ~ idx, data = df_pos_exons, median)
  # medianRho_exons$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_exons, length)$wmean.rho
  # medianRho_exons$condition = "observed"
  # medianRho_exons$hotoverlap = "both"
  # # medianRho_exons.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos_exons, median)
  # medianRho_exons.control = aggregate(mean.rho.control ~ idx, data = df_pos_exons, median)
  # medianRho_exons.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_exons, length)$mean.rho.control
  # medianRho_exons.control$condition = "control"
  # medianRho_exons.control$hotoverlap = "both"
  # colnames(medianRho_exons) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  # colnames(medianRho_exons.control) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  # medianRho_exons = rbind(medianRho_exons, medianRho_exons.control)
  # 
  # 
  # medianRho_introns = aggregate(wmean.rho ~ idx, data = df_pos_introns, median)
  # medianRho_introns$n.intervals = aggregate(wmean.rho ~ idx, data = df_pos_introns, length)$wmean.rho
  # medianRho_introns$condition = "observed"
  # medianRho_introns$hotoverlap = "both"
  # # medianRho_introns.control = aggregate(mean.rho.control ~ idx_resample, data = df_pos_introns, median)
  # medianRho_introns.control = aggregate(mean.rho.control ~ idx, data = df_pos_introns, median)
  # medianRho_introns.control$n.intervals = aggregate(mean.rho.control ~ idx, data = df_pos_introns, length)$mean.rho.control
  # medianRho_introns.control$condition = "control"
  # medianRho_introns.control$hotoverlap = "both"
  # colnames(medianRho_introns) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  # colnames(medianRho_introns.control) = c("idx", "medianRho", "n.intervals", "condition", "hotoverlap")
  # medianRho_introns = rbind(medianRho_introns, medianRho_introns.control)
  # 
  # 
  # medianRho.TTS$level = "all"
  # medianRho_exons$level = "exon" 
  # medianRho_introns$level = "intron"
  # 
  # medianRho_exons = medianRho_exons[which(medianRho_exons$idx >= 0),]
  # medianRho_introns = medianRho_introns[which(medianRho_introns$idx >= 0),]
  # 
  # medianRho.TTS = rbind(medianRho.TTS, medianRho_exons, medianRho_introns)
  # 
  # medianRho.TTS$position = "TSS"
  # 
  # 
  medianRho.TSS$atg_tss_dist = NA
  medianRho.TTS$atg_tss_dist = NA

  medianRho = rbind(medianRho.ATG, medianRho.TSS, medianRho.TTS)
  
  # Dataset
  medianRho$set = s
  
  if (exists("df_gradient")) {
    df_gradient = rbind(df_gradient, medianRho)
    rm(medianRho)
  } else {
    df_gradient = medianRho
    rm(medianRho)
  }
}
# Save the file ----
saveRDS(df_gradient, file = paste("Data/Recombination/Gradient/medianRho.rds", sep = ""))
rm(df_gradient)





