#!/usr/bin/env Rscript
source("Source/init.R")

chromosome_metadata = chromosome_metadata[which(!is.na(chromosome_metadata$ldmapname)),]


cat("======================================\n")
cat("Mean rho weighted per sequence length in genomic features\n")

# Compute mean rho ± 95% CI for each feature
features =c("intergenic", "intergenic_buffer",
            "buffer_upstream", "buffer_downstream",
            "gene", "mRNA", "promoter",
            "upstream1kb", "downstream1kb",
            "upstream2kb", "downstream2kb",
            "upstream3kb", "downstream3kb",
            "utr3", "utr5",
            "CDS", "exon", "intron")
df_meanrho = data.frame(species = rep(unique(data_all$species), each = length(features)),
                        feature = rep(features, length(unique(data_all$species))))
# TODO exon 1, 2, 3 and last

df_meanrho$mean = NA
df_meanrho$mean_lower_95 = NA
df_meanrho$mean_upper_95 = NA
df_meanrho$wmean = NA
df_meanrho$wmean_lower_95 = NA
df_meanrho$wmean_upper_95 = NA
# Log rho to reduce high variance of estimators
# df_meanrho$wmean_log10 = NA
# df_meanrho$wmean_log10_lower_95 = NA
# df_meanrho$wmean_log10_upper_95 = NA
df_meanrho$median = NA
df_meanrho$median_lower_95 = NA
df_meanrho$median_upper_95 = NA

# Sampling size
df_meanrho$n_features = NA
df_meanrho$n_bp = NA

# Hotspots
df_meanrho$hotspot_count = NA
df_meanrho$hotspot_countfiltered = NA
df_meanrho$hotspot_count_boot = NA
df_meanrho$hotspot_count_lower_95 = NA
df_meanrho$hotspot_count_upper_95 = NA
df_meanrho$hotspot_countfiltered_boot = NA
df_meanrho$hotspot_countfiltered_lower_95 = NA
df_meanrho$hotspot_countfiltered_upper_95 = NA

# Data and estimators quality
df_meanrho$variance = NA
df_meanrho$snp_density_kb = NA
df_meanrho$mean_size_bp = NA


# Bootstrap
n_boot = 1000
for (i in 1:nrow(df_meanrho)) {
  cat(round(i/nrow(df_meanrho)*100), "%\n")
  # Mean
  wmean = function(x) {
    idx = sample(which(data_all$species == df_meanrho$species[i] & data_all$feature == df_meanrho$feature[i]), replace = TRUE)
    m = data_all$mean.rho[idx]
    w = data_all$width[idx]
    weighted.mean(m, w, na.rm = TRUE)
  }
  boot = numeric(n_boot)
  boot = unlist(pbmclapply(1:n_boot, wmean))
  df_meanrho$mean[i] = mean(boot, na.rm = TRUE)
  df_meanrho$mean_lower_95[i] = quantile(boot, 0.025, na.rm = TRUE)
  df_meanrho$mean_upper_95[i] = quantile(boot, 0.975, na.rm = TRUE)
  # Weighted mean
  wmean = function(x) {
    idx = sample(which(data_all$species == df_meanrho$species[i] & data_all$feature == df_meanrho$feature[i]), replace = TRUE)
    m = data_all$weighted.mean.rho[idx]
    w = data_all$width[idx]
    weighted.mean(m, w, na.rm = TRUE)
  }
  boot = numeric(n_boot)
  boot = unlist(pbmclapply(1:n_boot, wmean))
  df_meanrho$wmean[i] = mean(boot, na.rm = TRUE)
  df_meanrho$wmean_lower_95[i] = quantile(boot, 0.025, na.rm = TRUE)
  df_meanrho$wmean_upper_95[i] = quantile(boot, 0.975, na.rm = TRUE)
  # # Weighted mean with log10 transform
  # wmeanlog = function(x) {
  #   idx = sample(which(data_all$species == df_meanrho$species[i] & data_all$feature == df_meanrho$feature[i]), replace = TRUE)
  #   m = log10(data_all$weighted.mean.rho[idx])
  #   w = data_all$width[idx]
  #   weighted.mean(m, w, na.rm = TRUE)
  # }
  # boot = numeric(n_boot)
  # boot = unlist(pbmclapply(1:n_boot, wmeanlog))
  # df_meanrho$wmean_log10[i] = mean(boot, na.rm = TRUE)
  # df_meanrho$wmean_log10_lower_95[i] = quantile(boot, 0.025, na.rm = TRUE)
  # df_meanrho$wmean_log10_upper_95[i] = quantile(boot, 0.975, na.rm = TRUE)
  
  # Median
  med = function(x) {
    idx = sample(which(data_all$species == df_meanrho$species[i] & data_all$feature == df_meanrho$feature[i]), replace = TRUE)
    m = data_all$weighted.mean.rho[idx]
    w = data_all$width[idx]
    median(m, na.rm = TRUE)
  }
  boot = numeric(n_boot)
  boot = unlist(pbmclapply(1:n_boot, med))
  df_meanrho$median[i] = mean(boot, na.rm = TRUE)
  df_meanrho$median_lower_95[i] = quantile(boot, 0.025, na.rm = TRUE)
  df_meanrho$median_upper_95[i] = quantile(boot, 0.975, na.rm = TRUE)
  
  # Sampling size
  # df_meanrho$n_features[i] = sum(data_all$species == df_meanrho$species[i] & data_all$feature == df_meanrho$feature[i], na.rm = TRUE)
  # The number of bp each feature represents
  df_meanrho$n_bp[i] = sum(data_all$width[which(data_all$species == df_meanrho$species[i] & data_all$feature == df_meanrho$feature[i])], na.rm = TRUE)
  # The number of instance of each feature
  df_meanrho$n_features[i] = length(which(data_all$feature == df_meanrho$feature[i] & data_all$species == df_meanrho$species[i]))
  
  # Hotspots
  # df_meanrho$hotspot_count[i] = sum(data_all$species == df_meanrho$species[i] & data_all$feature == df_meanrho$feature[i] & data_all$hotspot_overlap > 0, na.rm = TRUE)
  # df_meanrho$hotspot_meannb[i] = df_meanrho$hotspot_count[i]/df_meanrho$n_features[i] # hotspot_count/nb of features
  # df_meanrho$hotspot_proportion[i] = df_meanrho$hotspot_count[i]/sum(data_all$species == df_meanrho$species[i] & data_all$hotspot_overlap > 0, na.rm = TRUE) # hotspot_count/total number of hotspots
  df_meanrho$hotspot_count[i] = sum(data_all$hotspot_overlap_raw[which(data_all$species == df_meanrho$species[i] & data_all$feature == df_meanrho$feature[i])], na.rm = TRUE)
  boot = numeric(n_boot)
  boot = unlist(pbmclapply(1:n_boot, function(x) {sum(sample(data_all$hotspot_overlap_raw[which(data_all$species == df_meanrho$species[i] & data_all$feature == df_meanrho$feature[i])], replace = TRUE), na.rm = TRUE)}))
  df_meanrho$hotspot_count_boot[i] = mean(boot, na.rm = TRUE)
  df_meanrho$hotspot_count_lower_95[i] = quantile(boot, 0.025, na.rm = TRUE)
  df_meanrho$hotspot_count_upper_95[i] = quantile(boot, 0.975, na.rm = TRUE)
  
  df_meanrho$hotspot_count[i] = sum(data_all$hotspot_overlap_intensity4[which(data_all$species == df_meanrho$species[i] & data_all$feature == df_meanrho$feature[i])], na.rm = TRUE)
  boot = numeric(n_boot)
  boot = unlist(pbmclapply(1:n_boot, function(x) {sum(sample(data_all$hotspot_overlap_intensity4[which(data_all$species == df_meanrho$species[i] & data_all$feature == df_meanrho$feature[i])], replace = TRUE), na.rm = TRUE)}))
  df_meanrho$hotspot_countfiltered_boot[i] = mean(boot, na.rm = TRUE)
  df_meanrho$hotspot_countfiltered_lower_95[i] = quantile(boot, 0.025, na.rm = TRUE)
  df_meanrho$hotspot_countfiltered_upper_95[i] = quantile(boot, 0.975, na.rm = TRUE)
  # df_meanrho$hotspot_meannb[i] = df_meanrho$hotspot_count[i]/df_meanrho$n_features[i] # hotspot_count/nb of features
  # df_meanrho$hotspot_proportion[i] = df_meanrho$hotspot_count[i]/sum(data_all$hotspot_overlap[which(data_all$species == df_meanrho$species[i])], na.rm = TRUE) # hotspot_count/total number of hotspots
  
  # Data and estimators quality
  df = data_all[which(data_all$species == df_meanrho$species[i] & data_all$feature == df_meanrho$feature[i]),]
  df_meanrho$variance[i] = var(df$weighted.mean.rho, na.rm = TRUE)
  df_meanrho$snp_density_kb[i] = mean(df$snp_count/(df$width * 10^3), na.rm = TRUE)
  df_meanrho$mean_size_bp[i] = mean(df$width, na.rm = TRUE)
}

saveRDS(df_meanrho, "Data/Recombination/weightedmeanrho_per_feature.rds")
