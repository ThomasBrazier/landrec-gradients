#!/usr/bin/env Rscript

source("Source/init.R")

#============================================================================#
# Mean Rho in genomic features ----
#============================================================================#
cat("======================================\n")
cat("Mean rho in genomic features\n")

# table(data_all$feature[which(data_all$set == "Arabidopsis_thaliana_1001genomes")])
# table(data_all$feature[which(data_all$set == "Malus_domestica_Sun2020")])

# mean(df_tmp$mean.rho[which(df_tmp$feature %in% c("gene"))], na.rm = TRUE)
# mean(df_tmp$mean.rho[which(df_tmp$feature %in% c("exon"))], na.rm = TRUE)
# mean(df_tmp$mean.rho[which(df_tmp$feature %in% c("intron"))], na.rm = TRUE)
# mean(df_tmp$mean.rho[which(df_tmp$feature %in% c("utr5"))], na.rm = TRUE)
# mean(df_tmp$mean.rho[which(df_tmp$feature %in% c("utr3"))], na.rm = TRUE)
# mean(df_tmp$mean.rho[which(df_tmp$feature %in% c("five_prime_UTR"))], na.rm = TRUE)
# mean(df_tmp$mean.rho[which(df_tmp$feature %in% c("three_prime_UTR"))], na.rm = TRUE)
# mean(df_tmp$mean.rho[which(df_tmp$feature %in% c("exon", "intron", "utr5", "utr3"))], na.rm = TRUE)

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


# Bootstrap
n_boot = 1000

for (i in 1:nrow(df_meanrho)) {
  cat(round(i/nrow(df_meanrho)*100), "%\n")
  # Mean
  boot = numeric(n_boot)
  boot = unlist(pbmclapply(1:n_boot, function(x) {mean(sample(data_all$mean.rho[which(data_all$species == df_meanrho$species[i] & data_all$feature == df_meanrho$feature[i])], replace = TRUE), na.rm = TRUE)}))
  df_meanrho$mean[i] = mean(boot, na.rm = TRUE)
  df_meanrho$mean_lower_95[i] = quantile(boot, 0.025, na.rm = TRUE)
  df_meanrho$mean_upper_95[i] = quantile(boot, 0.975, na.rm = TRUE)
  # Weighted mean
  boot = numeric(n_boot)
  boot = unlist(pbmclapply(1:n_boot, function(x) {mean(sample(data_all$weighted.mean.rho[which(data_all$species == df_meanrho$species[i] & data_all$feature == df_meanrho$feature[i])], replace = TRUE), na.rm = TRUE)}))
  df_meanrho$wmean[i] = mean(boot, na.rm = TRUE)
  df_meanrho$wmean_lower_95[i] = quantile(boot, 0.025, na.rm = TRUE)
  df_meanrho$wmean_upper_95[i] = quantile(boot, 0.975, na.rm = TRUE)
  # Median
  boot = numeric(n_boot)
  boot = unlist(pbmclapply(1:n_boot, function(x) {median(sample(data_all$weighted.mean.rho[which(data_all$species == df_meanrho$species[i] & data_all$feature == df_meanrho$feature[i])], replace = TRUE), na.rm = TRUE)}))
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
  
  df_meanrho$hotspot_countfiltered[i] = sum(data_all$hotspot_overlap_intensity4[which(data_all$species == df_meanrho$species[i] & data_all$feature == df_meanrho$feature[i])], na.rm = TRUE)
  boot = numeric(n_boot)
  boot = unlist(pbmclapply(1:n_boot, function(x) {sum(sample(data_all$hotspot_overlap_intensity4[which(data_all$species == df_meanrho$species[i] & data_all$feature == df_meanrho$feature[i])], replace = TRUE), na.rm = TRUE)}))
  df_meanrho$hotspot_countfiltered_boot[i] = mean(boot, na.rm = TRUE)
  df_meanrho$hotspot_countfiltered_lower_95[i] = quantile(boot, 0.025, na.rm = TRUE)
  df_meanrho$hotspot_countfiltered_upper_95[i] = quantile(boot, 0.975, na.rm = TRUE)
  # df_meanrho$hotspot_meannb[i] = df_meanrho$hotspot_count[i]/df_meanrho$n_features[i] # hotspot_count/nb of features
  # df_meanrho$hotspot_proportion[i] = df_meanrho$hotspot_count[i]/sum(data_all$hotspot_overlap[which(data_all$species == df_meanrho$species[i])], na.rm = TRUE) # hotspot_count/total number of hotspots
  
}

saveRDS(df_meanrho, "Data/Recombination/meanrho_per_feature.rds")





