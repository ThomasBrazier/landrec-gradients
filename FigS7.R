# Loading env
source("Source/init_knitr.R")
source("Source/set_ggplot_theme.R")



# Rho ~ number of exons
df_subset = data_all[which(data_all$feature == "gene" & data_all$nb_exons <= max.exons),]

df_subset$species = factor(df_subset$species,
                           levels = list_species_ordered)

# meanRho = aggregate(weighted.mean.rho ~ nb_exons + species, data = df_subset, mean)
# meanRho$mean = NA
# meanRho$mean_lower_95 = NA
# meanRho$mean_upper_95 = NA
# n_boot = 1000

# for (i in 1:nrow(meanRho)) {
#   cat(i/nrow(meanRho)*100, "%\n")
#   # Mean
#   boot = numeric(n_boot)
#   boot = unlist(pbmclapply(1:n_boot, function(x) {mean(sample(df_subset$weighted.mean.rho[which(df_subset$species == meanRho$species[i] & df_subset$nb_exons == meanRho$nb_exons[i])], replace = TRUE), na.rm = TRUE)}, mc.cores = 1))
#   meanRho$mean[i] = mean(boot, na.rm = TRUE)
#   meanRho$mean_lower_95[i] = quantile(boot, 0.025, na.rm = TRUE)
#   meanRho$mean_upper_95[i] = quantile(boot, 0.975, na.rm = TRUE)
# }

### SAVE BOOTSTRAP
# saveRDS(meanRho, "Output/meanRho_genes_12species.rds")

fontsize = 16
dotsize = 0.2
linesize = 1


meanRho = readRDS("Output/meanRho_genes_12species.rds")


p1 = ggplot(data = meanRho, aes(x = nb_exons, y = mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(x = nb_exons, ymin = mean_lower_95, ymax = mean_upper_95)) +
  xlab("Number of exons") + ylab("Mean ρ/kb") +
  scale_x_continuous(breaks = c(1, 5, 10, 14)) +
  facet_wrap(~ species, ncol = 4, scale = "free") +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        axis.text.x=element_text(),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize, face = "italic"),
        legend.title=element_text(size=fontsize),
        legend.position='right')
p1




# medianRho = aggregate(weighted.mean.rho ~ nb_exons + species, data = df_subset, median)
# medianRho$median = NA
# medianRho$median_lower_95 = NA
# medianRho$median_upper_95 = NA
# n_boot = 1000

# for (i in 1:nrow(meanRho)) {
#   cat(i/nrow(meanRho)*100, "%\n")
#   # Mean
#   boot = numeric(n_boot)
#   boot = unlist(pbmclapply(1:n_boot, function(x) {median(sample(df_subset$weighted.mean.rho[which(df_subset$species == medianRho$species[i] & df_subset$nb_exons == medianRho$nb_exons[i])], replace = TRUE), na.rm = TRUE)}, mc.cores = 1))
#   medianRho$median[i] = mean(boot, na.rm = TRUE)
#   medianRho$median_lower_95[i] = quantile(boot, 0.025, na.rm = TRUE)
#   medianRho$median_upper_95[i] = quantile(boot, 0.975, na.rm = TRUE)
# }

### SAVE BOOTSTRAP
# saveRDS(medianRho, "Output/medianRho_genes_12species.rds")

fontsize = 16
dotsize = 0.2
linesize = 1


medianRho = readRDS("Output/medianRho_genes_12species.rds")


p2 = ggplot(data = medianRho, aes(x = nb_exons, y = median)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(x = nb_exons, ymin = median_lower_95, ymax = median_upper_95)) +
  xlab("Number of exons") + ylab("Median ρ/kb") +
  scale_x_continuous(breaks = c(1, 5, 10, 14)) +
  facet_wrap(~ species, ncol = 4, scale = "free") +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        axis.text.x=element_text(),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize, face = "italic"),
        legend.title=element_text(size=fontsize),
        legend.position='right')
p2

p = ggarrange(p1, p2, nrow = 2, labels = "AUTO")

ggsave(file = paste("Figure/Paper/FigS7.jpeg", sep = ""), p, width = 12, height = 18, dpi = 300)
