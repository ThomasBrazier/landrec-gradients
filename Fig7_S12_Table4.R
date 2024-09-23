# Loading env
source("Source/init_knitr.R")
source("Source/set_ggplot_theme.R")

library(modelr)

# # The gradient to fit with the toy model -----
# # Median of rho, median position of exon and median length of exon
# exons = data_all[which(data_all$feature == "CDS" & data_all$nb_exons <= max.exons),]
# # exons$dist_atg[which(exons$dist_atg < 0)] = NA
# # exons$dist_atg[which(exons$species == "Citrullus lanatus")] = abs(exons$dist_atg[which(exons$species == "Citrullus lanatus")])
# exons$dist_atg = abs(exons$dist_atg)

# # BUGFIX. Bug with G. max 14 exons and C. lanatus
# # Re-calculate distance to ATG
# cat("Relative position from ATG\n")
# gff_ranges = makeGRangesFromDataFrame(data_all[which(data_all$species == "Citrullus lanatus" & data_all$feature %in% c("gene", "mRNA", "exon", "CDS", "intron")),], keep.extra.columns = T)
# gff_ranges



# dist_atg = function(x) {
#   # Compute the distance to the ATG position
#   # Take an index of a row in a gff_ranges object in argument
#   atg = start(gff_ranges)[which(gff_ranges$gene_id == gff_ranges$gene_id[x] &
#                                   gff_ranges$feature == "CDS" & gff_ranges$rank == 1)]
#   atg
#   if (length(atg) > 0) {
#     # atg = max(atg) # Take the longest distance (longest transcript)
#     st = start(gff_ranges)[x]
#     atg
#     st
#     if (as.character(strand(gff_ranges)[x]) == "+") {
#       d_atg = st - atg
#       d_atg = min(abs(d_atg))
#       d_atg
#     }
#     if (as.character(strand(gff_ranges)[x]) == "-") {
#       d_atg = atg - st
#       d_atg = min(abs(d_atg))
#       d_atg
#     }
#   } else {
#     d_atg = NA
#   }
#   return(d_atg)
# }
# x = 1
# dist_atg(x)

# # Test
# atg = unlist(pbmclapply(1:1000, function(x) {dist_atg(x)}, mc.cores = 1))
# atg
# View(as.data.frame(gff_ranges[1:1000]))
# # All intervals
# atg = pbmclapply(1:length(gff_ranges), function(x) {dist_atg(x)}, mc.cores = 1)
# dist_atg = unlist(atg)
# dist_atg
# data_all$dist_atg[which(data_all$species == "Citrullus lanatus" & data_all$feature %in% c("gene", "mRNA", "exon", "CDS", "intron"))] = dist_atg



# genes = data_all[which(data_all$feature == "mRNA" & data_all$nb_exons <= max.exons),]

# colnames(genes)

# df = aggregate(weighted.mean.rho ~ rank + nb_exons + species, exons, median)
# df$nb_exons = as.factor(df$nb_exons)
# width = aggregate(width ~ rank + nb_exons + species, exons, mean)
# pos = aggregate(dist_atg ~ rank + nb_exons + species, exons, mean)
# df$exon.length = width$width
# df$exon.start.position = pos$dist_atg
# df$exon.centre.position = df$exon.start.position + df$exon.length/2

# ggplot(df, aes(x = exon.centre.position, y = weighted.mean.rho, group = nb_exons, color = nb_exons)) +
#   geom_point() +
#   geom_line() +
#   xlab("Exon centre position (bp)") + ylab("Median rho/kb") +
#   facet_wrap(~ species, scales = "free", ncol = 4) +
#   scale_color_viridis_d() +
#   theme_bw()



# View(exons[which(exons$species == "Citrullus lanatus"),])
# View(df[which(df$species == "Citrullus lanatus"),])


# # Compute CDS total length
# df$cds.length = NA

# for (i in 1:nrow(df)) {
#   idx = which(exons$species == df$species[i] & exons$nb_exons == df$nb_exons[i])
#   # idx = which(exons$species == df$species[i] & exons$nb_exons == df$nb_exons[i] & exons$rank == df$rank[i])
#   subset = exons[idx,]
#   exon.start = aggregate(start ~ parent, subset, min)
#   exon.end = aggregate(end ~ parent, subset, max)
#   cds = exon.start
#   cds$end = exon.end$end
#   cds$width = cds$end - cds$start
#   df$cds.length[i] = mean(cds$width, na.rm = TRUE)
# }



# # ggplot(df, aes(x = exon.centre.position, y = cds.length, group = nb_exons, color = nb_exons)) +
# #   geom_point() +
# #   geom_line() +
# #   xlab("Exon centre position (bp)") + ylab("Median rho/kb") +
# #   facet_wrap(~ species, scales = "free", ncol = 4) +
# #   scale_color_viridis_d() +
# #   theme_bw()
  
# write.table(df, file = "Output/gradients_fit.csv", sep = "\t",
# quote = F, col.names = T, row.names = F)



# Script to fit a gene gradient recombination model to rho estimates along genes
# Sylvain Glémin Februrary 2024 - sylvain.glemin@univ-rennes.fr


# library(ggplot2)
# library(gridExtra)


################# #
#### FUNCTIONS ####
################ #

# Expected recombination rate in an exon
# Ltot: total gene length
# Lexon: exon length
# Xstart: position of the begining of the exon
# par: vector with the parameters of the model
# rpred <- function(Ltot,Lexon,Xstart,par) {
#   r0 <- par[1]
#   p <- par[2]
#   a <- par[3]
#   b <- par[4]
#   x0 <- Xstart
#   x1 <- Xstart + Lexon
#   N1 <- (exp(-a*x0)-exp(-a*x1))*p/a
#   N2 <- ( exp(-b*(Ltot-x1)) -exp(-b*(Ltot-x0)))*(1-p)/b
#   return(r0*(N1+N2)/Lexon)
# }
rpred <- function(Ltot,Lexon,Xstart,par) {
  r0 <- par[1]
  r5 <- par[2]
  r3 <- par[3]
  a <- par[4]
  x0 <- Xstart
  x1 <- Xstart + Lexon
  N1 <- r5*(exp(-a*x0)-exp(-a*x1))/a
  N2 <- r3*(exp(-a*(Ltot-x1))-exp(-a*(Ltot-x0)))/a
  return(r0+(N1+N2)/Lexon)
}


# Sum of square difference between expected and observed recombination rates
# Function to be minimized
# Ltot: total gene length
# Lexon: exon length
# Xstart: position of the begining of the exon
# Robs: observed recombination rate
# par: vector with the parameters of the model
sumofsquare <- function(Ltot,Lexon,Xstart,Robs,par) {
  sum((rpred(Ltot,Lexon,Xstart,par)-Robs)^2)
}


# Function to fit the model by minimizinf least-square "by end"
# species: dataset for one species

# fitmodel <- function(species,ZERO=10^(-10),MAX=10^2) {
#   # data
#   Ltot <- species$cds.length
#   Lexon <- species$exon.length
#   Xstart <- species$exon.start.position
#   Robs <- species$weighted.mean.rho
#   # control parameters
#   inf <- c(-MAX,ZERO,ZERO,ZERO)
#   sup <- c(MAX,MAX,MAX,1)
#   # re-writing of the function
#   fn <- function(par) {
#     sumofsquare(Ltot,Lexon,Xstart,Robs,par)
#   }
#   #optimization (3 starting points)
#   r5init <- species[species$nb_exons==10 & species$rank==1,]$weighted.mean.rho
#   r3init <- species[species$nb_exons==10 & species$rank==10,]$weighted.mean.rho
#   init <- c(min(Robs),r5init,r3init,1/Ltot[1])
#   scaling <- init
#   mss <- optim(init,fn,lower=inf,upper=sup,method="L-BFGS-B",control=list(parscale=scaling,maxit=1000,factr=10^7,lmm=20,trace=1))
#   # init2 <- c(2*max(Robs),0.8,1/Ltot[1],1/Ltot[1])
#   # mss <- mss1
#   # scaling2 <- init2
#   # mss2 <- optim(init2,fn,lower=inf,upper=sup,method="L-BFGS-B",control=list(parscale=scaling2,maxit=1000,factr=10^7,lmm=20,trace=1))
#   # if(mss$value>mss2$value) mss <- mss2
#   # init3 <- c(2*max(Robs),0.2,1/Ltot[1],1/Ltot[1])
#   # scaling3 <- init3
#   # mss3 <- optim(init3,fn,lower=inf,upper=sup,method="L-BFGS-B",control=list(parscale=scaling3,maxit=1000,factr=10^7,lmm=20,trace=1))
#   # if(mss$value>mss3$value) mss <- mss3
#   # Goodness of fit
#   SSfit <- mss$value
#   SStot <- sum((Robs-mean(Robs))^2)
#   r2 <- 1 - SSfit/SStot
#   # result to return
#   result <- list("param"=mss$par,"r2"=r2,"pred"=rpred(species$cds.length,species$exon.length,species$exon.start.position,mss$par))
#   return(result)
# }

# Function to fit the model using the nls function
# species: dataset for one species

nlsfit <- function(species,ZERO=10^(-8),MAX=10^2) {
  # data
  Ltot <- species$cds.length
  Lexon <- species$exon.length
  Xstart <- species$exon.start.position
  Robs <- species$weighted.mean.rho
  df <- data.frame("Robs"=Robs,"Ltot"=Ltot,"Lexon"=Lexon,"Xtart"=Xstart)
  # loop to insure convergence: incrementing the ZERO value if error
  repeat{
    # Control parameters
    inf <- c(-MAX,ZERO,ZERO,ZERO)
    sup <- c(MAX,MAX,MAX,1)
    DONE <- F
    tryCatch({
      # Starting points for optimization
      r5init <- species[species$nb_exons==10 & species$rank==1,]$weighted.mean.rho
      r3init <- species[species$nb_exons==10 & species$rank==10,]$weighted.mean.rho
      par_init <- c(min(Robs),r5init,r3init,1/Ltot[1])
      mss <- nls(Robs ~ rpred(Ltot,Lexon,Xstart,par),data = df,start = list(par=par_init),algorithm = "port",lower=inf,upper=sup,
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024, printEval = TRUE, warnOnly = FALSE))
      DONE <- T
    },error=function(e){})
    if(DONE) break
    ZERO <- ZERO*2
  }
  # result to return
  r2 <- rsquare(mss,df)
  fitted_par <- summary(mss)$par[c(1:4),1]
  se_par <- summary(mss)$par[c(1:4),2]
  result <- list("param"=fitted_par,"se"=se_par,"r2"=r2,"pred"=rpred(species$cds.length,species$exon.length,species$exon.start.position,fitted_par))
  return(result)
}


# Function to plot the result
# species: dataset for one species
plotresult <- function(species,result) {
  g1 <- ggplot(species, aes(x = exon.centre.position, y = weighted.mean.rho, group = as.factor(nb_exons), color = as.factor(nb_exons))) +
    geom_point() +
    geom_line() +
    xlab("Exon centre position (bp)") + ylab("Median rho/kb") +
    facet_wrap(~ sp_name, scales = "free") +
    scale_color_viridis_d() +
    theme_bw()
  g2 <- ggplot(species, aes(x = exon.centre.position, y = rp, group = as.factor(nb_exons), color = as.factor(nb_exons))) +
    geom_point() +
    geom_line() +
    xlab("Exon centre position (bp)") + ylab("Median rho/kb") +
    facet_wrap(~ sp_name, scales = "free") +
    scale_color_viridis_d() +
    theme_bw()
  g3 <- ggplot(data = species,aes(x = rp,y=weighted.mean.rho)) +
    geom_point() +
    geom_abline(slope = 1,intercept = 0) +
    xlab("Predicted rho/kb") + ylab("Observed rho/kb") +
    ggtitle(paste("R2 = ",round(result$r2,2))) +
    theme_bw()
  tab <- data.frame("param"=c("r0","r5","r3","a"),"estimate"=result$param,"se"=result$se)
  t1 <- ggplot() +
    geom_text(data = tab,aes(x=1,y=param, label=paste(round(estimate,5),"   (",round(se,5),")",sep="")),size=5) +
    ggtitle("Parameters estimation") +
    theme_classic(base_size=10) +
    #Remove axis line, ticks and axis label to get the look of a table
    theme(
      axis.line = element_blank(),
      axis.ticks  = element_blank(),
      axis.title.y  = element_blank(),
      axis.title.x  = element_blank(),
      axis.text.x = element_text(color="white"),
      text = element_text(size=16)
    )
  pdf(file = paste("figures/",species$genus[1],"_",species$species[1],".pdf",sep=""),width = 8,height = 6)
    grid.arrange(g1,g2,g3,t1,ncol=2,heights=c(2,1))
  dev.off()
}


############ #
#### MAIN ####
############ #

# mydata <- read.table("Output/gradients_fit.csv", header = T, sep = "\t")
# mydata_corrected <- read.table("Data/FitGradient/gradients_fit_v1.csv", header = T, sep = " ")

# # mydata$sp_name <- paste(mydata$genus, mydata$species, sep=" ")
# colnames(mydata)
# colnames(mydata_corrected)
# which(mydata$species == "Citrullus lanatus")
# which(mydata_corrected$species == "lanatus")

# which(mydata$species == "Glycine max")
# which(mydata_corrected$species == "max")

# mydata$exon.start.position[which(mydata$species == "Citrullus lanatus")] = mydata_corrected$exon.start.position[which(mydata_corrected$species == "lanatus")]
# mydata$exon.centre.position[which(mydata$species == "Citrullus lanatus")] = mydata_corrected$exon.centre.position[which(mydata_corrected$species == "lanatus")]
# mydata$exon.length[which(mydata$species == "Citrullus lanatus")] = mydata_corrected$exon.length[which(mydata_corrected$species == "lanatus")]
# mydata$cds.length[which(mydata$species == "Citrullus lanatus")] = mydata_corrected$cds.length[which(mydata_corrected$species == "lanatus")]

# mydata$exon.start.position[which(mydata$species == "Glycine max")] = mydata_corrected$exon.start.position[which(mydata_corrected$species == "max")]
# mydata$exon.centre.position[which(mydata$species == "Glycine max")] = mydata_corrected$exon.centre.position[which(mydata_corrected$species == "max")]
# mydata$exon.length[which(mydata$species == "Glycine max")] = mydata_corrected$exon.length[which(mydata_corrected$species == "max")]
# mydata$cds.length[which(mydata$species == "Glycine max")] = mydata_corrected$cds.length[which(mydata_corrected$species == "max")]



# # Save corrected data
# write.table(mydata, "Output/gradients_fit_corrected.csv", col.names = T, row.names = F, sep = "\t", quote = F)

mydata <- read.table("Output/gradients_fit_corrected.csv", header = T, sep = "\t")


mydata$nb_exons = factor(mydata$nb_exons)
ggplot(mydata, aes(x = exon.centre.position, y = weighted.mean.rho, group = nb_exons, color = nb_exons)) +
  geom_point() +
  geom_line() +
  xlab("Exon centre position (bp)") + ylab("Median rho/kb") +
  facet_wrap(~ species, scales = "free", ncol = 4) +
  scale_color_viridis_d() +
  theme_bw()


# Not working for Citrullus
# mydata = mydata[which(mydata$species != "Citrullus lanatus"),]
# sp = "Citrullus lanatus"
# table(species$nb_exons)

rm(model_fit)
for(sp in levels(as.factor(mydata$species))) {
  species <- mydata[mydata$species == sp,]
  result <- nlsfit(species)
  result$param
  species$r0 = result$param[1]
  species$r5 = result$param[2]
  species$r3 = result$param[3]
  species$a = result$param[4]
  species$rp <- result$pred
  species$r2 = result$r2
  species$r0_se = result$se[1]
  species$r5_se = result$se[2]
  species$r3_se = result$se[3]
  species$a_se = result$se[4]
  # plotresult(species, result)
  if (exists("model_fit")) {
    model_fit = rbind(model_fit, species)
  } else {
    model_fit = species
  }
}


View(model_fit)



### Figure 7 - Gradient Fitting
fontsize = 22
dotsize = 0.2
linesize = 1


# toymodel_file = system.file("extdata", "Figure/ConceptModel.jpg", package = "cowplot")
p1 = ggplot() +
  ggtitle("Conceptual model") +
  draw_image("Figure/ConceptModel_v2.jpeg") +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        axis.text.x=element_text(),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize - 2, face = "italic"),
        legend.title=element_text(size=fontsize - 2),
        legend.position='right')
p1

p2 = ggplot(data = model_fit[which(model_fit$species %in% list_species),], aes(x = exon.centre.position/1000, y = rp, group = as.factor(nb_exons), color = as.factor(nb_exons))) +
  geom_point() +
  geom_line() +
  scale_color_viridis_d() +
  xlab("Exon centre position (kb)") + ylab("Median ρ/kb") +
  labs(colour = "# exons   ") +
  ggtitle("\nFitted model") +
  # geom_line(data = sim_gradient[which(sim_gradient$nb_exons == "all"),], aes(x = position/100, y = rho), color = "black") +
  facet_wrap(. ~ species, scales = "free", ncol = 3) +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        axis.text.x=element_text(),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(1,"line"),
        legend.key.width = unit(2,"line"),
        legend.text=element_text(size=fontsize - 6, face = "italic"),
        legend.title=element_text(size=fontsize - 4),
        legend.position='right')
p2


p2b = ggplot(data = model_fit[which(model_fit$species %in% list_species),], aes(x = rp,y=weighted.mean.rho)) +
  geom_point() +
  geom_abline(slope = 1,intercept = 0) +
  xlab("Predicted ρ/kb") + ylab("Observed ρ/kb") +
  facet_wrap(. ~ species, scales = "free", ncol = 3) +
  ggtitle("\nGoodness of fit") +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        axis.text.x=element_text(),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize - 2, face = "italic"),
        legend.title=element_text(size=fontsize - 2),
        legend.position='right')
p2b


ggpubr::ggarrange(p1, p2, p2b, nrow = 3, heights = c(2,3,2.5), labels = "AUTO", font.label = list(size = 24)) + bgcolor("White") + border(color = "White")

ggsave(file = paste("Figure/Paper/Fig7.jpeg", sep = ""), width = 15, height = 12, dpi = 300)

ggsave(file = paste("Figure/Paper/Fig7.jpeg", sep = ""), width = 15, height = 12, dpi = 600)


# Supplementary Figure
p3 = ggplot(data = model_fit, aes(x = exon.centre.position/1000, y = rp, group = as.factor(nb_exons), color = as.factor(nb_exons))) +
  geom_point() +
  geom_line() +
  scale_color_viridis_d() +
  xlab("Exon centre position (kb)") + ylab("Median ρ/kb") +
  labs(colour = "# exons   ") +
  # geom_line(data = sim_gradient[which(sim_gradient$nb_exons == "all"),], aes(x = position/100, y = rho), color = "black") +
  facet_wrap(. ~ species, scales = "free", ncol = 3) +
  theme(plot.title = element_text(color="black", size=fontsize, face="bold",hjust = 0.5),
        axis.title.x = element_text(color="black", size=fontsize),
        axis.title.y = element_text(color="black", size=fontsize),
        axis.text=element_text(size=fontsize, colour="black"),
        axis.text.x=element_text(),
        strip.text.x = element_text(size = fontsize, face = "italic"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(3,"line"),
        legend.text=element_text(size=fontsize-2),
        legend.title=element_text(size=fontsize-2),
        legend.position='right')
p3

ggsave(file = paste("Figure/Paper/FigS12.jpeg", sep = ""), p3, width = 15, height = 10, dpi = 300)


# Save TABLE ----
head(model_fit)
res = aggregate(r0 ~ species, model_fit, unique)
res$r0_se = aggregate(r0_se ~ species, model_fit, unique)$r0
res$r5 = aggregate(r5 ~ species, model_fit, unique)$r5
res$r5_se = aggregate(r5_se ~ species, model_fit, unique)$r5_se
res$r3 = aggregate(r3 ~ species, model_fit, unique)$r3
res$r3_se = aggregate(r3_se ~ species, model_fit, unique)$r3_se
res$ratio = res$r5/res$r3
res$a = aggregate(a ~ species, model_fit, unique)$a
res$a_se = aggregate(a_se ~ species, model_fit, unique)$a_se
res$tract = 1/res$a
res$r_squared = aggregate(r2 ~ species, model_fit, unique)$r2


res$r0 = round(res$r0, 2)
res$r0_se = round(res$r0_se, 2)
res$r5 = round(res$r5, 2)
res$r5_se = round(res$r5_se, 2)
res$r3 = round(res$r3, 2)
res$r3_se = round(res$r3_se, 2)
res$ratio = round(res$ratio, 2)
res$a = round(res$a, 6)
res$a_se = round(res$a_se, 6)
res$tract = round(res$tract, 0)
res$r_squared = round(res$r_squared, 2)

# View(res)

write.xlsx(x = res, file = paste("Table/Table4.xls", sep = ""), sheetName = "Model fit", row.names = FALSE, append = FALSE)
