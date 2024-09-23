rm(list=ls(all=TRUE))
source("Source/init.R")
gc()

summary(data_all)
colnames(data_all)

gene = data_all[which(data_all$feature == "gene" & data_all$nb_exons < 15),]


#============================================================================#
# Number of exons ----
#============================================================================#
sample_size = data.frame(species = levels(gene$species),
                         table(gene$species, gene$nb_exons)[,1],
                         table(gene$species, gene$nb_exons)[,2],
                         table(gene$species, gene$nb_exons)[,3],
                         table(gene$species, gene$nb_exons)[,4],
                         table(gene$species, gene$nb_exons)[,5],
                         table(gene$species, gene$nb_exons)[,6],
                         table(gene$species, gene$nb_exons)[,7],
                         table(gene$species, gene$nb_exons)[,8],
                         table(gene$species, gene$nb_exons)[,9],
                         table(gene$species, gene$nb_exons)[,10],
                         table(gene$species, gene$nb_exons)[,11],
                         table(gene$species, gene$nb_exons)[,12],
                         table(gene$species, gene$nb_exons)[,13],
                         table(gene$species, gene$nb_exons)[,14])
colnames(sample_size) = c("Species", as.character(1:14))
sample_size


write.xlsx(x = sample_size, file = "Table/TableS3.xls", sheetName = "Gene sample size", row.names = FALSE, append = FALSE)
