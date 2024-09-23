# Loading env
dataset_metadata = readODS::read_ods("Data/Metadata/dataset.ods")

### Table S1
tableS1 = dataset_metadata[c("dataset", "ref", "year", "species", "mating.system", "nb.SNP.raw", "nb.SNP.trim", "nb.individuals.total", "nb.individuals.sampled", "accession", "n.chromosome", "n.genes", "ref.genome", "public.database", "link_data", "doi")]
tableS1

write.xlsx(x = as.data.frame(tableS1), file = "Table/TableS1.xls", sheetName = "Metadata", row.names = FALSE, append = FALSE)
