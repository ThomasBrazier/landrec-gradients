require(vcfR)
require(adegenet)
require(poppr)
require(ggplot2)
require(yaml)
require(GenomicRanges)

source("Source/init.R")


for (dataset in list_dataset) {
  vcf_file = paste("Data/Polymorphism/", dataset, "/", dataset, ".pop.vcf.gz", sep = "")
  cat("Processing dataset", dataset, "\n")
  # Read the vcf
  vcf = read.vcfR(vcf_file, verbose = FALSE)
  
  # Extract SNP coordinates and REF/ALT alleles
  head(vcf@gt)
  
  chromosome = vcf@fix[,1]
  start = as.integer(vcf@fix[,2])
  end = (as.integer(vcf@fix[,2]) + 1)
  ref = vcf@fix[,4]
  alt = vcf@fix[,5]

  df = data.frame(chromosome = chromosome,
                  start = start,
                  end = end,
                  ref = ref,
                  alt = alt)
  
  snpranges = makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
  
  # Load GFF parsed and keep a subset of columns
  gff = readRDS(file = paste("Data/Genomic_landscapes/Rho/gff_rho_", dataset, ".rds", sep = ""))
  colnames(gff)
  
  gff = gff[,c("chromosome", "start", "end", "width", "strand", "feature",
               "id", "parent", "gene_biotype", "gene_id", "rank",
               "nb_exons")]
  
  gffranges = makeGRangesFromDataFrame(gff, keep.extra.columns = TRUE)
  
  
  # Count the number of SNPs
  gffranges$n_snp = countOverlaps(gffranges, snpranges, type = "any")
  snpranges_GC = snpranges[which(snpranges$ref %in% c("G", "C") & snpranges$alt %in% c("G", "C"))]
  gffranges$n_snp_GC = countOverlaps(gffranges, snpranges_GC, type = "any")
  snpranges_AT = snpranges[which(snpranges$ref %in% c("A", "T") & snpranges$alt %in% c("A", "T"))]
  gffranges$n_snp_AT = countOverlaps(gffranges, snpranges_AT, type = "any")
  
  # SAVE
  gffranges = as.data.frame(gffranges)
  
  saveRDS(gffranges, file = paste("Data/Genomic_landscapes/Polymorphism/gff_polymorphism_", dataset, ".rds", sep = ""))
  
}