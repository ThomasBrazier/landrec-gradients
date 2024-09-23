require(vcfR)
require(adegenet)
require(poppr)
require(ggplot2)
require(yaml)
require(GenomicRanges)

source("Source/init.R")

ncores = 16

args = commandArgs(trailingOnly=TRUE)
dataset = args[1]

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
             "nb_exons", "nb_snp", "missing_nucleotide", "mean.rho", "median.rho", "mean.rho.control",
             "median.rho.control", "weighted.mean.rho", "weighted.mean.rho.control",
             "mean.rho.startafter", "weighted.mean.rho.startafter", "relative.meanrho",
             "relative.wmeanrho", "dist_atg")]

gff = gff[which(gff$feature %in% c("gene", "mRNA", "CDS", "intron", "exon")),]

gffranges = makeGRangesFromDataFrame(gff, keep.extra.columns = TRUE)


# Count the number of SNPs
gffranges$n_snp = countOverlaps(gffranges, snpranges, type = "any")
snpranges_GC = snpranges[which(snpranges$ref %in% c("G", "C") & snpranges$alt %in% c("G", "C"))]
gffranges$n_snp_GC = countOverlaps(gffranges, snpranges_GC, type = "any")
snpranges_AT = snpranges[which(snpranges$ref %in% c("A", "T") & snpranges$alt %in% c("A", "T"))]
gffranges$n_snp_AT = countOverlaps(gffranges, snpranges_AT, type = "any")

# Mean Pi per site
sites_pi = read.table(paste0("Data/Polymorphism/", dataset, "/", dataset, ".sites.pi"),
                      header = T,
                      sep = "\t")
colnames(sites_pi) = c("chromosome", "start", "pi")
sites_pi$end = sites_pi$start + 1
sites_pi$strand = "*"
sites_pi = makeGRangesFromDataFrame(sites_pi, keep.extra.columns = T)

hits = findOverlaps(gffranges, sites_pi, type = "any")

idx = c(1:length(gffranges))

pi = unlist(pbmclapply(idx, function(x) {if (x %in% queryHits(hits)) {sum(sites_pi$pi[subjectHits(hits)[queryHits(hits) == x]], na.rm = TRUE)} else {NA}},
                       mc.cores = ncores))

gffranges$pi_bp = pi/width(gffranges)

# SAVE
gffranges = as.data.frame(gffranges)

saveRDS(gffranges, file = paste("Data/Genomic_landscapes/Polymorphism/gff_pi_", dataset, ".rds", sep = ""))
