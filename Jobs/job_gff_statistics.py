#!/usr/bin/python
"""
Estimate GC content and genomic statistics with PiSlice
Dataset must be passed in arguments to the command line script
"""
from datetime import datetime

print("===============================")
print("Parse a GFF with PiSlice")
print(datetime.now())
print("===============================")

import os
import PiSlice.input as input
import PiSlice.core as core
import pandas as pd
import numpy as np
import sys
import allel
from itertools import cycle
import multiprocessing

try:
    dataset = sys.argv[1]
    print(dataset)
except:
    print("job_gff_statistics.py <dataset_name>")
    sys.exit(2)

# ncpus = multiprocessing.cpu_count()
ncpus = 16

marey_chromosome_metadata = pd.read_csv("Data/Recombination/marey_chromosome_metadata.csv", sep="\t")
marey_maps = pd.read_csv("Data/Recombination/Marey/marey_maps.csv", sep="\t")
marey_dataset = pd.read_csv("Data/Recombination/marey_dataset.csv", sep="\t")
chromosome_metadata = pd.read_csv("Data/Genome/genome_chromosome_metadata.csv", sep="\t")
polymorphism_metadata = pd.read_csv("Data/Polymorphism/polymorphism_metadata.csv", sep="\t")

if any(marey_dataset["set"] == dataset):
    species = (marey_dataset.loc[marey_dataset["set"] == dataset]["species"]).to_string(index=False, header=False).replace(' ', '')
    accession = (marey_dataset.loc[marey_dataset["set"] == dataset]["accession"]).to_string(index=False).replace(' ', '')
elif any(polymorphism_metadata["set"] == dataset):
    species = (polymorphism_metadata.loc[polymorphism_metadata["set"] == dataset]["species"]).to_string(index=False, header=False).replace(' ', '')
    accession = (polymorphism_metadata.loc[polymorphism_metadata["set"] == dataset]["accession"]).to_string(index=False).replace(' ', '')
else:
    print("No dataset")


print(species + " " + accession)

gff_file = "Data/Genomic_landscapes/GFF_parsed/" + dataset + ".csv.gz"
gff_parsed = input.read_gff(gff_file)


"""
Estimate GC content and other genomic statistics (SNP density, Pi) for each gene and CDS in gff
"""
fasta_file = "Data/Genome/" + species + "_" + accession + "/" + species + "_" + accession + ".fna.gz"
print(fasta_file)
genome = input.fasta(fasta_file) # Fasta must be bgzipped
chrnames = genome.seqname()
chrnames = set(chrnames)
print("Fasta chromosome names")
print(chrnames)

print("Translate FASTA chromosome names")
# [f(x) if condition else g(x) for x in sequence]
chr_metadata = chromosome_metadata[chromosome_metadata["set"] == dataset]
old_key = list(genome.seq.keys())
new_key = [chr_metadata[chr_metadata["annotname"] == o]["ldmapname"].item() if o in chr_metadata["annotname"].unique() else o for o in old_key]
for ok, nk in zip(old_key, new_key):
    print(f'{ok} -> {nk}')
    genome.seq[nk] = genome.seq.pop(ok)

chrnames = genome.seqname()
chrnames = set(chrnames)
print("NEW Fasta chromosome names")
print(chrnames)

chrnames = gff_parsed["seqname"]
chrnames = set(chrnames)
print("GFF chromosome names")
print(chrnames)

vcf_file = "Data/Polymorphism/" + dataset + "/" + dataset + ".pop.vcf.gz"
print(vcf_file)
vcf_file = allel.read_vcf(vcf_file)
chrnames = vcf_file["variants/CHROM"]
chrnames = set(chrnames)
print("VCF chromosome names")
print(chrnames)

#results = core.piSlice(windows=gff_parsed, statistics=["gc", "gc_noncoding", "gc3exon1", "gc_codon", "gc_intron", "cpg",
#                                                          "gene_length", "intron_length", "exon_length", "gene_nbexons",
#                                                          "missing_nucleotide", "gc_count", "at_count", "snp_count",
#                                                          "snp_count_at", "snp_count_gc"],
#                       fasta=genome, gff=gff_parsed, vcf=vcf_file, n_cpus=ncpus)

# results = core.piSlice(windows=gff_parsed, statistics=["gc", "gc_codon", "gc_intron",
#                                                          "exon_length", "intron_length", "gene_nbexons",
#                                                          "missing_nucleotide", "gc_count", "at_count", "snp_count",
#                                                          "snp_count_at", "snp_count_gc"],
#                       fasta=genome, gff=gff_parsed, vcf=vcf_file, n_cpus=ncpus)

results = core.piSlice(windows=gff_parsed, statistics=["missing_nucleotide", "gc_count", "at_count", "snp_count",
                                                        "snp_count_at", "snp_count_gc"],
                    fasta=genome, gff=gff_parsed, vcf=vcf_file, n_cpus=ncpus)



# Save results
filename = "Data/Genomic_landscapes/GFF_parsed/" + dataset + ".csv.gz"
input.write_gff2csv(results, filename)
