#!/usr/bin/python
"""
Estimate GC content with PiSlice
Set, species, accession and chromosome must be passed in arguments to the command line script
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
import sys
import allel
from itertools import cycle

try:
    dataset = sys.argv[1]
    print(dataset)
except:
    print("job_gff_parse.py <dataset_name>")
    sys.exit(2)



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

"""
Estimate GC/GC123 for each gene and CDS in gff, with ranks and parent/children relationships
"""
#marey = marey_maps[(marey_maps['set'] == set)].copy()
#gff = input.read_gff(gff_file, parse=True, parse_introns=True)
# Estimating GC123 requires a gff and a genome
prefix = species + "_" + accession + "/" + species + "_" + accession
gff_file = "Data/Genome/" + prefix + ".gff.gz"
gff = input.read_gff(gff_file, parse=True, infer_rank=False, parse_introns=False, parse_utr=False, n_cpus=16)
#gff = input.read_gff(gff_file, parse=False)
fasta_file = "Data/Genome/" + prefix + ".fna.gz"
genome = input.fasta(fasta_file) # Fasta must be bgzipped

# Subset only chromosomes
#list_chromosomes = list(marey_maps[(marey_maps['set'] == set)]["chromosome.genome"].unique())
list_chromosomes = list(chromosome_metadata[(chromosome_metadata['set'] == dataset)]["annotname"].unique())
gff = gff.loc[gff["seqname"].isin(list_chromosomes),]

#gff = gff.loc[gff["feature"].isin(["region", "gene", "CDS", "exon"]),]

# Keep only genes that are not ambiguous: protein_coding or "" (assume protein coding)
# See https://m.ensembl.org/info/genome/genebuild/biotypes.html
# ambiguous_genes = gff.loc[(gff["feature"] == "gene") & ~((gff["gene_biotype"] == "protein_coding") | (gff["gene_biotype"] == ""))]
# ambiguous_genes_ids = list(ambiguous_genes["id"])
# ambiguous_genes_ids = [x for x in ambiguous_genes_ids if x != '']
# # Get ids of children
# attributes = list(gff["attribute"])
# attr = re.match(str(x), attributes)
# indices = [i for i, y in enumerate(attributes) if re.match(x, y)]
# children_ambiguous_genes = [gff.gff.children(x, all=True) for x in list(ambiguous_genes_ids)]
# children_ambiguous_genes = gff.gff.children(list(ambiguous_genes_ids), all=True)
# idx = list(ambiguous_genes.index) + list(children_ambiguous_genes.index)
# gff_trimmed = gff.drop(idx)
# filename = "Data/Genomic_landscapes/GFF_parsed/" + set + "_trimmed.csv.gz"
# input.write_gff2csv(gff_trimmed, filename)


# Parse introns and utrs
gff_parsed = gff.gff.parse_attributes(infer_rank=True, parse_introns=True, parse_utr=True, n_cpus=16)
#gff_parsed = gff.gff.parse_attributes(infer_rank=True, parse_introns=True, parse_utr=False, n_cpus=16)

# Add flanking regions up/downstream (+- 1,2,3kb)
# Take all genes
# Replicate 6 times each gene
# Tag as upstream1,2,3/downstream1,2,3
# Change coordinates start/end
# # TODO consider strand == +-

# Strand +
genes = gff_parsed.loc[(gff_parsed["feature"] == "gene") & (gff_parsed["strand"] == "+")].copy()
flanking = pd.concat([genes]*6).sort_index()
seq = cycle(["upstream3kb", "upstream2kb", "upstream1kb", "downstream1kb", "downstream2kb", "downstream3kb"])
flanking['feature'] = [next(seq) for count in range(flanking.shape[0])]
flanking["parent"] = flanking["id"]

flanking_tmp = flanking.copy()
# # not elegant but works
flanking.loc[flanking['feature'] == "upstream3kb", "end"] = flanking_tmp.loc[flanking_tmp['feature'] == "upstream3kb", "start"] - 2001
flanking.loc[flanking['feature'] == "upstream2kb", "end"] = flanking_tmp.loc[flanking_tmp['feature'] == "upstream2kb", "start"] - 1001
flanking.loc[flanking['feature'] == "upstream1kb", "end"] = flanking_tmp.loc[flanking_tmp['feature'] == "upstream1kb", "start"]
flanking.loc[flanking['feature'] == "upstream3kb", "start"] = flanking_tmp.loc[flanking_tmp['feature'] == "upstream3kb", "start"] - 3000
flanking.loc[flanking['feature'] == "upstream2kb", "start"] = flanking_tmp.loc[flanking_tmp['feature'] == "upstream2kb", "start"] - 2000
flanking.loc[flanking['feature'] == "upstream1kb", "start"] = flanking_tmp.loc[flanking_tmp['feature'] == "upstream1kb", "start"] - 1000

flanking.loc[flanking['feature'] == "downstream3kb", "start"] = flanking_tmp.loc[flanking_tmp['feature'] == "downstream3kb", "end"] + 2001
flanking.loc[flanking['feature'] == "downstream2kb", "start"] = flanking_tmp.loc[flanking_tmp['feature'] == "downstream2kb", "end"] + 1001
flanking.loc[flanking['feature'] == "downstream1kb", "start"] = flanking_tmp.loc[flanking_tmp['feature'] == "downstream1kb", "end"]
flanking.loc[flanking['feature'] == "downstream3kb", "end"] = flanking_tmp.loc[flanking_tmp['feature'] == "downstream3kb", "end"] + 3000
flanking.loc[flanking['feature'] == "downstream2kb", "end"] = flanking_tmp.loc[flanking_tmp['feature'] == "downstream2kb", "end"] + 2000
flanking.loc[flanking['feature'] == "downstream1kb", "end"] = flanking_tmp.loc[flanking_tmp['feature'] == "downstream1kb", "end"] + 1000
flanking.loc[flanking["start"] < 0, "start"] = 0
flanking.loc[flanking["end"] < 0, "end"] = 0

# gff_cds_flank = pd.concat([gff_parsed, flanking])
gff_cds_flank = gff_parsed.append(flanking).copy()
gff_cds_flank = gff_cds_flank.reset_index(drop=True)
gff_cds_flank.loc[gff_cds_flank["start"] <= 0, "start"] = 1
gff_cds_flank.loc[gff_cds_flank["end"] <= 0, "end"] = 1

# Strand -
genes = gff_parsed.loc[(gff_parsed["feature"] == "gene") & (gff_parsed["strand"] == "-")].copy()
flanking = pd.concat([genes]*6).sort_index()
seq = cycle(["upstream3kb", "upstream2kb", "upstream1kb", "downstream1kb", "downstream2kb", "downstream3kb"])
flanking['feature'] = [next(seq) for count in range(flanking.shape[0])]
flanking["parent"] = flanking["id"]

flanking_tmp = flanking.copy()
# # not elegant but works
flanking.loc[flanking['feature'] == "upstream3kb", "end"] = flanking_tmp.loc[flanking_tmp['feature'] == "upstream3kb", "end"] + 3000
flanking.loc[flanking['feature'] == "upstream2kb", "end"] = flanking_tmp.loc[flanking_tmp['feature'] == "upstream2kb", "end"] + 2000
flanking.loc[flanking['feature'] == "upstream1kb", "end"] = flanking_tmp.loc[flanking_tmp['feature'] == "upstream1kb", "end"] + 1000
flanking.loc[flanking['feature'] == "upstream3kb", "start"] = flanking_tmp.loc[flanking_tmp['feature'] == "upstream3kb", "end"] + 2001
flanking.loc[flanking['feature'] == "upstream2kb", "start"] = flanking_tmp.loc[flanking_tmp['feature'] == "upstream2kb", "end"] + 1001
flanking.loc[flanking['feature'] == "upstream1kb", "start"] = flanking_tmp.loc[flanking_tmp['feature'] == "upstream1kb", "end"]

flanking.loc[flanking['feature'] == "downstream3kb", "start"] = flanking_tmp.loc[flanking_tmp['feature'] == "downstream3kb", "start"] - 3000
flanking.loc[flanking['feature'] == "downstream2kb", "start"] = flanking_tmp.loc[flanking_tmp['feature'] == "downstream2kb", "start"] - 2000
flanking.loc[flanking['feature'] == "downstream1kb", "start"] = flanking_tmp.loc[flanking_tmp['feature'] == "downstream1kb", "start"] - 1000
flanking.loc[flanking['feature'] == "downstream3kb", "end"] = flanking_tmp.loc[flanking_tmp['feature'] == "downstream3kb", "start"] - 2001
flanking.loc[flanking['feature'] == "downstream2kb", "end"] = flanking_tmp.loc[flanking_tmp['feature'] == "downstream2kb", "start"] - 1001
flanking.loc[flanking['feature'] == "downstream1kb", "end"] = flanking_tmp.loc[flanking_tmp['feature'] == "downstream1kb", "start"]
flanking.loc[flanking["start"] < 0, "start"] = 0
flanking.loc[flanking["end"] < 0, "end"] = 0

#gff_cds_flank = pd.concat([gff_cds_flank, flanking])

gff_cds_flank = gff_parsed.append(flanking).copy()
gff_cds_flank = gff_cds_flank.reset_index(drop=True)
gff_cds_flank.loc[gff_cds_flank["start"] <= 0, "start"] = 1
gff_cds_flank.loc[gff_cds_flank["end"] <= 0, "end"] = 1

# Save results
filename = "Data/Genomic_landscapes/GFF_parsed/" + dataset + "_parsed.csv.gz"
#gff_parsed.to_csv(filename, sep="\t", compression="gzip", index=False)
input.write_gff2csv(gff_cds_flank, filename)

print("End of GFF parsing: success.")
