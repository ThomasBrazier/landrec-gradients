#!/bin/bash

# Parce a VCF file fo extract the ref and alt alleles at each position
# + missing data per site
. /local/env/envvcftools-0.1.16.sh

dataset=$1

rm Data/Polymorphism/$dataset/$dataset.alleles

zcat Data/Polymorphism/$dataset/$dataset.pop.vcf.gz | grep -v ^# | awk '{printf ("%s\t%s\t%s\t%s\t%s\n", ds, $1, $2, $4, $5) }' ds=$dataset > Data/Polymorphism/$dataset/$dataset.tmp
# Add missing data per site with vcftools --missing-site
vcftools --missing-site --gzvcf Data/Polymorphism/$dataset/$dataset.pop.vcf.gz --out Data/Polymorphism/$dataset/$dataset

cat Data/Polymorphism/$dataset/$dataset.lmiss | tail -n +2 | awk '{ print $6 }' > Data/Polymorphism/$dataset/$dataset.freqmissing

paste Data/Polymorphism/$dataset/$dataset.tmp Data/Polymorphism/$dataset/$dataset.freqmissing > Data/Polymorphism/$dataset/$dataset.alleles

rm Data/Polymorphism/$dataset/$dataset.tmp
rm Data/Polymorphism/$dataset/$dataset.freqmissing

gzip --force Data/Polymorphism/$dataset/$dataset.alleles

