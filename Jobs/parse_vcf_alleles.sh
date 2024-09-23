#!/bin/bash

echo "Dataset   Chromosome  Position    Ref Alt Freq_missing" > Data/Polymorphism/vcf_alleles

while read -r line;
do
   echo "$line"
   bash Source/parse_vcf_alleles.sh $line
   zcat Data/Polymorphism/$line/$line.alleles.gz >> Data/Polymorphism/vcf_alleles
done < Data/list_sets

gzip Data/Polymorphism/vcf_alleles
