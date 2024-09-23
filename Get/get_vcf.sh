#!/bin/bash

while read -r line;
do
   echo "$line"
   scp LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$line/${line}.vcf.gz Data/Polymorphism/${line}/
   scp LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$line/K*.pop*/${line}.pop.vcf.gz Data/Polymorphism/${line}/
   scp LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$line/K*.pop*/${line}.ldhat.vcf.gz Data/Polymorphism/${line}/
   scp LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$line/K*.pop*/poplist Data/Polymorphism/${line}/
   scp LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$line/K*.pop*/subsetpop Data/Polymorphism/${line}/
done < Data/list_sets

