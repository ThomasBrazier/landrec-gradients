#!/bin/bas

while read -r line;
do
   echo "$line"
   scp -v  LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$line/K*.pop*/${line}.pop.vcf.gz  Data/Polymorphism/${line}.pop.vcf.gz
   scp -v  LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$line/K*.pop*/${line}.pop.vcf.gz.csi  Data/Polymorphism/${line}.pop.vcf.gz.csi
   scp -v  LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$line/K*.pop*/${line}.pop.vcf.gz.tbi  Data/Polymorphism/${line}.pop.vcf.gz.tbi
   #rsync -rv  LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$line/${line}.vcf.gz  Data/Polymorphism/${line}.vcf.gz
   #rsync -rv  LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$line/${line}.vcf.gz.csi  Data/Polymorphism/${line}.vcf.gz.csi
   #rsync -rv  LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$line/${line}.vcf.gz.tbi  Data/Polymorphism/${line}.vcf.gz.tbi
done < Data/list_sets
