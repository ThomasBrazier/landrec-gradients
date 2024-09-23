#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes.fr
#SBATCH --mail-type=all
#SBATCH --mem=100GB
#SBATCH --cpus-per-task=16
#SBATCH --time=25-60:00:00
#SBATCH --job-name=PiPerSite

. /local/env/envconda.sh

source activate vcftools

# Calculate pi with monomorphic sites. No maf filter but rm indels and multiallelics. 
# Calculate pi also for the LDhat dataset

while read -r line;
do
   echo "$line"
   #Raw data, all individuals of the pop
   vcftools --gzvcf Data/Polymorphism/${line}/${line}.vcf.gz --keep Data/Polymorphism/${line}/poplist --remove-indels --max-alleles 2 --out Data/Polymorphism/${line}/${line}.noMAF --site-pi
   vcftools --gzvcf Data/Polymorphism/${line}/${line}.vcf.gz --keep Data/Polymorphism/${line}/poplist --remove-indels --max-alleles 2 --out Data/Polymorphism/${line}/${line}.noMAF --site-quality
   vcftools --gzvcf Data/Polymorphism/${line}/${line}.vcf.gz --keep Data/Polymorphism/${line}/poplist --remove-indels --max-alleles 2 --out Data/Polymorphism/${line}/${line}.noMAF --hardy
   
   # Raw data, only subset of the pop used for LDhat
   vcftools --gzvcf Data/Polymorphism/${line}/${line}.vcf.gz --keep Data/Polymorphism/${line}/subsetpop --remove-indels --max-alleles 2 --out Data/Polymorphism/${line}/${line}.noMAF.subset --site-pi
   vcftools --gzvcf Data/Polymorphism/${line}/${line}.vcf.gz --keep Data/Polymorphism/${line}/subsetpop --remove-indels --max-alleles 2 --out Data/Polymorphism/${line}/${line}.noMAF.subset --site-quality
   vcftools --gzvcf Data/Polymorphism/${line}/${line}.vcf.gz --keep Data/Polymorphism/${line}/subsetpop --remove-indels --max-alleles 2 --out Data/Polymorphism/${line}/${line}.noMAF.subset --hardy

   # Population data after filtering
   vcftools --gzvcf Data/Polymorphism/${line}/${line}.pop.vcf.gz --out Data/Polymorphism/${line}/${line}.pop --site-pi
   vcftools --gzvcf Data/Polymorphism/${line}/${line}.pop.vcf.gz --out Data/Polymorphism/${line}/${line}.pop --hardy
   vcftools --gzvcf Data/Polymorphism/${line}/${line}.pop.vcf.gz --out Data/Polymorphism/${line}/${line}.pop --site-quality

   # Subset of the population data used for LDhat, and after phasing and pseudodiploid
   vcftools --gzvcf Data/Polymorphism/${line}/${line}.ldhat.vcf.gz --out Data/Polymorphism/${line}/${line}.pop --site-pi
   vcftools --gzvcf Data/Polymorphism/${line}/${line}.ldhat.vcf.gz --out Data/Polymorphism/${line}/${line}.pop --hardy
   vcftools --gzvcf Data/Polymorphism/${line}/${line}.ldhat.vcf.gz --out Data/Polymorphism/${line}/${line}.pop --site-quality
done < Data/list_sets

conda deactivate
