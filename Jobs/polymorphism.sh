#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --mem=200GB
#SBATCH --cpus-per-task=16
#SBATCH --time=25-60:00:00
#SBATCH --job-name=polymorphism

# . /local/env/envvcftools-0.1.16.sh
# 
# while read -r line;
# do
#    echo "$line"
#    vcftools --gzvcf Data/Polymorphism/${line}/${line}.pop.vcf.gz --out Data/Polymorphism/${line}/${line}  --site-pi
# done < Data/list_sets


. /local/env/envconda.sh

conda activate Renv

Rscript Jobs/polymorphism.R

conda deactivate
