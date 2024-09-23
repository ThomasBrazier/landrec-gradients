#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --mem=600GB
#SBATCH --cpus-per-task=4
#SBATCH --time=25-60:00:00
#SBATCH --job-name=pi

dataset=$1

. /local/env/envvcftools-0.1.16.sh

vcftools --gzvcf Data/Polymorphism/${dataset}/${dataset}.pop.vcf.gz --out Data/Polymorphism/${dataset}/${dataset}  --site-pi

. /local/env/envconda.sh

conda activate Renv

Rscript Jobs/pi.R $dataset

conda deactivate
