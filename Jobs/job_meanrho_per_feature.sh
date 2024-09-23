#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=120GB
#SBATCH --time=20-60:00:00
#SBATCH --job-name=rho

#source ~/miniconda3/etc/profile.d/conda.sh
. /softs/local/env/envconda.sh
conda activate Renv

Rscript -e Source/meanrho_per_feature.R $PWD
Rscript -e Source/weightedmeanrho_per_feature.R $PWD

conda deactivate
