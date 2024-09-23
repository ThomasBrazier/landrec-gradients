#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes.fr
#SBATCH --mail-type=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=220GB
#SBATCH --time=20-60:00:00
#SBATCH --job-name=rhodeciles

#source ~/miniconda3/etc/profile.d/conda.sh
. /softs/local/env/envconda.sh
conda activate Renv

Rscript Jobs/job_gradients_deciles.R

echo "==========================="
echo "End of script"
echo "==========================="

conda deactivate
