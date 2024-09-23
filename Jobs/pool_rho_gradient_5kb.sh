#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --mem=60GB
#SBATCH --cpus-per-task=4
#SBATCH --time=25-60:00:00
#SBATCH --job-name=poolGradients

. /local/env/envconda.sh

echo Pool gradients of recombination

conda activate Renv

Rscript Jobs/pool_rho_gradient_5kb.R

Rscript Jobs/pool_rho_gradient_5kb_median.R

conda deactivate
