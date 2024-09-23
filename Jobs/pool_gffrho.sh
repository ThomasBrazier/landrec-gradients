#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --mem=120GB
#SBATCH --cpus-per-task=16
#SBATCH --time=25-60:00:00
#SBATCH --job-name=poolGFFRho

. /local/env/envconda.sh

echo Pool GFF Rho results

conda activate Renv

Rscript Jobs/pool_gffrho.R

conda deactivate
