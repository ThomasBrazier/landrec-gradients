#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --mem=60GB
#SBATCH --cpus-per-task=4
#SBATCH --time=25-60:00:00
#SBATCH --job-name=poolLDhot

. /local/env/envconda.sh

echo Pool LDhot results

conda activate Renv

Rscript Jobs/pool_ldhot.R

conda deactivate
