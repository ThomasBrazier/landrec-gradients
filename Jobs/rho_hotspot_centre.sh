#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --cpus-per-task=1
#SBATCH --mem=120GB
#SBATCH --time=20-60:00:00
#SBATCH --job-name=rhoHotspot

. /softs/local/env/envconda.sh
conda activate Renv

Rscript Jobs/rho_hotspot_centre.R

conda deactivate
