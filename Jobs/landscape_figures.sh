#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --mem=60GB
#SBATCH --cpus-per-task=4
#SBATCH --time=25-60:00:00
#SBATCH --job-name=fig_maps

. /local/env/envconda.sh

echo Generate LD landscapes of all chromosomes

conda activate Renv

Rscript Jobs/landscape_figures.R

conda deactivate
