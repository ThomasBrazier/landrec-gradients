#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --mem=60GB
#SBATCH --cpus-per-task=8
#SBATCH --time=25-60:00:00
#SBATCH --job-name=ld_windows

. /local/env/envconda.sh

echo Smooth LD maps in 10 ad 100 kb windows

conda activate Renv

Rscript Jobs/ldmaps_windows_10_100kb.R

conda deactivate
