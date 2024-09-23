#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --cpus-per-task=8
#SBATCH --mem=120GB
#SBATCH --time=20-60:00:00
#SBATCH --job-name=knitHotspots

echo How to define hotspots

. /local/env/envconda.sh 

conda activate knitRenv

Rscript -e 'library(rmarkdown); rmarkdown::render("define_hotspots.Rmd", "html_document")' 

conda deactivate
