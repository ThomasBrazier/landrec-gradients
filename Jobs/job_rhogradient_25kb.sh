#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --cpus-per-task=16
#SBATCH --mem=120GB
#SBATCH --time=20-60:00:00
#SBATCH --job-name=grad25kb

. /softs/local/env/envconda.sh
conda activate Renv

dataset=$1
chromosome=$2

Rscript Jobs/job_rhogradient_25kbATG.R $dataset $chromosome
Rscript Jobs/job_rhogradient_25kbTSS.R $dataset $chromosome
Rscript Jobs/job_rhogradient_25kbTTS.R $dataset $chromosome

conda deactivate
