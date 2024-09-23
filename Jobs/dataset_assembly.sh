#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --cpus-per-task=16
#SBATCH --mem=200GB
#SBATCH --time=20-60:00:00
#SBATCH --job-name=df_assembly

echo Compute the final dataset from LDhat and LDhot outputs

. /local/env/envconda.sh 

conda activate Renv

Rscript Jobs/job_dataset_assembly_LD.R

conda deactivate
