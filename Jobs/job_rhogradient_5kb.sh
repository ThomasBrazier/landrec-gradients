#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --cpus-per-task=16
#SBATCH --mem=320GB
#SBATCH --time=20-60:00:00
#SBATCH --job-name=gradient5kb

. /softs/local/env/envconda.sh
conda activate Renv

dataset=$1
# chromosome=$2

# Rscript Jobs/job_rhogradient_5kbATG.R $dataset $chromosome

# Run job on all chromosomes at once
list_chromosomes=$(cat Data/Genome/genome_chromosome_metadata.csv | grep $dataset | awk '{ print $5 }')

for chromosome in $list_chromosomes
do
if [ -f Data/Recombination/LD/ldhat/${dataset}.${chromosome}.bpen5.res.txt.gz ]
then
  if [ -f Data/Recombination/LD/ldhot/${dataset}.${chromosome}.bpen5.hot_summary.txt.gz ]
  then
    echo Dataset: $dataset
    echo Chromosome: $chromosome
    echo "ESTIMATE RHO GRADIENT AROUND ATG"
    Rscript Jobs/job_rhogradient_5kbATG.R $dataset $chromosome

    echo "ESTIMATE RHO GRADIENT AROUND TSS"
    Rscript Jobs/job_rhogradient_5kbTSS.R $dataset $chromosome

    echo "ESTIMATE RHO GRADIENT AROUND TTS"
    Rscript Jobs/job_rhogradient_5kbTTS.R $dataset $chromosome
  fi
fi
done

conda deactivate
