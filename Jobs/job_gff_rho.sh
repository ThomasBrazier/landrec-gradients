#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --cpus-per-task=16
#SBATCH --mem=120GB
#SBATCH --time=20-60:00:00
#SBATCH --job-name=gffrho

dataset=$1

#source ~/miniconda3/etc/profile.d/conda.sh
. /softs/local/env/envconda.sh
conda activate mamba_env
. ~/pislice_env/bin/activate

# Add a header with data/time
dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "==========================="
echo "$dataset"
echo "$dt"
echo "==========================="

echo "Parse the GFF and compute genomic statistics"

echo "==========================="
echo "Parse the <dataset>"
echo "==========================="
# python Jobs/job_gff_parse.py $dataset


echo "==========================="
echo "Add intergenic regions"
echo "==========================="
# Keep only protein-coding genes
# Add intergenic regions
# Add intergenic - buffer 5kb
# Add intergenic - buffer 10kb
# TODO Flag flanking overlapping a gene or another flanking (redundant)
# Rscript Jobs/job_add_intergenic.R $dataset
# Replace the gff_parsed with a new one with intergenic regions

cp -f LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$dataset/K*.pop*/$dataset.pop.vcf.gz* Data/Polymorphism/$dataset/

echo "==========================="
echo "Estimate GC content and other genomic statistics (SNP density, Pi) for each genomic feature in gff"
echo "==========================="
# Import the gff_parsed and output a new gff file with enhanced annotation and statistics
# python Jobs/job_gff_statistics.py $dataset
conda deactivate

echo "==========================="
echo "Compute GFF Rho for <dataset> and <chromosome>"
echo "==========================="
# Split chromosomes for computational issues
conda activate Renv
list_chromosomes=$(cat Data/Genome/genome_chromosome_metadata.csv | grep $dataset | awk '{ print $5 }')
for chromosome in $list_chromosomes
do
if [ -f Data/Recombination/LD/ldhat/${dataset}.${chromosome}.bpen5.res.txt.gz ]
then
  if [ -f Data/Recombination/LD/ldhot/${dataset}.${chromosome}.bpen5.hot_summary.txt.gz ]
  then
  Rscript Jobs/job_gff_rho.R $dataset $chromosome
  fi
fi
done

echo "==========================="
echo "End of script"
echo "==========================="

conda deactivate
