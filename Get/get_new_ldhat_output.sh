#!/bin/bash
remote_path=$(cat Get/remote_path_genossh)

while IFS= read -r dataset
do
# Polymorphism
scp $remote_path/LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$dataset/K*.pop*/$dataset.pop.vcf.gz* $PWD/Data/Polymorphism/$dataset/
scp $remote_path/LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$dataset/K*.pop*/poplist $PWD/Data/Polymorphism/$dataset/
scp $remote_path/LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$dataset/K*.pop*/$dataset.*.quality.html $PWD/Data/Polymorphism/$dataset/

scp $remote_path/LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$dataset/$dataset.vcf.gz* $PWD/Data/Polymorphism/$dataset/
scp $remote_path/LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$dataset/$dataset.popstatistics.* $PWD/Data/Polymorphism/$dataset/
scp $remote_path/LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$dataset/$dataset.quality.html $PWD/Data/Polymorphism/$dataset/
scp $remote_path/LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$dataset/indlist $PWD/Data/Polymorphism/$dataset/
scp $remote_path/LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$dataset/samplelist $PWD/Data/Polymorphism/$dataset/

# LDHAT
scp $remote_path/LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$dataset/K*.pop*/ldhat/*.res.txt.gz $PWD/Data/Recombination/LD/ldhat/
# MCMC and Quality
scp $remote_path/LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$dataset/K*.pop*/*.quality.html $PWD/Data/Recombination/LD/ldhat_mcmc/
scp $remote_path/LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$dataset/K*.pop*/MCMC/*.html $PWD/Data/Recombination/LD/ldhat_mcmc/
scp $remote_path/LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$dataset/K*.pop*/MCMC/*.jpeg $PWD/Data/Recombination/LD/ldhat_mcmc/
# LDHOT
scp $remote_path/LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$dataset/K*.pop*/ldhot/*.hot_summary.txt.gz $PWD/Data/Recombination/LD/ldhot/
scp $remote_path/LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$dataset/K*.pop*/ldhot/*.hotspots.txt.gz $PWD/Data/Recombination/LD/ldhot/
done < Data/list_sets