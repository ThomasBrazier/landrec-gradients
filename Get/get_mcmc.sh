#!/bin/bash
remote_path=$(cat Get/remote_path_scratch)

rm $PWD/Data/Recombination/LD/ldhat_mcmc/*

scp $remote_path/ldhat-snakemake-pipeline/data/*/K*.pop*/MCMC/*.html $PWD/Data/Recombination/LD/ldhat_mcmc/
scp $remote_path/ldhat-snakemake-pipeline/data/*/K*.pop*/MCMC/*.jpeg $PWD/Data/Recombination/LD/ldhat_mcmc/

