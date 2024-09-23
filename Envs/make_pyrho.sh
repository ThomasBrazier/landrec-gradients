#!/bin/bash
# MUST EXECUTE within Envs/ or where you wish to install it
conda create -n pyrho
conda activate pyrho

conda install -c conda-forge -c bioconda hdf5 openssl gsl cython msprime numba cyvcf2 vcftools samtools bedtools

conda install -c conda-forge singularity

mkdir Pyrho
cd Pyrho

singularity pull smcpp.sif docker://terhorst/smcpp:latest

git clone https://github.com/popgenmethods/ldpop.git ldpop
pip install ldpop/

git clone https://github.com/popgenmethods/pyrho.git pyrho
pip install pyrho/
