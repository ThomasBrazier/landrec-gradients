#!/bin/bash
conda init bash
conda create -n herho && \
conda activate herho && \
conda install pandas numpy matplotlib && \
conda install -c conda-forge nlopt && \
pip install pybedtools tqdm docopt scikit-allel && \
conda install -c bioconda tabix && \
conda install -c defaults -c conda-forge -c bioconda bedops

mamba install -c conda-forge -c bioconda r-vcfr

cd Envs
git clone https://github.com/samebdon/heRho.git
