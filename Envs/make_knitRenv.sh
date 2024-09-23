#!/bin/bash
conda env create -f knitRenv.yaml

conda activate knitRenv

Rscript -e 'devtools::install_github("daijiang/phyr")'

conda deactivate
