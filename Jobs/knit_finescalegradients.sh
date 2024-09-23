#!/bin/bash
wd=$PWD
ncores=8

echo Compute the final dataset from LDhat and LDhot outputs

. /local/env/envconda.sh 

conda activate knitRenv

cd Report/

Rscript -e 'library(rmarkdown); rmarkdown::render("fine_scale_gradients.Rmd", "html_document")' 
#Rscript ../Report/fine_scale_gradients.Rmd $wd

conda deactivate
