#!/bin/bash

echo Recompute hotspot overlap in genomic features

conda activate Renv

Rscript $PWD/Source/job_hotspot_overlap.R

conda deactivate

