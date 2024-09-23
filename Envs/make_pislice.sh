#!/bin/bash

python3 -m venv pislice_env
cd pislice_env

source bin/activate

git clone https://github.com/ThomasBrazier/PiSlice.git

cd PiSlice
python3 -m pip install -e .
cd ..

python3 -m pip install pandas numpy biopython intervaltree multiprocess cyvcf2 mapply scikit-allel pyfaidx egglib pysam


