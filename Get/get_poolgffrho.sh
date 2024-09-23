#!/bin/bash
remote_path=$(cat Get/remote_path)

scp $remote_path/Data/Recombination/Gradient/gff_rho_all.rds $PWD/Data/Recombination/Gradient/
#scp $remote_path/Data/Genomic_landscapes/Rho/gff_rho_*.rds $PWD/Data/Genomic_landscapes/Rho/
