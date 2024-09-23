#!/bin/bash
remote_path=$(cat Get/remote_path)

#scp $remote_path/Data/Recombination/Gradient/gff_rho_all.rds $PWD/Data/Recombination/Gradient/
#scp $remote_path/Data/Recombination/Gradient/meanRho_rank.rds $PWD/Data/Recombination/Gradient/
#scp $remote_path/Data/Recombination/Gradient/meanRho_rank_nbexons.rds $PWD/Data/Recombination/Gradient/
#scp $remote_path/Data/Genomic_landscapes/nb_exons.rds $PWD/Data/Genomic_landscapes/
#scp $remote_path/Data/Recombination/Hotspot/hotspot_gradient_rank.rds $PWD/Data/Recombination/Hotspot/
scp $remote_path/Figure/Paper/* $PWD/Figure/Paper/
scp $remote_path/Table/* $PWD/Table/
#scp $remote_path/Data/Recombination/Hotspot/hotspot_gradient_rank_filtered.rds $PWD/Data/Recombination/Hotspot/
