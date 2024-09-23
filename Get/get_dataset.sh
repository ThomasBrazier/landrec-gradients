#!/bin/bash
remote_path=$(cat Get/remote_path)

scp $remote_path/Data/Recombination/Gradient/gff_rho_all.rds $PWD/Data/Recombination/Gradient/
scp $remote_path/Data/Recombination/rho_hotspot_center.*.rds $PWD/Data/Recombination/

#scp $remote_path/Data/Recombination/Gradient/gff_rho_all.rds $PWD/Data/Recombination/Gradient/
#scp $remote_path/Data/Recombination/Gradient/meanRho_rank.rds $PWD/Data/Recombination/Gradient/
#scp $remote_path/Data/Recombination/Gradient/meanRho_rank_nbexons.rds $PWD/Data/Recombination/Gradient/
#scp $remote_path/Data/Genomic_landscapes/nb_exons.rds $PWD/Data/Genomic_landscapes/
#scp $remote_path/Data/Recombination/Hotspot/hotspot_gradient_rank.rds $PWD/Data/Recombination/Hotspot/
#scp $remote_path/Data/Recombination/Hotspot/hotspot_gradient_rank_filtered.rds $PWD/Data/Recombination/Hotspot/
scp $remote_path/Data/Recombination/weightedmeanrho_per_feature.rds $PWD/Data/Recombination/
scp $remote_path/Data/Recombination/meanrho_per_feature.rds $PWD/Data/Recombination/
scp $remote_path/Data/Recombination/LD/ldhat/LD_maps.csv.gz $PWD/Data/Recombination/LD/ldhat/
scp $remote_path/Data/Recombination/LD/ldhat/LD_widths.csv.gz $PWD/Data/Recombination/LD/ldhat/


scp $remote_path/Data/Recombination/ldmap_windows_100kb.rds $PWD/Data/Recombination/
scp $remote_path/Data/Recombination/ldmap_windows_10kb.rds $PWD/Data/Recombination/

scp $remote_path/Data/Recombination/Gradient/meanRho*.rds $PWD/Data/Recombination/Gradient/
scp $remote_path/Data/Recombination/Gradient/gradient*.rds $PWD/Data/Recombination/Gradient/

scp $remote_path/Output/df_meanrho.rds $PWD/Data/Output/

scp $remote_path/Data/Recombination/LD/ldhotspots_raw.csv.gz $PWD/Data/Recombination/LD/