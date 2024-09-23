#!/bin/bash

remote_path=$(cat Get/remote_path)

source $PWD/Get/get_ldhat_output.sh

rm $PWD/Data/Recombination/LD/ldhat_mcmc/*

scp $remote_path/ldhat-snakemake-pipeline/data/*/K*.pop*/MCMC/*.html $PWD/Data/Recombination/LD/ldhat_mcmc/
scp $remote_path/ldhat-snakemake-pipeline/data/*/K*.pop*/MCMC/*.jpeg $PWD/Data/Recombination/LD/ldhat_mcmc/


scp $remote_path/Figure/Paper/* $PWD/Figure/Paper/
scp $remote_path/Table/* $PWD/Table/

rsync -ave ssh $remote_path/Figure/Landscapes/* $PWD/Figure/Landscapes/

scp $remote_path/Data/Recombination/LD/r2/* $PWD/Data/Recombination/LD/r2/

scp $remote_path/Data/Recombination/ldmap_windows_100kb.rds $PWD/Data/Recombination/
scp $remote_path/Data/Recombination/ldmap_windows_10kb.rds $PWD/Data/Recombination/

rsync -ave $remote_path/Data/Recombination/LD/ldhat/*.csv.gz $PWD/Data/Recombination/LD/ldhat/
rsync -ave $remote_path/Data/Recombination/LD/ldhat/LD_maps.csv.gz $PWD/Data/Recombination/LD/ldhat/
rsync -ave $remote_path/Data/Recombination/LD/ldhat/LD_widths.csv.gz $PWD/Data/Recombination/LD/ldhat/

rsync -ave $remote_path/Data/Recombination/LD/ldhotspots_raw.csv.gz $PWD/Data/Recombination/LD/
rsync -ave $remote_path/Data/Recombination/LD/ldhotspots_filtered.csv.gz $PWD/Data/Recombination/LD/

rsync -ave $remote_path/Data/Recombination/Gradient/gff_rho_*.rds $PWD/Data/Recombination/
rsync -ave $remote_path/Data/Recombination/Gradient/gff_rho_all.rds $PWD/Data/Recombination/
rsync -ave $remote_path/Data/Genomic_landscapes/nb_exons.rds $PWD/Data/Genomic_landscapes/

rsync -ave $remote_path/Data/Recombination/Hotspot/dist2hotspot_filtered.rds $PWD/Data/Recombination/Hotspot/
rsync -ave $remote_path/Data/Recombination/Hotspot/dist2hotspot_raw.rds $PWD/Data/Recombination/Hotspot/

rsync -ave $remote_path/Data/Recombination/Gradient/meanRho.rds $PWD/Data/Recombination/Gradient/
rsync -ave $remote_path/Data/Recombination/Gradient/meanRho_rank.rds $PWD/Data/Recombination/Gradient/
rsync -ave $remote_path/Data/Recombination/Gradient/meanRho_rank_nbexons.rds $PWD/Data/Recombination/Gradient/
rsync -ave $remote_path/Data/Recombination/Gradient/meanRho_rank.rds $PWD/Data/Recombination/Gradient/
rsync -ave $remote_path/Data/Recombination/Gradient/medianRho.rds $PWD/Data/Recombination/Gradient/

rsync -ave $remote_path/Data/Recombination/Gradient/meanRho_hotoverlap.rds $PWD/Data/Recombination/Gradient/
rsync -ave $remote_path/Data/Recombination/Gradient/meanRho_hotoverlap_TTS.rds $PWD/Data/Recombination/Gradient/

rsync -ave $remote_path/Data/Recombination/Gradient/RhoGradient_ATG_*.rds $PWD/Data/Recombination/Gradient/
rsync -ave $remote_path/Data/Recombination/Gradient/RhoGradient_5kbTTS_*.rds $PWD/Data/Recombination/Gradient/
rsync -ave $remote_path/Data/Recombination/Gradient/RhoGradient_5kbTSS_*.rds $PWD/Data/Recombination/Gradient/
rsync -ave $remote_path/Data/Recombination/Gradient/gradient_5kbTSS_TTS.rds $PWD/Data/Recombination/Gradient/
rsync -ave $remote_path/Data/Recombination/Gradient/gradient_TTS.rds $PWD/Data/Recombination/Gradient/

rsync -ave $remote_path/Data/Recombination/Hotspot/hotspot_rank_nbexons.rds $PWD/Data/Recombination/Hotspot/
rsync -ave $remote_path/Data/Recombination/Hotspot/hotspot_rank_nbexons_filtered.rds $PWD/Data/Recombination/Hotspot/
rsync -ave $remote_path/Output/df_meanrho.rds $PWD/Output/

scp $remote_path/Data/Recombination/rho_hotspot_center.*.rds $PWD/Data/Recombination/


rsync -ave $remote_path/Output/gff_rho_*_*.Rda $PWD/Output/

