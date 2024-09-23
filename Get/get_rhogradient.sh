#!/bin/bash
remote_path=$(cat Get/remote_path)

scp $remote_path/Data/Recombination/Gradient/meanRho.rds $PWD/Data/Recombination/Gradient/
#scp $remote_path/Data/Recombination/Gradient/*.rds $PWD/Data/Recombination/Gradient/
#scp $remote_path/Output/RhoGradient_*.Rda $PWD/Output/
