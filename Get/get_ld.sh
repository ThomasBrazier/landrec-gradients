#!/bin/bash
remote_path=$(cat Get/remote_path)

#scp $remote_path/Data/Recombination/LD/ldhotspots_filtered.csv.gz $PWD/Data/Recombination/LD/
scp $remote_path/Data/Recombination/LD/r2/* $PWD/Data/Recombination/LD/r2/

