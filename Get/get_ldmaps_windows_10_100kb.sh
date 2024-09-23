#!/bin/bash
remote_path=$(cat Get/remote_path)

scp $remote_path/Data/Recombination/ldmap_windows_100kb.rds $PWD/Data/Recombination/
scp $remote_path/Data/Recombination/ldmap_windows_10kb.rds $PWD/Data/Recombination/

