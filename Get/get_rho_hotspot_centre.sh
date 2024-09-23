#!/bin/bash
remote_path=$(cat Get/remote_path)

scp $remote_path/Data/Recombination/rho_hotspot_center.*.rds $PWD/Data/Recombination/
