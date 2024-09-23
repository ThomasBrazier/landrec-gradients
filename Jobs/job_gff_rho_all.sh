#!/bin/bash

while read -r line;
do
   echo "$line"
   sbatch -p genouest,ecobio Jobs/job_gff_rho.sh $line
done < Data/list_sets
