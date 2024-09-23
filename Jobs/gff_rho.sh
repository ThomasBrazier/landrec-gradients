#!/bin/bash
while read -r line;
do
   echo "$line"
   sbatch -p genouest,bigmem,ecobio Jobs/job_gff_rho.sh $line
done < Data/list_sets
