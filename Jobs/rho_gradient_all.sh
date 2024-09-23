#!/bin/bash

while read -r line;
do
   echo "$line"
   sbatch -p ecobio,genouest,bigmem Jobs/job_rhogradient_5kb.sh $line
done < Data/list_sets
