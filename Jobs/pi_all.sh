#!/bin/bash

while read -r line;
do
   echo "$line"
   sbatch -p genouest,ecobio Jobs/pi.sh $line
done < Data/list_sets
