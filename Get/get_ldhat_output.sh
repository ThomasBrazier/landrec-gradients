#!/bin/bash
remote_path=$(cat Get/remote_path)

scp $remote_path/Data/Recombination/LD/ldhat/* $PWD/Data/Recombination/LD/ldhat/
scp $remote_path/Data/Recombination/LD/ldhat/* $PWD/Data/Recombination/LD/ldhat/

#while IFS= read -r dataset
#do
#scp $remote_path/LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$dataset/K*.pop*/ldhat/*.res.txt.gz $PWD/Data/Recombination/LD/ldhat/
#scp $remote_path/LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$dataset/K*.pop*/ldhot/*.hot_summary.txt.gz $PWD/Data/Recombination/LD/ldhot/
#scp $remote_path/LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$dataset/K*.pop*/ldhot/*.hotspots.txt.gz $PWD/Data/Recombination/LD/ldhot/
#done < Data/list_sets


