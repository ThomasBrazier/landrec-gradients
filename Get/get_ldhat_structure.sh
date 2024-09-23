#!/bin/bash
remote_path=$(cat Get/remote_path)

while IFS= read -r dataset
do
rsync -ave ssh $remote_path/LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$dataset/structure/distruct.*.svg $PWD/Data/Recombination/LD/structure/$dataset/
rsync -ave ssh $remote_path/LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$dataset/structure/chooseK $PWD/Data/Recombination/LD/structure/$dataset/
rsync -ave ssh $remote_path/LDhat-snakemake-pipeline/ldhat-snakemake-pipeline/data/$dataset/$dataset.popstatistics.* $PWD/Data/Recombination/LD/structure/$dataset/
done < Data/list_sets

