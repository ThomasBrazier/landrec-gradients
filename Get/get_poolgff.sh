#!/bin/bash
remote_path=$(cat Get/remote_path)

scp $remote_path/Data/Genome/gff_all.rds $PWD/Data/Genome/
scp $remote_path/Data/Genomic_landscapes/GFF_parsed/*.csv.gz $PWD/Data/Genomic_landscapes/GFF_parsed/

