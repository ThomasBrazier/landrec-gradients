#!/bin/bash
remote_path=$(cat Get/remote_path)

mkdir -p $PWD/Data/Genomic_landscapes/Polymorphism
rsync -ave ssh $remote_path/Data/Genomic_landscapes/Polymorphism/* $PWD/Data/Genomic_landscapes/Polymorphism/
