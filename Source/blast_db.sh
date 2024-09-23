#!/usr/sh


# Build local BLAST database
echo "Building local database for blast."
cd $dir/fasta
gunzip -k *.fna.gz
## Making db
makeblastdb -in *.fna -parse_seqids -blastdb_version 5 -title "Local_Blast_$accession" -dbtype nucl
rm *.fna
cd ../../..
