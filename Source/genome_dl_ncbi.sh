#!/bin/sh

# Dependencies
# biopython.convert
# pip install biopython-convert

species=$1
accession=$2

#----------------------------------------------------------------------
# Download
#----------------------------------------------------------------------
# the genomes to download are presented in a csv file with three columns:
# Species_name
# Accession, the id to be used for a request to the database
# Version

dir=$(echo "./Data/Genome/"$species"_"$accession)
mkdir -p $dir # Make a directory for the species
# Now, download from FTP
# rsync [OPTION]... SRC [SRC]... DEST
# Build the path of the request
# Use the RefSeq FTP, contains more information of better quality
request=$(rsync rsync://ftp.ncbi.nlm.nih.gov/genomes/all/${accession:0:3}/${accession:4:3}/${accession:7:3}/${accession:10:3}/ | grep -e 'GC')
# Multiple versions, get the one matching accession
version=$(echo $request | tr ' ' '\n' | grep $accession)

# Download metadata
rsync --copy-links --recursive --times --verbose --delete rsync://ftp.ncbi.nlm.nih.gov/genomes/all/${accession:0:3}/${accession:4:3}/${accession:7:3}/${accession:10:3}/$version/*.txt $dir/

# Download FASTA
echo "Downloading fasta..."
rsync --copy-links --recursive --times --verbose --delete rsync://ftp.ncbi.nlm.nih.gov/genomes/all/${accession:0:3}/${accession:4:3}/${accession:7:3}/${accession:10:3}/$version/${version}_genomic.fna.gz $dir/
# Change fasta name
mv -f $dir/${version}_genomic.fna.gz $dir/$species"_"$accession.fna.gz
# BGzip
gunzip $dir/$species"_"$accession.fna.gz && bgzip $dir/$species"_"$accession.fna


# DOWNLOAD THE GFF ANNOTATION FILE
echo "Downloading gff..."
rsync --copy-links --recursive --times --verbose --delete rsync://ftp.ncbi.nlm.nih.gov/genomes/all/${accession:0:3}/${accession:4:3}/${accession:7:3}/${accession:10:3}/$version/${version}_genomic.gff.gz $dir/
mv -f $dir/${version}_genomic.gff.gz $dir/$species"_"$accession.gff.gz

# Alternative strategy if gff is not available
# FILE=$dir/${version}_genomic.gff.gz
# if test -f "$FILE";
#   then
#      echo "${req}_genomic.gff.gz exists."
#   else
#      echo "${req}_genomic.gff.gz does not exist. Downloading the gbff file.."
#      rsync --copy-links --recursive --times --verbose --delete rsync://ftp.ncbi.nlm.nih.gov/genomes/all/${accession:0:3}/${accession:4:3}/${accession:7:3}/${accession:10:3}/$version/${version}_genomic.gbff.gz $dir/
#      # Convert gbff to gff
#      gunzip $dir/*.gbff.gz
#      echo "Convert gbff to gff3..."
#      biopython.convert $dir/*.gbff genbank $dir/$version.gff gff3
#      # Gunzip gff
#      gzip $dir/*.gff
#      gzip $dir/*.gbff
