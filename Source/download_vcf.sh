#!/bin/bash
globalpath=$1
dataset=$2

# DOWNLOADING RESEQUENCING DATA
# Recoded trimmed VCF is renamed after the population name in our study
# i.e. Genus_species_AuthorYear
echo "DOWNLOADING RESEQUENCING DATA - VCF FILES"
echo "Processing dataset $dataset"

echo "Creating directory"
cd $globalpath/data/polymorphism_data
mkdir -p $dataset
cd $dataset

link=$(cat ../vcf_links.csv | grep $dataset | grep .vcf.gz | awk '{print $2}' | tr "\t" ";")

echo "Downloading data"
cat ../vcf_links.csv | grep $dataset | while read line
do
	echo "Download $line..."
	wget -N $(echo $line | awk '{print $2}')
	# Some links do not have any extension instead of expected .vcf.gz
	filename=$(basename $(echo $line | awk '{print $2}'))
	ext="${filename##*.}"
	if [[ -z $ext ]]; then mv $filename $filename.vcf.gz; fi
done

# Some downloaded files have an extension after .vcf.gz
if [ $(ls -l | grep vcf.gz[A-Za-z0-9?]+ | wc -l) -ge 1  ]
then
	"Correcting filename extension"
	filename=$(ls | grep vcf.gz[A-Za-z0-9?]+ | cut -f 1 -d '.')
	mv $filename.vcf.gz* $filename.vcf.gz
fi

# Some files may be compressed as .tar archives
if [ $(ls -l | grep *.tar | wc -l) -ge 1  ]
then
	echo "Untar"
	tar --strip-components=1 -xf *.tar
fi

echo "Indexing vcf(s)"
# Index vcf files
gunzip *.vcf.gz
for f in *.vcf
do
	echo "... $f"
	bgzip $f
	tabix -p vcf $f.gz
done

# If multiple vcfs (e.g. splitted by chromosomes), merge them
if [ $(ls -l | grep -E "*.vcf.gz$" | wc -l) -gt 1  ]
then
	echo "Merge multiple vcfs"
	. /softs/local/env/envbcftools-1.9.sh
	bcftools concat -o $dataset.vcf *.vcf.gz
	#vcf-concat *.vcf.gz | gzip -c > out.vcf.gz
	# remove chromosome vcfs
	# For Zea_mays_QiSun2018
	#rm merged_flt_*
	rm *.vcf.gz
	rm *.vcf.gz.tbi
	bgzip $dataset.vcf
	tabix -p vcf $dataset.vcf.gz
fi

# Compute x -p vcf $f.gzSHA1 checksum
echo "Computing SHA1 checksum"
vcfchecksum=$(sha1sum *.vcf.gz | awk '{print $1}')
echo $vcfchecksum

# Rename vcf
mv *.vcf.gz $dataset.vcf.gz
mv *.vcf.gz.csi $dataset.vcf.gz.csi
# Individuals will be filtered in further analyses; during trimming in the LDhat pipeline
# User have to put a .ind file in the dataset directory to specify which individuals will be imported
wget -N https://thomasbrazier.github.io/dev/Polymorphism_data/ind/$dataset.ind

# Standardization of chromosomes names
# Chromosome names must be integers
echo "Correcting chromosome names"
bgzip -d $dataset.vcf.gz
# Malus sp has chromosome naming like Chr01
sed -i -E "s/Chr0?//g" $dataset.vcf
# Citrullus lanatus Guo 2019
sed -i -E "s/Cla97Chr0?//g" $dataset.vcf
# Brachypodium distachyon Stritt 2017
sed -i -E "s/^Bd//g" $dataset.vcf
bgzip $dataset.vcf

# Compute summary statistics for quality control
#echo "Computing summary statistics for quality control"
#vcftools --gzvcf *.vcf.gz --missing-indv
#vcftools --gzvcf *.vcf.gz --missing-sites
#vcftools --gzvcf *.vcf.gz --site-quality
#vcftools --gzvcf *.vcf.gz --depth
#vcftools --gzvcf *.vcf.gz --site-depth
#vcftools --gzvcf *.vcf.gz --site-mean-depth
# Use 'vcfstats' for further analysis
# https://github.com/pwwang/vcfstats

# Clean up
rm *.vcf


# Save logs
cd $globalpath/data/polymorphism_data
echo "$dataset\t$link\t$vcfchecksum\t$(date -u)\t$(whoami)" >> polymorphism_data.csv
echo "End of download for $dataset"

