#!/bin/bash
dataset=$1
path="Data/Polymorphism"

# DOWNLOADING RESEQUENCING DATA
# Recoded trimmed VCF is renamed after the population name in our study
# i.e. Genus_species_AuthorYear
echo "DOWNLOADING RESEQUENCING DATA - VCF FILES"
echo "Processing dataset $dataset"

echo "Creating directory"
mkdir -p $path/$dataset

cd $path/$dataset

echo "Downloading data"
cat ../vcf_links.csv | grep $dataset | while read line
do
	echo "Download $line..."
	wget -N $(echo $line | awk '{print $2}') --no-check-certificate
	# Some links do not have any extension instead of expected .vcf.gz
	filename=$(basename $(echo $line | awk '{print $2}'))
	ext="${filename##*.}"
	if [[ -z $ext ]]; then mv $filename $filename.vcf.gz; fi
done

# Copy metadata in a sub-directory
mkdir metadata
cp *.doc *.docx *.xls *.xlsx *.txt *.csv metadata/

# Some downloaded files have an extension after .vcf.gz
if [ $(ls -l | grep vcf.gz[A-Za-z0-9?]+ | wc -l) -ge 1  ]
then
	"Correcting filename extension"
	filename=$(ls | grep vcf.gz[A-Za-z0-9?]+ | cut -f 1 -d '.')
	mv $filename.vcf.gz* $filename.vcf.gz
fi

unrar e *.rar
tar --strip-components=1 -xf *.tar
# Some files may be compressed as .tar archives
if [ $(ls -l | grep *.tar | wc -l) -ge 1  ]
then
	echo "Untar"
	tar --strip-components=1 -xf *.tar
fi

# If multiple vcfs (e.g. splitted by chromosomes), merge them
#if [ $(ls -l | grep -E "*.vcf.gz$" | wc -l) -gt 1  ]; then
echo "Merge multiple vcfs"
gunzip *.vcf.gz
for F in *.vcf ; do bgzip "${F}" done	
for F in *.vcf.gz ; do tabix -f -p vcf "${F}" done
bcftools concat -o $dataset.vcf *.vcf.gz
rm *.vcf.gz
rm *.vcf.gz.tbi
#fi

mv *.vcf $dataset.vcf
mv *.vcf.gz $dataset.vcf.gz
gunzip $dataset.vcf.gz

echo "Indexing vcf"
# Index vcf files
bgzip $dataset.vcf
tabix -p vcf $dataset.vcf.gz


# Compute x -p vcf $f.gzSHA1 checksum
echo "Computing SHA1 checksum"
vcfchecksum=$(sha1sum $dataset.vcf.gz | awk '{print $1}')
echo $vcfchecksum > metadata/md5checksum

echo "End of download for $dataset"

