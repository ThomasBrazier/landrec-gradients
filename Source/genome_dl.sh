#!/bin/sh

file=$1

#----------------------------------------------------------------------
# Download
#----------------------------------------------------------------------
# the genomes to download are presented in a txt file with three columns:
# Species_name
# Accession, the id to be used for a request to the database
# Database (Ensembl, NCBI, Phytozome), giving the method to use for downloading

while IFS= read -r line || [[ -n "$line" ]];
do
    #echo $line
    # If db = Ensembl
    if (echo $line) | grep 'ensembl'
        then echo 'ensembl'
            echo $line | awk '{print $2}' # accession
            species=$(echo $line | awk '{print $1}' | tr '[:upper:]' '[:lower:]')
            accession=$(echo $line | awk '{print $2}')
            version=$(echo $line | awk '{print $4}')
            echo "Version explicitly required: $version"
            dir=$(echo $line | awk '{print $1 "/" $2}')
            mkdir $species # Make a directory for the species
            mkdir $dir # Make a directory for the accession
            # Make sub-directories
            #   ../fasta: Fasta sequences
            mkdir $dir/fasta
            #mkdir $dir/fasta/dna
            #mkdir $dir/fasta/cds
            #mkdir $dir/fasta/cdna            
            #       ../../DNA
            #       ../../CDS
            #       ../../cDNA
            #   ../gff3
            mkdir $dir/gff3
            #   ../gtf: Gene sets for each species. These files include annotations of both coding and non-coding genes
            mkdir $dir/gtf
            # Now, download from FTP in each directory for each species
            # rsync [OPTION]... SRC [SRC]... DEST
            rsync --copy-links --recursive --times --verbose --delete rsync://ftp.ensemblgenomes.org/all/pub/current/plants/fasta/$species/dna/ $dir/fasta/
            #wget -r -P $dir/fasta/dna ftp://ftp.ensemblgenomes.org/pub/plants/current/fasta/$species/dna/
            #rsync --copy-links --recursive --times --verbose rsync://ftp.ensemblgenomes.org/all/pub/current/plants/fasta/$species/cdna/ $dir/fasta/cdna/
            #rsync --copy-links --recursive --times --verbose rsync://ftp.ensemblgenomes.org/all/pub/current/plants/fasta/$species/cds/ $dir/fasta/cds/
            rsync --copy-links --recursive --times --verbose --delete rsync://ftp.ensemblgenomes.org/all/pub/current/plants/gff3/$species/ $dir/gff3/
            rsync --copy-links --recursive --times --verbose --delete rsync://ftp.ensemblgenomes.org/all/pub/current/plants/gtf/$species/ $dir/gtf/
            # Build local database
            ## Concatenate all chromosomes
            cd $dir/fasta
            cat $(ls | grep -e 'dna.chromosome\|dna.nonchromosomal') > chromosomes.fa.gz
            gunzip chromosomes.fa.gz
            ## Making db
            makeblastdb -in chromosomes.fa -parse_seqids -blastdb_version 5 -title "Local_Blast_$accession" -dbtype nucl
            #rm chromosomes.fa.gz
            rm chromosomes.fa
            cd ../../..
    fi
    # If db = NCBI
    if (echo $line) | grep 'ncbi'
        then echo 'ncbi'
            echo $line | awk '{print $2}' # accession
            species=$(echo $line | awk '{print $1}'| tr '[:upper:]' '[:lower:]')
            accession=$(echo $line | awk '{print $2}')
            version=$(echo $line | awk '{print $4}')
            echo $version
            dir=$(echo $line | awk '{print $1 "/" $2}')
            mkdir $species # Make a directory for the species
            mkdir $dir # Make a directory for the accession
            # Make sub-directories
            #   ../fasta: Fasta sequences
            mkdir $dir/fasta
            #mkdir $dir/fasta/dna
            #mkdir $dir/fasta/cds
            #mkdir $dir/fasta/cdna            
            #       ../../DNA
            #       ../../CDS
            #       ../../cDNA
            #   ../gff3
            mkdir $dir/gff3
            #   ../gtf: Gene sets for each species. These files include annotations of both coding and non-coding genes
            mkdir $dir/gtf
            # Now, download from FTP in each directory for each species
            # rsync [OPTION]... SRC [SRC]... DEST

            # Build the path of the request
            # Use the RefSeq FTP, contains more information of better quality
            request=$(rsync rsync://ftp.ncbi.nlm.nih.gov/genomes/all/${accession:0:3}/${accession:4:3}/${accession:7:3}/${accession:10:3}/ | grep -e 'GC')
            # If there is multiple versions, you need to choose the version you want
            # Otherwise, if only one version (one element) take the first one

            # If version is manually given as argument
            if [ -n "$version" ]
            then
                echo "You asked for a version: $version"
                req=$version
            else
                if [ $(echo $request | awk '{print $10}' | wc -l) = 0 ]
                then echo "Only one version"
                     req=$(echo $request | awk '{print $5}')
                else
                    # If you want v1
                    if [ $(echo ${accession:14:1}) = 1 ]
                    then echo 'v1'
                         req=$(echo $request | awk '{print $5}')
                    fi
                    # Else, if you want v2
                    if [ $(echo ${accession:14:1}) = 2 ]
                    then echo 'v2'
                         req=$(echo $request | awk '{print $10}')
                    fi
                    if [ $(echo ${accession:14:1}) = 3 ]
                    then echo 'v3'
                         req=$(echo $request | awk '{print $10}')
                    fi
                    if [ $(echo ${accession:14:1}) = 4 ]
                    then echo 'v4'
                         req=$(echo $request | awk '{print $10}')
                    fi
                fi
            fi


            rsync --copy-links --recursive --times --verbose --delete rsync://ftp.ncbi.nlm.nih.gov/genomes/all/${accession:0:3}/${accession:4:3}/${accession:7:3}/${accession:10:3}/$req/${req}_genomic.fna.gz $dir/fasta/

            # DOWNLOAD THE GFF ANNOTATION FILE IF IT EXISTS, OTHERWISE RETRIEVE AND CONVERT THE GBFF FILE
            echo "Downloading the gff file..."
            rsync --copy-links --recursive --times --verbose --delete rsync://ftp.ncbi.nlm.nih.gov/genomes/all/${accession:0:3}/${accession:4:3}/${accession:7:3}/${accession:10:3}/$req/${req}_genomic.gff.gz $dir/gff3/

            FILE=$dir/gff3/${req}_genomic.gff.gz
            if test -f "$FILE";
            then
                echo "${req}_genomic.gff.gz exists."
            else
                echo "${req}_genomic.gff.gz does not exist. Downloading the gbff file.."
                rsync --copy-links --recursive --times --verbose --delete rsync://ftp.ncbi.nlm.nih.gov/genomes/all/${accession:0:3}/${accession:4:3}/${accession:7:3}/${accession:10:3}/$req/${req}_genomic.gbff.gz $dir/gbff/
                # Convert gbff to gff
                # Clean gff dir before
                rm $dir/gff3/*
                # Retrieve the name of the genbank file
                gbname=$(sed "s/.gbff.gz//g" <<<$(ls $dir/gbff/))
                #perl /Users/tbrazier/Qsync/PhD/Analyses/sources/genbank2gff3.pl --dir $dir/genbank/*.gbf.gz --outdir $dir/gff3/ -z
                gunzip $dir/gbff/$gbname.gbff.gz
                echo "Convert gbff to gff3..."
                biopython.convert $dir/gbff/$gbname.gbff genbank $dir/gff3/$gbname.gff gff3
                # Gunzip gff
                gzip $dir/gff3/$gbname.gff
                gzip $dir/gbff/$gbname.gbff
                # Remove unzipped gbff in 'genbank'
                #rm -r $dir/gbff
            fi

            # Build local database
            echo "Building local database for blast."
            cd $dir/fasta
            gunzip -k *.fna.gz
            ## Making db
            makeblastdb -in *.fna -parse_seqids -blastdb_version 5 -title "Local_Blast_$accession" -dbtype nucl
            rm *.fna
            cd ../../..
    fi
done < $file

