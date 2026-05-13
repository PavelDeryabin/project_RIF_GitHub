#!/bin/bash 

# salmon=1.9.0 was installed using conda
# GENCODE Human Release 41 (GRCh38.p13) files were downloaded from https://www.gencodegenes.org/human/release_41.html


# index was generated as follows

wget ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.transcripts.fa.gz
wget ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz

grep "^>" <(gunzip -c GRCm38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt

sed -i.bak -e 's/>//g' decoys.txt

cat gencode.v41.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gentrome.fa.gz

salmon index -t gentrome.fa.gz -d decoys.txt -p 16 -i [PATH_TO_INDEX] --gencode


# txp2gene was generated as follows

wget ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.primary_assembly.annotation.gtf.gz

bioawk -c gff '$feature=="transcript" {print $group}' gencode.v41.primary_assembly.annotation.gtf.gz | awk -F ' ' '{print substr($4,2,length($4)-3) "\t" substr($2,2,length($2)-3)}' > [PATH_TO_INDEX]/txp2gene.tsv
