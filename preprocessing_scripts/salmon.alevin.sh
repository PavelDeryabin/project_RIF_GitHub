#!/bin/bash 

# salmon=1.9.0 was installed using conda
# 'SraAccList.txt' files were downloaded from NCBI GEO
# 'txp2gene.tsv' and index files were generated as descibed in salnom.index.sh 

cd [PATH_TO_DATASET]

while read id; do

	cd ./"$id"
	
	echo "Mapping reads for $id"
	salmon alevin \
	--libType ISR \
	--mates1 *1.fastq \
	--mates2 *2.fastq \
	--chromiumV3 \
	--index [PATH_TO_INDEX] \
	--threads 16 \
	--output [PATH_TO_DATASET]/quants/"$id" \
	--tgMap [PATH_TO_INDEX]/txp2gene.tsv \
	--expectCells 10000
	
	cd ..
	rm -r "$id"
	
done < SRR_Acc_List.txt
