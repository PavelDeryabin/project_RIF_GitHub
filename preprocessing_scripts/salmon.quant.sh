#!/bin/bash 

# salmon=1.9.0 was installed using conda
# 'SraAccList.txt' files were downloaded from NCBI GEO
# index file was generated as descibed in salnom.index.sh 

cd [PATH_TO_DATASET]

while read id; do

	cd ./"$id"
	
	echo "Mapping reads for $id"
	salmon quant \
	--index [PATH_TO_INDEX] \
	--libType A \
	--mates1 *1.fastq \
	--mates2 *2.fastq \
	--threads 16 \
	--output [PATH_TO_DATASET]/quants/"$id" \
	--seqBias \
	--gcBias \
	--validateMappings \
	--writeUnmappedNames
	
	cd ..
	rm -r "$id"
	
done < SRR_Acc_List.txt
