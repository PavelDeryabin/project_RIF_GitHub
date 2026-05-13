#!/bin/bash 

# sra-tools=2.11.0 was installed using conda 
# 'SraAccList.txt' files were downloaded from NCBI GEO


cd [PATH_TO_DATASET]

while read id; do

	echo "Downloading .sra for $id"
	prefetch -p --max-size 10000000000000000 "$id"
	
	cd ./"$id"
	
	echo "Dumping .sra for $id"
	fasterq-dump --threads 16 -p --split-files *.sra
	rm *.sra
	
	cd ..
	
done < SRR_Acc_List.txt
