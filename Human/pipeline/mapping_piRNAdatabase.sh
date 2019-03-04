#!/bin/bash


#Generated piRNA database in splin02
#cat piRBase/hsa_V2.fa > human_piRNA_database.f
#cat piRNAQuest/hsa_pirq.fa >> human_piRNA_database.fa
#awk '{if(NR%3!=0){print}}' piRNAQuest/human_NCBI_piRNA.fa >> human_piRNA_database.fa
#bowtie-build --thread 8 human_piRNA_database.fa bt-index/human_piRNA_database

piRDatabase=/data/usrhome/LabSPLin/splin02/Genomes/database/piRNA_database/bt-index/human_piRNA_database
inputDir=$HOME/human_gonad_piRNA/bed/20140214/bam_to_fastq

#Map piRNA candidates identified to piRNA database

for i in $(ls $inputDir|grep 'fastq'|sed 's/_piRNA_candidate.fastq//g')

do
	echo Start Mapping $i piRNA candidates
	bowtie -q --norc -S -v 0 -p 16 --al $i"_piRNA.fastq" --un $i"_non_piRNA.fastq" $piRDatabase $inputDir/$i"_piRNA_candidate.fastq" > $i".sam" 2> $i"_log.txt"
done
