#!/bin/bash

#Program
#       Generate fasta file from bed file, then make the graph of sequence logo via Weblogo.
#Usage
#	Modify variation "inputDir" to where the bed file you want to make sequence logo locating.
#       Modify variation "desDir" to where you want to store output file.
#       Modify for loop input to incorporate all file you want.
#       Execute.
#History
#2019/01/21, HHL
#       Add description.

inputDir=$HOME/human_gonad_piRNA/data/20190212_PD_sample/piRNA_origin/bed_file
desDir=./
refDir=/data/usrhome/LabSPLin/splin02/Genomes

for input in $(ls $inputDir|grep 'bed12'|sed 's/.bed12//g')

do
	echo Generate $input fasta file;
	fastaFromBed -fi $refDir/UCSC/hg38/hg38.fa -split -bed $inputDir/$input.bed12 -fo $desDir/$input.fasta;
	wait;
	echo $input fasta has been generated;

	echo Generate $input weblogo
	awk -F "\t" 'NR%2==0{$1=toupper(substr($1,1,20))}{print}' $desDir/$input.fasta| sed 's/T/U/g'|
	seqlogo -T 0.5 -E 1 -F PNG -w 30 -o $input -k 1 -c -e -n -Y -f -
done
