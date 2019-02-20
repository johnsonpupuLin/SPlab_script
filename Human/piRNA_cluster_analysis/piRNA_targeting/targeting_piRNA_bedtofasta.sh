#!/bin/bash

#Program
#	Get the sense and antisense fasta file of piRNAs generated from specific clusters.
#Usage
#	Copy this file to where you want to output the result.
#	Modify "inputDir" to where your piRNA candidates bed file locating.
#	Modify "clusterDir" to where your specific cluster bed file locating.
#	Modify "refDir" and related files if you want to change species.
#	Execute.
#History
#2019/01/31, HHL
#	This script is made by KW. Modify and add description.



inputDir=$HOME/human_gonad_piRNA/bed/20130502
clusterDir=$HOME/human_gonad_piRNA/data/20190128/20130502/piRPKM/merged_cluster
refDir=/data/usrhome/LabSPLin/splin02/Genomes

for input in $(ls $inputDir|grep 'bed12'|sed 's/.bed12//g')
do 

echo Start $input
#Align whole piRNA candidates to specific clusters to identify the piRNAs which are generated from these clusters.
	intersectBed -split -wa -s -a $inputDir/$input".bed12" -b $clusterDir/$input"_clusters.split.bed"  > temp.intersectBed; 
	wait;

#Length filtering
	awk '{if(($3-$2)<=34){print}}' temp.intersectBed > temp2.intersectBed;
	wait;

#Get the origin sense sequences of piRNA candidates
	fastaFromBed -split -s -name -fi $refDir/UCSC/hg38/hg38.fa -bed temp2.intersectBed -fo $input"_sense.fasta"; 
	wait;

#Get the antisense sequences of piRNA candidates for targeting
	awk 'BEGIN{FS=OFS="\t"}{if($6=="+"){$6="-"}else{$6="+"}{print}}' temp2.intersectBed |
		fastaFromBed -split -s -name -fi $refDir/UCSC/hg38/hg38.fa -bed - -fo $input"_asense.fasta"; 
	wait; 

done 


rm temp*

