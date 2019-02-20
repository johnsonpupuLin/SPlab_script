#!/bin/bash

#Program
#	Calculate the piRPKM of merged clusters in datasets.
#Usage
#	Copy this script to where you want to output the result.
#	Modify "clusterDir" to where merged cluster bed file locating.
#	Modify "inputDir" to where piRNA bed file locating.
#	Execute.
#History
#2019/01/29, HHL
#	First version.

inputDir=$HOME/human_gonad_piRNA/bed/20130502
clusterDir=

for input in $(ls $inputDir|grep 'bed12'|sed 's/.bed12//g')

do

#Calculating the reads number aligned to merged clusters.
	echo start $input;
	intersectBed -s -c -a $clusterDir/merged.cluster.bed -b $inputDir/$input.bed12 > $input"_cluster_piRPM";
	wait;
	echo finish $input;

# Counting piRPKM (piRPKM = (piRNA numReads)/(geneLength/1000 * totalNumReads/1000000))
	totalreads=`wc -l $inputDir/$input".bed12"`
	echo $totalreads
	awk -F "\t" 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,($7/((($3-$2)/1000)*(totalreads/1000000)))}' totalreads="$totalreads" $input"_cluster_piRPM" > $input"_cluster_piRPKM"

done

