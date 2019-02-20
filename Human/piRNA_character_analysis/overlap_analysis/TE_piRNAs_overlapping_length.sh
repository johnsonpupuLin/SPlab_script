#!/bin/bash

#Program
#	Counting the overlapping length distribution of pairable piRNA generated from repeat region.
#Usage
#       Link this script to where you want to store your output data.
#	Modify variation "inputDir" to where uniq piRNA bed file locating. (Generated via uniq_and_strand.sh)
#	Modify for loop input to incorporate all file you want.
#	Execute.
#History
#2018/11/17, HHL
#	First version
#2019/01/22, HHL
#	Add description.

refDir=/data/usrhome/LabSPLin/splin02/Genomes
inputDir=$HOME/human_gonad_piRNA/bed/20130502

for input in $(ls $inputDir|grep 'bed12'|sed 's/.bed12//g')
do

#Isolated piRNAs generated from repeat regions
	echo Isolated sample: $input piRNAs generated from repeat regions
	intersectBed -s -u -wa -a $inputDir/$input.uniq.plus.bed -b $refDir/UCSC/hg38/RepeatMasker/hg38_RepeatMasker_annotated.bed > plus.temp.bed
	intersectBed -s -u -wa -a $inputDir/$input.uniq.minus.bed -b $refDir/UCSC/hg38/RepeatMasker/hg38_RepeatMasker_annotated.bed > minus.temp.bed

#Start overlapping piRNAs
	echo Start aligning $input + and - reads.
	intersectBed -wo -a plus.temp.bed -b minus.temp.bed > temp.bed
	cat temp.bed|awk '$8<=$2&&$2<=$9&&$9<=$3'|awk '{if ($13 <= 34){print $13, "\t", $5*$11}}' |sort -n >overlap.bed
	wait;

#Start counting length distribution
	echo Start counting overlapping reads amounts.

	for ((i=1;i<=34;i++));
	do
		cat overlap.bed|awk '{if ($1==i){print}}' i=$i|awk '{sum+=$2}END{print i, "\t", sum}' i=$i >> $input"_length_distribution"
	done
	wait;

	echo $input counting is finished, delete intermediates...
	rm *temp.bed overlap.bed
done
