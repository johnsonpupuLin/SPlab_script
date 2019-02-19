#!/bin/bash

#Program
#	Calculating the proportion of piRNA originated among different RefSeq regions.
#Usage
#	Link this script to where you want to store output file.
#       Modify variation "inputdir" to where piRNA candidate bed file locating. (Generated via mapping SP lab pipeline)
#       Modify variation "refdir" to where the RefSeq file locating.
#       Modify for loop input to incorporate all file you want.
#       Execute.
#History
#2018/11/07, HHL
#	Modified for human fetal gonad dataset
#2019/01/23, HHL
#	Add description.

#Gene composition: whole gene = intron + exons; exons= 5'UTR exons + conding exons + 3'UTR exons

inputdir=/data/usrhome/LabSPLin/splin01/human_gonad_piRNA/bed/20130502
refdir=/data/usrhome/LabSPLin/splin02/Genomes/UCSC/hg38/

for input in $(ls $inputdir/|grep '.bed12'|sed 's/.bed12//g')
do
	echo "Input $input small RNA-seq data and analyzing the distribution from 5'UTR to 3'UTR among RefSeq gene.";
	echo "" >> RefSeq_distribution.txt;
	echo --$input-- >> RefSeq_distribution.txt;

#Counting reads mapping to RefSeq whole genes.
	echo Start aligning $input to RefSeq whole genes;
	intersectBed -s -u -wa -a $inputdir/$input.bed12 -b $refdir/RefSeq/hg38_RefSeq_wholegene.bed > wholegene.temp;
	wholegene=`cat wholegene.temp|wc -l`;
	echo -e Whole_gene '\t' $wholegene >> RefSeq_distribution.txt;
	echo Reads mapping to RefSeq whole genes = $wholegene;
	wait;

#Counting reads mapping to RefSeq introns.
	echo Start aligning $input to RefSeq introns;
	intersectBed -s -u -wa -a $inputdir/$input.bed12 -b $refdir/RefSeq/hg38_RefSeq_introns.bed > introns.temp;
	introns=`cat introns.temp|wc -l`;
	echo -e Introns '\t' $introns >> RefSeq_distribution.txt;
	echo Reads mapping to RefSeq introns = $introns;
	wait;

#Counting reads mapping to RefSeq exons.
	echo Start aligning $input to RefSeq exons;
	intersectBed -s -u -wa -a $inputdir/$input.bed12 -b $refdir/RefSeq/hg38_RefSeq_exons.bed > exons.temp;
	exons=`cat exons.temp|wc -l`;
	echo -e Exons '\t' $exons >> RefSeq_distribution.txt;
	echo Reads mapping to RefSeq exons = $exons;
	wait;

#Counting reads mapping to RefSeq 5'UTR exons.
	echo Start aligning $input to RefSeq 5\'UTR exons;
	intersectBed -s -u -wa -a $inputdir/$input.bed12 -b $refdir/RefSeq/hg38_RefSeq_5-UTRexons.bed > 5UTRexons.temp;
	UTR_5_exons=`cat 5UTRexons.temp|wc -l`;
	echo -e 5\'UTR_Exons '\t' $UTR_5_exons >> RefSeq_distribution.txt;
	echo Reads mapping to RefSeq 5\'UTR exons = $UTR_5_exons;
	wait;

#Counting reads mapping to RefSeq coding exons.
	echo Start aligning $input to RefSeq coding exons;
	intersectBed -s -u -wa -a $inputdir/$input.bed12 -b $refdir/RefSeq/hg38_RefSeq_codingexons.bed > codingexons.temp;
	codingexons=`cat codingexons.temp|wc -l`;
	echo -e Coding_exons '\t' $codingexons >> RefSeq_distribution.txt;
	echo Reads mapping to RefSeq Coding exons = $codingexons;
	wait;

#Counting reads mapping to RefSeq 3'UTR exons.
	echo Start aligning $input to RefSeq 3\'UTR exons;
	intersectBed -s -u -wa -a $inputdir/$input.bed12 -b $refdir/RefSeq/hg38_RefSeq_3-UTRexons.bed > 3UTRexons.temp;
	UTR_3_exons=`cat 3UTRexons.temp|wc -l`;
	echo -e 3\'UTR_Exons '\t' $UTR_3_exons >> RefSeq_distribution.txt;
	echo Reads mapping to RefSeq 3\'UTR exons = $UTR_3_exons;
	wait;

done

#Remove intermediate
echo $input analysis is finished, remove intermediate
rm *temp


