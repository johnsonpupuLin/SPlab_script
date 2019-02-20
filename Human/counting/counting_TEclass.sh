#!/bin/bash

#Program
#	Analyze piRNAs which originate from TE belong to which TE class.
#Usage
#       Link this script to where you want to store output file.
#       Modify variation "inputdir" to where piRNA candidate bed file locating. (Generated via mapping SP lab pipeline)
#	Modify variation "refdir" to where the RepeatMasker file locating.
#       Modify for loop input to incorporate all file you want.
#       Execute.
#History
#2018/11/07, HHL
#	Modified for human fetal gonad
#2019/01/23, HHL
#	Add description.

inputdir=/data/usrhome/LabSPLin/splin01/human_gonad_piRNA/bed/20130502/
refdir=/data/usrhome/LabSPLin/splin02/Genomes/UCSC/hg38/RepeatMasker

for input in $(ls $inputdir|grep '.bed12'|sed 's/.bed12//g')
do

echo Start aligning $input to repeat region and counting reads;
intersectBed -s -wa -a $refdir/hg38_RepeatMasker_annotated.bed -b $inputdir/$input.bed12 > repeatmasker.temp;

echo Start classifing TEs;
cat repeatmasker.temp|awk '{split ($4,a,"|");print a[3]}'|sort -V|uniq -c|sort -r -n > $input.repeatmasker.Class;

echo $input is finished;
done

rm repeatmasker.temp
