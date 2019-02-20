#!/bin/bash

#Program
#       Using intersectBed to overlap uniq-standed piRNAs generated via uniq_and_strand.sh to investigate pairable piRNAs. Outputing pairable piRNA as bed file.
#Usage
#       Link this script to where you want to store your output data.
#	Modify variation "inputdir" to where uniq-stranded piRNA bed file locating and execute directly.
#History
#       2019/01/21, HHL
#	Modified and Added description.


inputdir=$HOME/human_gonad_piRNA/bed/20130502

for input in $(ls $inputdir|grep 'bed12'|sed 's/\.bed12//g')
do

	intersectBed -wo -a $inputdir/$input.uniq.plus.bed -b $inputdir/$input.uniq.minus.bed > $input"_uniq_pairable.bed"

done
