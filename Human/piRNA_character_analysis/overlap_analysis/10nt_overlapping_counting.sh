#!/bin/bash

#Program
#       Counting the amounts of 10nt-overlapping pairable piRNA.
#Usage
#	Link this script to where you want to store your output data.
#       Modify variation "inputdir" to where uniq piRNA bed file locating. (Generated via uniq_and_strand.sh)
#       Modify for loop input to incorporate all file you want.
#       Execute.
#History
#2019/01/22, HHL
#       Add description


inputdir=$HOME/human_gonad_piRNA/bed/20130502

#Overlap uniq piRNA candidate
for input in $(ls $inputdir|grep 'bed12'|sed 's/\.bed12//g')
do
	intersectBed -wo -a $inputdir/$input.uniq.plus.bed -b $inputdir/$input.uniq.minus.bed > $input"_uniq_pairable.bed"
done


#Counting the total number of piRNA candidates that are able to overlap 10nt with others
for input in $(ls ./|grep '_uniq_pairable.bed'|sed 's/_uniq_pairable.bed//g');
do
a=`cat $inputdir/$input.bed12|wc -l`
	echo $input total reads = $a
	awk '{if ($13==10&&$8<=$2&&$2<=$9&&$9<=$3){print}}' $input"_uniq_pairable.bed" |
		awk 'BEGIN{print "+","\t","-","\t","total"}{a+=$5;b+=$11}END{print a,"\t",b,"\t",a+b}' >> 10nt_overlapping_count;
	echo $input total reads = $a >> 10nt_overlapping_count
	echo "" >> 10nt_overlapping_count
done


