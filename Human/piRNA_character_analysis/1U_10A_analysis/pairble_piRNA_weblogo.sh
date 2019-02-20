#!/bin/bash

#Program
#	Sort 10nt-overlapping reads generated via overlap_uniq_reads.sh, then make weblogo.
#Usage
#	Modify variation "inputDir" to where pairable piRNA bed file locating.
#       Modify variation "desDir" to where you want to store output file.
#	Modify for loop input to incorporate all file you want.
#	Execute.
#History
#2017/01/19, HHL
#	First version
#2018/11/10, HHL
#	Second version, modified and designed for human fetal gonad piRNA candidate analysis.
#2019/01/21, HHL
#	Add description.


inputDir="$HOME/human_gonad_piRNA/data/20181110/piRNA_overlap_analysis/20130502"
desDir=
echo $inputDir;

for INPUT in $(ls $inputDir|grep '_uniq_pairable.bed'|sed 's/\_uniq_pairable.bed//g');

do
        echo sorting $INPUT;
	awk '$13==10' $inputDir/$INPUT"_uniq_pairable.bed"|awk '$8<=$2&&$2<=$9&&$9<=$3' > $desDir/$INPUT"_10ntoverlapping.bed";
	wait;
	echo $INPUT 10nt overlapped pairs are output;
	awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4, $5, $6}' $INPUT"_10ntoverlapping.bed"|uniq >> $desDir/temp;
	awk 'BEGIN{OFS="\t"}{print $7, $8, $9, $10, $11, $12}' $INPUT"_10ntoverlapping.bed"|uniq >> $desDir/temp;
	awk '{for (i=1;i<=$5;i++){print}}' temp > $desDir/$INPUT"_10ntoverlapping_whole.bed"
	wait;
	echo $INPUT 10nt overlapped reads were isolated
	awk '$3-$2<35' $INPUT"_10ntoverlapping_whole.bed" |fastaFromBed -name -s -fi $HOME/genomes/UCSC/hg38/hg38.fa -bed - -fo $desDir/$INPUT.fasta
	wait;
	echo $INPUT.fasta was generated;
	awk -F "\t" 'NR%2==0{$1=toupper(substr($1,1,20))}{print}' $desDir/$INPUT.fasta| sed 's/T/U/g'| 
		weblogo -D fasta --format PNG -U bits -S 1 --ticmarks 0.5 -o $desDir/$INPUT.pingpong.weblogo.png -A rna -W 30 -c classic;
	wait;
	echo weblogo complete, delete intermediate
	#rm temp $INPUT"_10ntoverlapping.bed" $INPUT"_10ntoverlapping_whole.bed" $INPUT.fasta
	rm $desDir/temp
done
