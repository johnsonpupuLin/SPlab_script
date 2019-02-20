#!/bin/bash

#Program
#       Compress same reads and print the number of reads at $5 of bed file. Afterward, Isolate + and - strand reads.
#Usage
#       Modify variation "inputdir" to where bed file locating. (Such as piRNA bed file after pipeline mapping)
#	Modify variation "desDir" to where you want to store output file.
#	Modify for loop input to incorporate all file you want.
#	Execute.
#History
#       2019/01/21, HHL
#	Modified and Added description.



inputDir="$HOME/human_gonad_piRNA/bed/20140214";
desDir="$HOME/human_gonad_piRNA/bed/20140214";

for INPUT in $(ls $inputDir|grep '.bed12'|sed 's/\.bed12//g')
do
	echo "Processing $INPUT ...";
	awk '{print $1"\t"$2"\t"$3"\t"$6}' $inputDir/$INPUT.bed12 | sort -V | uniq -c | awk '$1 = $1' | tr ' ' '\t' | awk '{print $2"\t"$3"\t"$4"\tUniqRead_"NR"\t"$1"\t"$5}' > $desDir/$INPUT.uniq.bed;
	awk '$6=="+"' $desDir/$INPUT.uniq.bed > $desDir/$INPUT.uniq.plus.bed;
	awk '$6=="-"' $desDir/$INPUT.uniq.bed > $desDir/$INPUT.uniq.minus.bed;
done
