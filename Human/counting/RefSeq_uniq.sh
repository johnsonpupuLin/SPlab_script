#!/bin/bash

#Program
#	To print the isomer which is the longest among all isomers
#Usage
#	Execute directly.
#History
#2018/11/16, HHL
#	First version
#2019/01/22, HHL
#	Add description


#awk '{print $5}' /data/usrhome/LabSPLin/splin01/genomes/UCSC/hg38/RefSeq/hg38_RefSeq_annotated.bed |sort|uniq > gene_name.txt

#It is unclear why the grep function could not precisely print what I want.
#For example, grep '<\AAA\>' will give me AAA and AAA-AS1.
#To solve this problem, isolate all -ASn (antisense gene) from RefSeq bed file first.
#	grep '\-[a-zA-Z0-9]' gene_name.txt > antisense_gene.txt
#	grep -v '\-[a-zA-Z0-9]' gene_name.txt > sense_gene.txt
#	grep '\-[a-zA-Z0-9]' /data/usrhome/LabSPLin/splin01/genomes/UCSC/hg38/RefSeq/hg38_RefSeq_annotated.bed > antisense_gene.bed
#	grep -v '\-[a-zA-Z0-9]' /data/usrhome/LabSPLin/splin01/genomes/UCSC/hg38/RefSeq/hg38_RefSeq_annotated.bed > sense_gene.bed


#Start greping sense gene
echo Start greping sense gene

for i in $(cat sense_gene.txt)
do

	echo Now greping $i

	grep -w -i "$i" sense_gene.bed > temp
	M=`awk 'BEGIN {max=0} {if (($3-$2)+0 > max+0) max=($3-$2)} END{print max}' "temp"`
	awk '{if (($3-$2)==M){print}}' M=$M temp|uniq -f 4 >> sense_temp.bed

done


#Start greping anti-sense gene
echo Start greping anti-sense gene

for i in $(cat antisense_gene.txt)
do

	echo Now greping $i

	grep -w -i "$i" antisense_gene.bed > temp
	M=`awk 'BEGIN {max=0} {if (($3-$2)+0 > max+0) max=($3-$2)} END{print max}' "temp"`
	awk '{if (($3-$2)==M){print}}' M=$M temp|uniq -f 4 >> antisense_temp.bed

done

#Combine two sense RefSeq gene.
#As some longest RefSeq would be printed multiple times, uniq after sorting.
cat sense_temp.bed antisense_temp.bed|sort -k 5|uniq > hg38_RefSeq_uniq.bed

#Remove intermediate.
rm temp 
#rm sense_temp.bed antisense_temp.bed


