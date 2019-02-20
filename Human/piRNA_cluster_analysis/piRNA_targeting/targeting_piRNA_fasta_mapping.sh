#!/bin/bash

#Program
#	Identify the possible targets of piRNAs.	
#Usage
#	Copy this script to where you want to output the result.
#	Modify "fastaDir" to where your antiense piRNA fasta file locating.
#	Modify "refDir" to the species you are tackling.
#	Execute.
#History
#2019/01/31, HHL
#	This script is made by KW. Modify previous version and add description.


# get piR mapping, from motif_assoc piRClust

fastaDir=$HOME/human_gonad_piRNA/data/20190203/cluster_analysis/20140214/piRNA_targeting
refDir=/data/usrhome/LabSPLin/splin02/Genomes

for fname in $(ls $fastaDir |grep ".fasta"|sed 's/.fasta//g'); 

do

#Map the antisense piRNA sequences to human genome.
	mkdir $fname.Tophat2;
	tophat2 -o $fname.Tophat2 -p 16 -g 1 -N 0 --read-gap-length 0 --read-edit-dist 0 --read-realign-edit-dist 0 \
	--max-insertion-length 0 --max-deletion-length 0 --transcriptome-index=/data/usrhome/LabSPLin/splin02/Genomes/UCSC/hg38/RefSeq/bt2-index/hg38_RefSeq \
	-G $refDir/UCSC/hg38/RefSeq/hg38_RefSeq.gtf $refDir/UCSC/hg38/bt2-index/hg38 $fastaDir/$fname.fasta;
	wait;

#Get piRNA targets (unspliced)
	bamToBed -i $fname.Tophat2/accepted_hits.bam |intersectBed -wa -s -a $refDir/UCSC/hg38/RefSeq/hg38_RefSeq_annotated.bed -b -|sort -V|uniq -c|
		awk 'BEGIN{OFS="\t"}{$6=$1;sub($1"\t","");print}' | tee $fname.Tophat2/$fname.transc.bed| 
		awk 'BEGIN{FS=OFS="\t"}length(a[$4])==0{a[$4]=0}{a[$4]+=$5}END{for(i in a){print i,a[i]}}' |sort -V > $fname.Tophat2/$fname.transc.count;
	wait;
#Get piRNA targets (mRNA)
	bamToBed -i $fname.Tophat2/accepted_hits.bam |intersectBed -wa -s -F 1 -split -a $refDir/UCSC/hg38/RefSeq/hg38_RefSeq_annotated.bed -b -|sort -V|uniq -c|
		awk 'BEGIN{OFS="\t"}{$6=$1;sub($1"\t","");print}' | tee $fname.Tophat2/$fname.mRNA.bed| 
		awk 'BEGIN{FS=OFS="\t"}length(a[$4])==0{a[$4]=0}{a[$4]+=$5}END{for(i in a){print i,a[i]}}' |sort -V > $fname.Tophat2/$fname.mRNA.count;
	wait;
#Get piRNA TE targets
	bamToBed -i $fname.Tophat2/accepted_hits.bam |intersectBed -wa -s -a $refDir/UCSC/hg38/RepeatMasker/hg38_RepeatMasker_annotated.bed -b -|sort -V|uniq -c|
		awk 'BEGIN{OFS="\t"}{$6=$1;sub($1"\t","");print}' | tee $fname.Tophat2/$fname.TEs.bed| 
		awk 'BEGIN{FS=OFS="\t"}length(a[$4])==0{a[$4]=0}{a[$4]+=$5}END{for(i in a){print i,a[i]}}' |sort -V > $fname.Tophat2/$fname.TEs.count; 
done

