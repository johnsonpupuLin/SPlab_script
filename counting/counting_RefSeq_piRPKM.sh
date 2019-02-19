#!/bin/bash

#Program
#	Calculating the piRNA expression (piRPKM) of RefSeq genes. RefSeq gene sets are modified via "RefSeq_uniq.sh".
#Usage
#       Link this script to where you want to store output file.
#       Modify variation "inputdir" to where piRNA candidate bed file locating. (Generated via mapping SP lab pipeline)
#	Modify variation "refdir" to where modified RefSeq file locating.
#       Modify for loop input to incorporate all file you want.
#       Execute.
#History
#2018/11/07, HHL
#	Modified for human fetal gonad dataset
#2019/01/23, HHL
#	Add description.


inputdir=/data/usrhome/LabSPLin/splin01/human_gonad_piRNA/bed/20130502
refdir=/data/usrhome/LabSPLin/splin01/human_gonad_piRNA/data/20181115/piRNA_genic_origin/

for input in $(ls $inputdir/|grep '.bed12'|sed 's/.bed12//g')
do

Mappedreads=`samtools view $inputdir/$input/$input"_cutadapt_tophat_genome_transcriptom"/accepted_hits.bam|wc -l`
echo $Mappedreads

echo IntersectBed $input to hg38_RefSeq_uniq.bed;
intersectBed -s -wa -a $refdir/hg38_RefSeq_uniq.bed -b $inputdir/$input.bed12|
awk 'BEGINS{OFS="\t"}{print $5, ($3-$2)}' |sort|uniq -c|sort -k 1 -r -n|
awk 'BEGIN{OFS="\t"}{print $2, $3, $1}' > temp.txt

echo counting RPKM
cat temp.txt|awk '{print $1, $2, $3, ($3/((M/1000000)*($2/1000)))}' M=$Mappedreads |sort -k 4 -r -g|
awk 'BEGIN{OFS="\t";print "RefSeq_gene", "Length", "Counts", "piRPKM"}{print $1, $2, $3 ,$4}' > $input"_RefSeq_origin.txt"


done



