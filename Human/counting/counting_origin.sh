#!/bin/bash

#Program
#	Calculating the proportion of piRNA originated from different genome region
#Usage
#       Modify variation "inputdir" to where piRNA candidate bed file locating. (Generated via mapping SP lab pipeline)
#	Modify variation "desdir" to where you want to store output file.
#       Modify variation "refdir" to where the UCSC annotated gene files (including RefSeq and RepeatMasker) locating.
#       Modify for loop input to incorporate all file you want.
#       Execute.
#History
#2018/11/07, HHL
#	Modified for human fetal gonad dataset.
#2019/01/22, HHL
#	Add description.


inputdir=/data/usrhome/LabSPLin/splin01/human_gonad_piRNA/bed/20130502/
desdir=
refdir=/data/usrhome/LabSPLin/splin02/Genomes/UCSC/hg38/

for input in $(ls $inputdir/|grep '.bed12'|sed 's/.bed12//g')
do
	echo Input $input small RNA-seq data;
	echo "" >> $desdir/origin_proportion;
	echo --$input-- >> $desdir/origin_proportion;

#Counting total reads
	echo Start counting total reads;
	total=`cat $inputdir/$input.bed12|wc -l`;
	echo Total_reads $total >> $desdir/origin_proportion;
	echo Total reads = $total;
	wait;

#Counting reads originated from genic region
	echo Start aligning $input to genic region and counting reads;
	intersectBed -s -u -wa -a $inputdir/$input.bed12 -b $refdir/RefSeq/hg38_RefSeq_annotated.bed > $desdir/refseq.temp;
	genic=`cat refseq.temp|wc -l`;
	echo Genic $genic >> $desdir/origin_proportion;
	echo Reads originated from genic region = $genic;
	wait;

#Counting reads originated from repeat region
	echo Start aligning $input to repeat region and counting reads;
	intersectBed -s -u -wa -a $inputdir/$input.bed12 -b $refdir/RepeatMasker/hg38_RepeatMasker_annotated.bed > $desdir/repeatmasker.temp;
	repeat=`cat repeatmasker.temp|wc -l`;
	echo Repeat $repeat >> $desdir/origin_proportion;
	echo Reads originated from repeat region = $repeat;
	wait;

#Counting reads originated from genic&repeat fusion region 
	echo Start aligning $input to genic"&"repeat fusion region and counting reads;
	intersectBed -s -u -wa -a $desdir/refseq.temp -b $refdir/RepeatMasker/hg38_RepeatMasker_annotated.bed > fusion.temp;
	fusion=`cat fusion.temp|wc -l`;
	echo Genic"&"Repeat $fusion >> $desdir/origin_proportion;
	echo Reads originated from repeat region = $fusion;
	wait;

#Counting reads originated from intergenic region
	echo Start counting reads from intergenic region;
	echo Intergenic $(($total-$genic-$repeat+$fusion)) >> $desdir/origin_proportion;
	echo Reads originated from intergenic region = $(($total-$genic-$repeat+$fusion));
	wait;

	echo $input is finished,remove temp;
	rm *.temp;
	echo "" >> $desdir/origin_proportion;
	wait

done

