#!/bin/bash

#Program
#	Using bowtie and HTseq to charaterize human small RNA. It a method derived from the 2018 paper from Diana J. Laird. (https://doi.org/10.1101/364315)
#Usage
#       Modify variation "input_path" to where small RNA-seq raw datasets locating.
#	Modify variation "des_path" to where you want to output your result.
#       Modify for loop input to incorporate all file you want.
#       Execute.
#History
#	2019/1/14, HHL
#	Version 1
#	2019/1/23
#	Add description.


input_path="/data/usrhome/LabSPLin/splin01/human_gonad_piRNA/raw_data/20130502"
des_path="/data/usrhome/LabSPLin/splin01/human_gonad_piRNA/data/20190114_HTseq/20130502/"
ref_path="/data/usrhome/LabSPLin/splin02/genomes"
echo Input $(ls $input_path|grep "_cutadapt.fastq"|sed 's/_cutadapt.fastq//g')


for i in $(ls $input_path|grep "_cutadapt.fastq"|sed 's/_cutadapt.fastq//g')
do

echo Input and Generate $i directory
mkdir $des_path/$i
cd $des_path/$i

#Using Bowtie to map small RNA-seq data to Repbase database.
echo Start mapping $i to Repbase via Bowtie
bowtie --sam -p 8 -v 1 -M 1 --best --strata --al Rep_aligned.fastq --un Rep_unaligned.fastq \
$ref_path/database/RepBase/bt-index/humrep $input_path/$i"_cutadapt.fastq" > Rep_aligned.sam 2> Rep_aligned_log.txt

#Using Bowtie to unaligned data to hg38 reference.
echo Start mapping $i Repbase-unaligned reads to hg38 fasta via Bowtie.
bowtie --sam -p 8 -v 1 -M 1 --best --strata --al hg38_aligned.fastq --un hg38_unaligned.fastq \
$ref_path/UCSC/hg38/bt-index/hg38 Rep_unaligned.fastq > hg38_aligned.sam 2> hg38_aligned_log.txt

#Using HTseq to count reads in feature of hg38-aligned data.
echo Start counting $i hg38-aligned reads in feature via HTseq.
htseq-count -f sam -t exon -s yes -i gene_type -m intersection-nonempty --nonunique all -a 0 -o $i"_HTseq.sam" \
hg38_aligned.sam $ref_path/database/GENCODE/gencode.v29.whole.gtf > $i"_HTseq.txt" 2> $i"_HTseq_log.txt"

echo $i mapping is done
cd $des_path

done

