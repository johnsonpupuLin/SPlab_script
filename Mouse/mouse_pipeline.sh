#!/bin/bash

#Program
#	SP lab small RNA sequencing mapping pipeline.
#Usage
#	Link this script to where the small RNA-seq fastq file locating.
#	It's for human dataset now. Modify all the reference and database directories if using on other species.
#	Execute as below:
#		./human_pipeline.sh <filename without ".fastq">
#History
#2018/10/15, HHL	
#	Add variation, please input file name (no .fastq) as var1
#2019/01/23, HHL	
#	Add description.
#2019/01/25, HHL
#	Modified for mouse.

seq_dir=$PWD
seq_name=$1
#adapt_seq="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC"
ref_dir="/data/usrhome/LabSPLin/splin02/Genomes/UCSC/mm10"
rRNAtRNA_db_dir="/data/usrhome/LabSPLin/splin02/pipeline/database/mouse"
miRNA_db_dir="/data/usrhome/LabSPLin/splin02/Genomes/database/mirbase/22/mmu"
sh_dir="/data/usrhome/LabSPLin/splin02/pipeline/sh"
#gff_dir="/data/usrhome/LabSPLin/splin02/pipeline/gff/human"

cd $seq_dir

##FastQC
#module load FastQC
#fastqc $seq_name".fastq"

##cutadapt
#module load cutadapt
#cutadapt -a $adapt_seq -m 15 -M 45 -O 4 -o $seq_name"_cutadapt.fastq" ./$seq_name".fastq"

##FastQC_cutadapt
#module load FastQC
#fastqc $seq_name"_cutadapt.fastq"


##mapping to genome/transcriptome --mapping1
#module load TopHat
tophat2 -o $seq_name"_tophat_genome_transcriptom" -p 16 -g 1 -N 0 \ 
--read-gap-length 0 --read-edit-dist 0 --read-realign-edit-dist 0 --max-insertion-length 0 --max-deletion-length 0 \
--transcriptome-index=$ref_dir/RefSeq/bt2-index/mm10_RefSeq -G $ref_dir/RefSeq/mm10_RefSeq.gtf $ref_dir/bt2-index/mm10 $seq_name".fastq"


cd $seq_name"_tophat_genome_transcriptom"

##bam2sam
#module load SAMtools
samtools view accepted_hits.bam > mapping1.sam
wait

##U6_count
#module load HTSeq
#htseq-count -m intersection-strict -t exon -s no -i ID -o repname_U6.sam mapping1.sam $gff_dir/repname.gff > htseq-count_repname_U6.log 2> htseq-count_repname_U6.err;
#wait

##bam to fastq
#module load bam2fastq
bamToFastq -i accepted_hits.bam -fq $seq_name"_accepted_hit_mapping1.fastq"

##rRNAtRNA_filter
#module load Bowtie
bowtie --norc -S -v 0 -p 16 --al rRNAtRNA.fastq --un rRNAtRNA_removed.fastq $rRNAtRNA_db_dir/Rfam/rRNAtRNA_full $seq_name"_accepted_hit_mapping1.fastq" rRNAtRNA.sam > rRNAtRNA.log 2> rRNAtRNA.err

##miRDeep2 prediction
#module load miRDeep2
cat rRNAtRNA_removed.fastq | awk '{if (NR % 4 == 1){sub(/@/, ">");print;}; if(NR % 4 == 2) print}' > rRNAtRNA_removed.fa
mapper.pl rRNAtRNA_removed.fa -n -c -j -l 18 -m -p $ref_dir/bt-index/mm10 -s reads_collapsed.fa -t reads_collapsed_vs_genome.arf -v -o 16
miRDeep2.pl reads_collapsed.fa $ref_dir/mm10.fa reads_collapsed_vs_genome.arf $miRNA_db_dir/mature_mmu.fa $miRNA_db_dir/mature_rno_hsa.fa $miRNA_db_dir/hairpin_mmu.fa -t mmu 2> report.log
#use -g -1 in case default fail to work

##miRDeepPredictedPrecursors
cat mirdeep_runs/run_*/output.mrd | egrep "^pri_seq" | cut -c26- | tr "u" "t" | tr "atcg" "ATCG" | perl $sh_dir/seq2fasta.pl > miRDeepPredictedPrecursors.fasta


##filter miRNA
#module load Bowtie
cat $miRNA_db_dir/hairpin_mmu.fa > miRBase_pre.fasta
bowtie-build miRBase_pre.fasta miRBase_pre &
bowtie-build miRDeepPredictedPrecursors.fasta miRDeepPredictedPrecursors &
wait
bowtie -norc -S -v 0 -p 16 --al miRBase.fastq --un miRBase_removed.fastq miRBase_pre rRNAtRNA_removed.fastq miRBase.sam > miRBase_removed.log 2> miRBase_removed.err
bowtie -norc -S -v 0 -p 16 --al miRDeep2.fastq --un miRDeep2_removed.fastq miRDeepPredictedPrecursors miRBase_removed.fastq miRDeep2.sam > miRDeep2_removed.log 2> miRDeep2_removed.err


##filterByLen
cat miRDeep2_removed.fastq | perl $sh_dir/filterFastqByLen.pl 24 34 > filteredByLen.fastq 2> filteredByLen.err


##map to genome/transcriptome --mapping2
#module load TopHat
tophat2 -o tophat_mapping2 -p 8 -g 1 -N 0 \
--read-gap-length 0 --read-edit-dist 0 --read-realign-edit-dist 0 --max-insertion-length 0 --max-deletion-length 0 \
--transcriptome-index=$ref_dir/RefSeq/bt2-index/mm10_RefSeq -G $ref_dir/RefSeq/mm10_RefSeq.gtf $ref_dir/bt2-index/mm10 filteredByLen.fastq

cd ./tophat_mapping2


##bam2sam
#module load SAMtools
samtools view accepted_hits.bam > filteredByLen.sam


