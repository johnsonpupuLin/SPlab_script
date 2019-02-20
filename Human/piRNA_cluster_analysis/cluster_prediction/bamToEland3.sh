#!/bin/bash

#Program
#	Transfer bam file into eland3 file.
#Usage
#	Execute as below,
#	./bamToEland3.sh <bam file name> <genome file name>
#	Example,
#	./bamToEland3.sh M12_1 hg38
#	The output eland3 file will locate in bam file directory.
#History
#2019/01/28, HHL
#	This Script is from KW. Description is added by HHL.




#piRNAs to the genome. proTRAC requires a map file that contains the (sorted)
#coordinates of mapped piRNAs in the following format (ELAND3 tab-delimited table):

#Chr1	7400394	TTGCTACGTCAGATCGTGCGGGTAA	12	TTGCTACGTCAGATCGTGCGGGTCC	2	+

#where each column referres to:
#chromosome, coordinate, target sequence, number of reads, query sequence, number of mismatch, strand.

if [ $#	-ne 2 ]; then
	echo '$1 = bam file without suffix (<$1>.bam)'
	echo '$2 = genome file without suffix (<$2>.fa)'
	exit 1
fi

#Mapping data, from mapped data
#samtools view $1.bam | awk 'BEGIN{OFS="\t"}$2==0{print $1,$3,$4,$10,$6,"+"}$2==16{print $1,$3,$4,$10,$6,"-"}' - > $1.mapData

#Target sequence
bamToBed -bed12 -i $1.bam > $1.mapData
fastaFromBed -s -name -tab -split -fi $2.fa -bed $1.mapData -fo $1.fasta

wait;

#Join Table (assume no mismatch)
join -1 4 -2 1 $1.mapData $1.fasta | awk 'BEGIN{OFS="\t"}{print $2,$3,toupper($13),0,$6}' - |sort -V|uniq -c|awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$4,$5,$6}' > $1."eland3" ;

wait;
samtools view -c $1.bam | awk -v nm=$1 'END{print "\t"nm".bam:\n\t\ttotal lines(reads): "$1}' -
awk -v nm=$1 'END{print "\t"nm".mapData:\n\t\ttotal lines(reads): "NR}' $1.mapData
awk -v nm=$1 'END{print "\t"nm".fasta:\n\t\ttotal lines(reads): "NR}' $1.fasta
awk -v nm=$1 'BEGIN{s=0}{s=s+$4}END{print "\t"nm".eland3:\n\t\ttotal lines: "NR"\n\t\ttotal reads: "s}' $1.eland3

wait;
rm $1.mapData
rm $1.fasta
