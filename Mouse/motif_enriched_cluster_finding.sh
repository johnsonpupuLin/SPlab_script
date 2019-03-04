#!/bin/bash

#Program
#	After recognizing highly represented motifs among cluster sets, the cluster possessing the specific motifs are identified.
#	Next, the piRNAs generated from these clusters are isolated, then their potential targeting sites are identified.
#Usage
#	It is a flow, but not a script of analysis. Please refer to the command and modify by yourself.
#History
#2019/03/04, HHL
#	Modified. Description added.


#Make directories and link needed file.
mkdir fasta
mkdir motif
mkdir HOMER


#Apply Homer script to identify the clusters possessing specific motifs.
homer2 find -i fasta/piRClust_upregulated.fasta -strand + -o HOMER/HOMER.piRClust_upregulated.zbtb12 -m motif/zbtb12.motif -offset 1 -p 4;
homer2 find -i fasta/piRClust_upregulated.fasta -strand + -o HOMER/HOMER.piRClust_upregulated.stat3 -m motif/stat3.motif -offset 1 -p 4;
homer2 find -i fasta/piRClust_upregulated.fasta -strand + -o HOMER/HOMER.piRClust_upregulated.foxp1 -m motif/foxp1.motif -offset 1 -p 4;
homer2 find -i fasta/piRClust_downregulated.fasta -strand + -o HOMER/HOMER.piRClust_downregulated.otx2 -m motif/otx2.motif -offset 1 -p 4;
homer2 find -i fasta/piRClust_downregulated.fasta -strand + -o HOMER/HOMER.piRClust_downregulated.klf6 -m motif/klf6.motif -offset 1 -p 4;
homer2 find -i fasta/piRClust_downregulated.fasta -strand + -o HOMER/HOMER.piRClust_downregulated.sp5 -m motif/sp5.motif -offset 1 -p 4;


#Identify the piRNAs generated from those clusters, then investigate their potential targeting sites.
mkdir HOMER.piRClust.beds
mkdir motif.piRC_class.piRNA

for i in $(ls HOMER/|grep 'upregulated');
do 
awk 'BEGIN{FS=OFS="\t"}FNR==NR&&length(a[$1])==0{a[$1]=1;next;}a[$4]==1{print}' HOMER/$i $HOME/mouse_gonad_piRNA/HHL/data/20170313/piRNA_cluster_RNA_RPKM/condensed_AB_merged_cluster.bed  > HOMER.piRClust.beds/$i.bed; 
intersectBed -split -wa -s -a bed/A_DL0.bed12 -b HOMER.piRClust.beds/$i.bed > temp.intersectBed; 
wait;
awk '{if(($3-$2)<=34){print}}' temp.intersectBed > temp2.intersectBed;
wait;
fastaFromBed -split -s -name -fi $HOME/genomes/UCSC/mm10/mm10.fa -bed temp2.intersectBed -fo motif.piRC_class.piRNA/$i.sense.fasta; 
wait;
awk 'BEGIN{FS=OFS="\t"}{if($6=="+"){$6="-"}else{$6="+"}{print}}' temp2.intersectBed |fastaFromBed -split -s -name -fi ~/genomes/UCSC/mm10/mm10.fa -bed - -fo motif.piRC_class.piRNA/$i.asense.fasta; 
wait; 
done 

for i in $(ls HOMER/|grep 'downregulated');
do
awk 'BEGIN{FS=OFS="\t"}FNR==NR&&length(a[$1])==0{a[$1]=1;next;}a[$4]==1{print}' HOMER/$i $HOME/mouse_gonad_piRNA/HHL/data/20170313/piRNA_cluster_RNA_RPKM/condensed_AB_merged_cluster.bed  > HOMER.piRClust.beds/$i.bed;
intersectBed -split -wa -s -a bed/A_WL0.bed12 -b HOMER.piRClust.beds/$i.bed > temp.intersectBed;
wait;
awk '{if(($3-$2)<=34){print}}' temp.intersectBed > temp2.intersectBed;
wait;
fastaFromBed -split -s -name -fi $HOME/genomes/UCSC/mm10/mm10.fa -bed temp2.intersectBed -fo motif.piRC_class.piRNA/$i.sense.fasta;
wait;
awk 'BEGIN{FS=OFS="\t"}{if($6=="+"){$6="-"}else{$6="+"}{print}}' temp2.intersectBed |fastaFromBed -split -s -name -fi ~/genomes/UCSC/mm10/mm10.fa -bed - -fo motif.piRC_class.piRNA/$i.asense.fasta;
wait;
done


rm temp.*

