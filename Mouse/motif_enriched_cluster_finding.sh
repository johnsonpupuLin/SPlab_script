### link needed file
mkdir fasta
mkdir motif
mkdir HOMER


###### run motif find
homer2 find -i fasta/piRClust_upregulated.fasta -strand + -o HOMER/HOMER.piRClust_upregulated.zbtb12 -m motif/zbtb12.motif -offset 1 -p 4;
homer2 find -i fasta/piRClust_upregulated.fasta -strand + -o HOMER/HOMER.piRClust_upregulated.stat3 -m motif/stat3.motif -offset 1 -p 4;
homer2 find -i fasta/piRClust_upregulated.fasta -strand + -o HOMER/HOMER.piRClust_upregulated.foxp1 -m motif/foxp1.motif -offset 1 -p 4;
homer2 find -i fasta/piRClust_downregulated.fasta -strand + -o HOMER/HOMER.piRClust_downregulated.otx2 -m motif/otx2.motif -offset 1 -p 4;
homer2 find -i fasta/piRClust_downregulated.fasta -strand + -o HOMER/HOMER.piRClust_downregulated.klf6 -m motif/klf6.motif -offset 1 -p 4;
homer2 find -i fasta/piRClust_downregulated.fasta -strand + -o HOMER/HOMER.piRClust_downregulated.sp5 -m motif/sp5.motif -offset 1 -p 4;

######### make piRNA targeting fasta
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

#for i in $(ls HOMER/); do awk 'length(a[$1])==0{print $1;a[$1]=1}' HOMER/$i | grep -w -f - ../20160616_combine_clusters/all.cluster.fpkm_adjust.bed > temp.bed; cat ~/chicken_gonad_oxidation_piRNA/bed/*.uniq.bed12| intersectBed -wa -s -a - -b temp.bed > temp.intersectBed; awk 'BEGIN{FS=OFS="\t"}{for(i=0;i<$13;i++){t=$4;$4=t"_rep"i;print;$4=t}}' temp.intersectBed |fastaFromBed -s -name -fi ~/genomes/UCSC/galGal4/galGal4.fa -bed - -fo motif.piRC_class.piRNA/$i.sense.fasta; awk 'BEGIN{FS=OFS="\t"}{if($6=="+"){$6="-"}else{$6="+"} ;for(i=0;i<$13;i++){t=$4;$4=t"_rep"i;print;$4=t}}' temp.intersectBed |fastaFromBed -s -name -fi ~/genomes/UCSC/galGal4/galGal4.fa -bed - -fo motif.piRC_class.piRNA/$i.asense.fasta; wait; done; rm temp.*;

