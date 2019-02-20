#!/bin/bash

#Program
#	Extract the result table genrated from proTRAC.
#Usage
#	Execute as below:
#	./extractCluster.sh <result table from proTRAC>
#History
#2019/01/29, HHL
#	This script is made by KW. HHL add description.



# $1 = output prefix 

wait;
# extract cluster information from proTRAC report
grep "^Cluster" results.table |	
	awk -v pref=$1 --posix -F "\t" 'BEGIN{OFS="\t"}{
		sub(" ","_",$1); $1=pref"_"$1; 
		$2=substr($2,length("Location: ")+1);
		$3=substr($3,length("Coordinates: ")+1); sub("-","\t",$3);
		$4=substr($4,length("Size [bp]: ")+1);
		$5=substr($5,length("Hits (absolute): ")+1);
		$6=substr($6,length("Hits (normalized): ")+1);
		$7=substr($7,length("Hits (normalized) per kb: ")+1);
		$8=substr($8,length("Normalized hits with 1T: ")+1);
		$9=substr($9,length("Normalized hits with 10A: ")+1);
		$10=substr($10,length("Normalized hits 24-34 nt: ")+1);
		$11=substr($11,length("Normalized hits on the main strand(s): ")+1);
		$12=substr($12,length("Predicted directionality: ")+1)}
		$13!=""{$13=substr($13,length("Binding sites: ")+1)}{print}' |
	awk -F "\t" 'BEGIN{OFS="\t";
			print "id", "chr", "start", "end", "length", "absHits",
			"normHits","normHitsPerKb","normHits_1T","normHits_10A",
			"normHits_24-34nt","normHits_mainStrand","directionality",
			"Binding_Sites"}
		{print}' > results.curated.table

# make results in bed format,
awk -F "\t" 'BEGIN{OFS="\t"}
	$13~/mono:p/{$13="+"}
	$13~/mono:m/{$13="-"}
	$13~/bi:p/{$13="+-"}
	$13~/bi:m/{$13="-+"}
	$13~/dual/{$13="="}
	NR>1{print $2,$3,$4,$1,$8,$13}' results.curated.table > clusters.bed

# make results in bed format, split bidirectionality
awk -F "\t" 'BEGIN{OFS="\t";out=""}
	$13~/mono:p/&&NR>1{split($13,arr," "); print $2,$3,$4,$1"_"arr[1],$8,"+"}
	$13~/mono:m/&&NR>1{split($13,arr," "); print $2,$3,$4,$1"_"arr[1],$8,"-"}
	$13~/bi:plu/&&NR>1{split($13,arr," "); print $2,$3,substr(arr[6],1,length(arr[6])-1),$1"_"arr[1]"(+)",$8,"+\n"$2,arr[4],$4,$1"_"arr[1]"(-)",$8,"-"}
	$13~/bi:min/&&NR>1{split($13,arr," "); print $2,$3,substr(arr[6],1,length(arr[6])-1),$1"_"arr[1]"(-)",$8,"-\n"$2,arr[4],$4,$1"_"arr[1]"(+)",$8,"+"}
	$13~/dual/&&NR>1{split($13,arr," "); print $2,$3,$4,$1"_"arr[1]"(+)",$8,"+\n"$2,$3,$4,$1"_"arr[1]"(-)",$8,"-"}' results.curated.table > clusters.split.bed

# make results in gtf format,
awk 'BEGIN{OFS="\t"}{print $1,"proTRAC","exon",$2+1,$3,$5,$6,".","gene_id \""$4"\"; transcript_id \""$4"\";"}' clusters.split.bed > clusters.split.gtf
