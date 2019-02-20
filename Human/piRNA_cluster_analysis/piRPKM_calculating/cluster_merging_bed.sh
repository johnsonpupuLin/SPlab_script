#!/bin/bash

#Program
#	Merge clusters generated from different datasets into a integrated cluster file.
#Usage
#	Copy this script to where the clusters you want to merge locating.
#	Execute as below:
#	./cluster_merging_bed.sh <cluster 1.bed> <cluster 2.bed> ...
#History
#2019/01/29, HHL
#	This script is made by KW. HHL add description.


#take bed as input
echo $@
cat $@ | sort -k6,6 -k1V,3 |
	awk -v bw=0 'BEGIN{OFS="\t";ch="";st="";ed="";dr=""}
	     NR==1{ch=$1;st=$2;ed=$3;dr=$6}
             NR>1{if(ch==$1&&dr==$6&&((ed>=$3&&$3>=st)||($3>=ed&&ed>=$2)||(ed+bw>=$2))) {
			ch=$1;dr=$6;st=(st<$2 ? st:$2);ed=(ed>$3 ? ed:$3);
		  } else {
			print ch,st,ed,dr;
			ch=$1;st=$2;ed=$3;dr=$6}
		  }
	     END{print ch,st,ed,dr}' | 
	sort -V | 
	awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"merged_cluster_"NR,0,$4}' |
	tee merged.cluster.bed |
        awk -v bidirTH=0.1 '
                BEGIN{FS=OFS="\t"}
                NR==1{ch=$1;st=$2;ed=$3;v1=$2;v2=$3;n=1;sd=$6;next}
                $1==ch&&st<=$2&&$2<=ed{
                        v1=$2;v2=ed;bidirKey=(((ed-$2)/(ed-st))<=bidirTH&&((ed-$2)/($3-$2))<=bidirTH);
                        st=(st<$2?st:$2);ed=(ed>$3?ed:$3);
                        if((sd=="+-"&&$6=="-")||(sd=="-+"&&$6=="+")){next};
                        if(sd=="="||sd==$6){v1=st;v2=ed;next};
                        if(length(sd)==1&&bidirKey){sd=sd""$6}else{v1=st;v2=ed;sd="="};next;}
                {print ch,st,ed,"merged_cluster_"n,v1"-"v2,sd;ch=$1;st=$2;ed=$3;v1=$2;v2=$3;n=n+1;sd=$6}
                END{print ch,st,ed,"merged_cluster_"n,v1"-"v2,sd;}' > merged.cluster.directioned.bed;

