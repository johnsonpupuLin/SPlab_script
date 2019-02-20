#!/bin/bash

#Program:
#	Run QuantTE to pseudoalign sample fasta to 7 split annotated coding gene fasta file
#History
#2018/10/23 HHL
#	Version 1, for repeat running


sample="O-ST-2"

echo "Start QuantTE of sample" $sample;

for i in 2 3 4 5 6 7

do

#Run kallsito
	echo "Start running kallisto and pseudoaligning to split fasta file #"$i;
	python3 $HOME/programs/QuantTE-master/main.py --kallisto --input json_file/"input_"$i".json";
	echo "#"$i" kallisto is finished"
	wait;

#Change directory name
	cd kallisto;
	mv $sample $sample"_"$i
	cd $HOME/Axolotl/QuantTE_mapping/annotated_coding_gene_mapping/$sample
	wait;

done
