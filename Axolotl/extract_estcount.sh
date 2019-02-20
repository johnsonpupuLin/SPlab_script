
#!/bin/bash

#Program:
#	Combine abundance.tsv file of split QuantTE result, followed by sample_loop.sh
#History:
#2018/11/05 HHL
#	Version 1


for i in "$@"

do
	echo "Combining "$i" split tpm"

	cd $i
	awk '(NR == 1){print}' kallisto/$i"_1"/abundance.tsv > $i"_estcount.tsv"
	awk '(NR > 1){print}' kallisto/$i"_"[1-7]/abundance.tsv|grep 'AMEXTC' >> $i"_estcount.tsv"
	cd /data/usrhome/LabSPLin/splin01/Axolotl/QuantTE_mapping/annotated_coding_gene_mapping

	wait;
done
