
#!/bin/bash

#Program:
#	Combine merge_abundance.tsv file of split QuantTE result, followed by sample_loop.sh
#History:
#2018/10/29 HHL
#	Version 1


for i in "$@"

do
	echo "Combining "$i" split tpm"

	cd $i
	awk '(NR == 1){print}' kallisto/$i"_1"/merge_abundance.tsv > $i"_tpm.tsv"
	awk '(NR > 1){print}' kallisto/$i"_"[1-7]/merge_abundance.tsv|grep 'AMEXTC' >> $i"_tpm.tsv"
	cd /data/usrhome/LabSPLin/splin01/Axolotl/QuantTE_mapping/annotated_coding_gene_mapping

	wait;
done
