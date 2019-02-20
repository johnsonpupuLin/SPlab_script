
#!/bin/bash


#Program:
#	For automatically running QuantTE by simply input the sample name
#History:
#	2018/10/29 HHL made


for i in $@

do
	echo "Start inputing "$@

	cd /data/usrhome/LabSPLin/splin01/Axolotl/QuantTE_mapping/annotated_coding_gene_mapping/$i;
	sh loop.sh;
	wait;
	
	cd /data/usrhome/LabSPLin/splin01/Axolotl/QuantTE_mapping/annotated_coding_gene_mapping;
done
