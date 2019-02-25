


# Instruction

Is's a series of small RNA-seq analysis in SP lab.
I will briefly describe the general flow and the function of scripts in below section.

I might update some scripts of other analysis in the future.


# Small RNA-seq

## General Flow

In general, the analysis is carried out in below order.
There are different scripts in those directories. Moreover, the detailed instruction are written in scripts respectively.

### 1. pipeline

In this step, you will map your sequencing data to genome and identify piRNA candidates.

There are two pipelines. One is Tophat2-based and is constructed by SP lab members. The other one is HTseq-based and is from Diana's lab, mimiced from their paper.

* human_pipeline.sh

	This is SP lab pipeline. 
	The raw fastq data are first mapped to genome via Tophat2. 
	Next, rRNA/tRNA, miRNA and predicted miRNA among mapped reads are eliminated. 
	After that, reads shorter than 24bp and longer than 34bp are discarded.
	In the end, the data are mapped to genome again.

* Diana_HTseq_pipeline.sh

	This is Diana lab pipeline. 
	The raw fastq data are first aligned to Repbase database to identify reads from repeat elements.
	Next, non-aligned reads are mapped to genome via Bowtie.
	Finally, genome-alignmed reads are analyzed by HTseq to identify their characteristics.


### 2. counting

In this step, you will like to clarity the origin of piRNA candidates.

* counting_origin.sh
	
	This script could figur out the reads originating from genic (RefSeq), repeat (RepeatMasker), genic&repeat fusion, or intergenic regions.
	Please not that the "genic&repeat fusion region" is defined as a region containing gene and repeat element simultaneously.
	In this case, this script could not identify whether the read is from gene or repeat.
	It could be revised in the future.

* counting_RefSeq_distribution.sh

	This script is for identifying the exact origin (from intro, exon, etc.) of genic-piRNA candidate.
	However, isoforms are not considered in this script.
	That is, a reads might be aligned to exons in one gene isoform and to intron in another isoform.
	It could be revised in the future.


* RefSeq_uniq.sh

	This script is to integrate the isoform in RefSeq.
	As there are lots of isoform present in UCSC RefSeq datasets, if you apply the origin file, one gene might be counted multiple times due to isoforms.
	Therefore, I wrote a script to pick out the longest isoform as representation.
	It might be rough, but I think it enough for priliminary analysis.

* counting_RefSeq_piRPKM.sh

	This script is to calculate the piRPKM of RefSeq genes. The gene list is modified via "RefSeq_uniq.sh" above.
	
* counting_TEclass.sh

	This script is to classifiy the family and subfamily that piRNA candidates locate.
	The RepeatMasker file has been modified by whole table to get family and subfamily.
	


### 3. piRNA_character_analysis

### 4. piRNA_cluster_analysis

