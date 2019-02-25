


# Instruction

Is's a series of small RNA-seq analysis in SP lab.
I will briefly describe the general flow and the function of scripts in below section.

I might update some scripts of other analysis in the future.


# Small RNA-seq

## General Flow

In general, the analysis is carried out in below order.
There are different scripts in those directories. Moreover, the detailed instruction are written in scripts respectively.

### 1. pipeline

In this step, you will map your sequencing data to genome.

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





### 3. piRNA_character_analysis

### 4. piRNA_cluster_analysis

