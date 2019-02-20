
#!/bin/bash

#Program
#	Apply piRNA eland3 file to proTRAC to predict piRNA clusters.
#Usage
#	Copy this script, link genome fasta to where your piRNA eland3 file locating.
#	Modify "inputDir" to where your piRNA eland3 file locating.
#	Execute.
#History
#2017/02/22, Hsien-Hen, Lin
#	Written.
#2019/01/28, HHL
#	Modified for human datasets.

inputDir="$HOME/human_gonad_piRNA/data/20190128/20130502/eland3_file/"
programDir=

for ELAND3 in $(ls $inputDir|grep '.eland3')

do
	echo "start proTRAC of $ELAND3";
	perl $programDir/proTRAC_2.4.2.kw_mod.pl -genome $inputDir/hg38.fa -map $inputDir/$ELAND3 -pdens 0.05 -1Tor10A 0.01 -1Tand10A 0.01 -distr 1-90 -clsize 1000 -clsplit 0.2 -pimin 24 -pimax 34 -clhitsn 100 -nofaspi;
	wait;
	echo "finish proTRAC of $ELAND3"
done



