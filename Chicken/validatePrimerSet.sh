#!/bin/bash

input=$1;
refBed=$2;

if [[ $# -ne 2 ]]; then
  echo '    Usage: validatePrimerSet.sh <input> <reference.bed>';
  echo '    input = file contain result text copied from UCSC in silico PCR';
  echo '    reference.bed = TE reference in bed format';
  exit;
fi

grep "^>" $input |sed -e 's/>//g' -e 's/:/ /g' -e 's/\(+\|-\)/ \1 /g'|awk 'BEGIN{OFS="\t"}{print $1,$2,$4,"Read_"NR,0,$3}'|intersectBed -s -a $refBed -b -|awk '{print $4}'|sort|uniq

