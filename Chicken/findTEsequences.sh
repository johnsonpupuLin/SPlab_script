#!/bin/bash

TE=$1;
refBED=$2;
refGenome=$3;

if [[ $# -ne 3 ]]; then
  echo '    Usage: findTEsequences.sh <TE name> <reference.bed> <genome.fa>';
  echo '    TE name = name of TE you wish to find';
  echo '    reference.bed = TE reference in bed format';
  echo '    genome.fa = Reference Genome in fasta format';
  exit;
fi

grep -e $TE"\|*" -e $TE -w $refBED | fastaFromBed -name -split -s -fi $refGenome -bed - -fo $TE.fasta;

if [[ -s $TE.fasta ]]; then
  echo '    Sequences of '$TE' are succesfully retrived in '$TE.fasta;  
else
  echo '    No identifiable '$TE' in the reference,';
  echo '    or Invalid inputs';
  rm $TE.fasta;
fi

