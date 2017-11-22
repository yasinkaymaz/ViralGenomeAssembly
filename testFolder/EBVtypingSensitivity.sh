#!/bin/bash

toolDir='/home/yk42w/codes/ViralGenomeAssembly'
module unload python/2.7.5
module load python/2.7.9
module load python/2.7.9_packages/pandas/0.17.1
for genome in `grep ">" sequences.aln.fasta|grep -v NC_|sed 's/>//g'`;
do
  echo $genome;
  python $toolDir/bin/MSA_parser_cleaner.py $genome NC_007605 sequences.aln.fasta 1;
  python $toolDir/bin/MSA_parser_cleaner.py $genome NC_009334 sequences.aln.fasta 2;
done
