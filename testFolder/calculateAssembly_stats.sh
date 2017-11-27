#!/bin/bash

toolDir='/home/yk42w/codes/ViralGenomeAssembly'

for i in `ls -1 *_genome_fixed_assembly_genome.fa`;
do
  echo -e "${i%_genome_fixed_assembly_genome.fa}" >> n50stats.txt
  python $toolDir/bin/Assembly_Splitter_from_Ns.py $i 1 test.fa;
  python $toolDir/bin/Fasta2N50.py test.fa >> n50stats.txt;
done
