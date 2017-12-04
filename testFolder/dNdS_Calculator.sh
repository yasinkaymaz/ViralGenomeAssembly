#!/bin/bash

toolDir='/home/yk42w/codes/ViralGenomeAssembly'
module unload python/2.7.5
module load python/2.7.9
module load python/2.7.9_packages/pandas/0.17.1

#Alignment file has to end with ".aln.fasta"
AlignmentFile=$1
RefGenesDir=$toolDir/resources/Annotation/Type1/ORFs

for gene in `ls -1 $RefGenesDir|sed 's/.bed//g'`;
do
  echo $gene;
  python $toolDir/bin/MSA_gene_extractor.py $AlignmentFile 1 $RefGenesDir/"$gene".bed;
  perl $toolDir/bin/fasta_to_tab.pl "$gene".bed.aln.fasta > "$gene".aln.tab;
  perl $toolDir/bin/SNAP.pl "$gene".aln.tab;
  grep all summary.* |awk -v gene="$gene" '{print gene"\t"$0 }' >> EBV_Genes_dNdS_summary_"${AlignmentFile%.aln.fasta}".txt;
  rm summary.*;
done
