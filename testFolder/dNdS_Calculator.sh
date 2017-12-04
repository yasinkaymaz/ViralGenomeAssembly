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
  RefORFlen=`grep NC_007605 "$gene".aln.tab |sed 's/[n|-]//g'|awk '{print length($2)}'`;
  #Exclude genome samples if the given gene is not covered more than %50:
  awk -v RefORFlen="$RefORFlen" '{seq=$2;gsub(/n|-/,"",seq); if ( 100*length(seq)/RefORFlen > 50) print}' "$gene".aln.tab > "$gene".aln.filtered.tab
  perl $toolDir/bin/SNAP.pl "$gene".aln.filtered.tab;
  grep all summary.* |awk -v gene="$gene" '{print gene"\t"$0 }' >> EBV_Genes_dNdS_summary_"${AlignmentFile%.aln.fasta}".txt;
  rm summary.*;
done
