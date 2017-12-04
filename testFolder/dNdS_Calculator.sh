#!/bin/bash

toolDir='/home/yk42w/codes/ViralGenomeAssembly'
module unload python/2.7.5
module load python/2.7.9
module load python/2.7.9_packages/pandas/0.17.1

#Alignment file has to end with ".aln.fasta"
AlignmentFile=$1
RefGenesDir=$toolDir/resources/Annotation/Type1/

#for gene in `ls -1 $RefGenesDir|sed 's/.bed//g'`;
for gene in `cut -f4 $RefGenesDir/EBV_Reference_genelist_Genenames_stranded.txt|uniq`;
do
  echo $gene;
  grep -w $gene $RefGenesDir/EBV_Reference_genelist_Genenames_stranded.txt > "$gene".bed;
  strand=`cut -f5 "$gene".bed|uniq`;
  python $toolDir/bin/MSA_gene_extractor.py $AlignmentFile 1 "$gene".bed;
  perl $toolDir/bin/fasta_to_tab.pl "$gene".bed.aln.fasta > "$gene".aln.tab;
  RefORFlen=`grep NC_007605 "$gene".aln.tab |sed 's/[n|-]//g'|awk '{print length($2)}'`;
  #Exclude genome samples if the given gene is not covered more than %50:
  awk -v RefORFlen="$RefORFlen" '{seq=$2;gsub(/n|-/,"",seq); if ( 100*length(seq)/RefORFlen > 50) print}' "$gene".aln.tab > "$gene".aln.filtered.tab
  if $strand == "-";
  then
    #if gene is from reverse strand:
    while read line;
    do
      sample_name=`echo $line|cut -d " " -f1`;
      seq=`echo $line|cut -d " " -f2`;
      rcseq=`echo $seq|rev | tr "ATGCNatgcn" "TACGNtacgn"`;
      echo -e "$sample_name $rcseq";
    done < "$gene".aln.filtered.tab > "$gene".aln.filtered.rc.tab
    testalignmentFile="$gene".aln.filtered.rc.tab
  else
    #else
    testalignmentFile="$gene".aln.filtered.tab
  fi

  perl $toolDir/bin/SNAP.pl "$testalignmentFile";
  grep all summary.* |awk -v gene="$gene" '{print gene"\t"$0 }' >> EBV_Genes_dNdS_summary_"${AlignmentFile%.aln.fasta}".txt;
  rm summary.*;
done
