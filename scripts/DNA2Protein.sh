#!/bin/bash
dir=`pwd`

#toolDir='/home/yk42w/codes/ViralGenomeAssembly'
toolDir='/Users/yasinkaymaz/Dropbox/codes/ViralGenomeAssembly'

type=1
nt=2
VariantFilterer=1
#EBV genome type
if [ "$type" = "1" ]
then
echo "Genome is Type I"
Type='type1'
EBVgenome="NC_007605"
VfatDIR="$toolDir/resources/Vfat"
RefDIR="$toolDir/resources/Annotation/Type1/Vfat/EBV"
AnnotationDir="$toolDir/resources/Annotation/Type1/"
RefGenesDir="$toolDir/resources/Annotation/Type1/"
else
echo "Genome is Type II"
Type='type2'
EBVgenome="NC_009334"
VfatDIR="$toolDir/resources/Vfat/Type2"
RefDIR="$toolDir/resources/Annotation/Type2/Vfat/EBV"
AnnotationDir="$toolDir/resources/Annotation/Type2/"
fi

#Take input dna multi fasta file. this file has to end with .fasta
InputFasta=$1

#Convert input fasta to tab file
perl $toolDir/bin/fasta_to_tab.pl "$InputFasta" > "$InputFasta".aln.tab;

while read line;
do
  #Transform lowercase bases to UPPERCASEs only sequences
  sample_name=`echo $line|cut -d " " -f1`;
  seq=`echo $line|cut -d " " -f2`;
  newseq=`echo $seq | tr "atgcn" "ATGCN"`;
  #echo -e "$sample_name $newseq";
  #Revert tab to fasta for a given line.
  echo -e ">$sample_name\n$newseq" > dna.fa
  #translate DNA to protein.
  perl $toolDir/bin/Vfat/translateDna.pl dna.fa protein.fa
  #write the output
  cat protein.fa >> "${InputFasta%.fasta}".protein.fasta

done < "$InputFasta".aln.tab;

rm "$InputFasta".aln.tab protein.fa dna.fa;

#~/Documents/Tools/weblogo/seqlogo -k 0 -M  -c -n -Y  -F PDF -f epitopes.bed.aln.fasta > test.pdf
#weblogo -A protein -c hydrophobicity -l 10 -u 20 -D fasta -F pdf -f my_ICed_EBNA3C_protein.fasta -o test2.pdf
