#!/bin/bash

#toolDir='/home/yk42w/codes/ViralGenomeAssembly'
toolDir='/Users/yasinkaymaz/Dropbox/codes/ViralGenomeAssembly'

InputMultiFasta=$1
strand=$2
#How to run this: bash RC_multiFasta.sh inputfasta.fasta +


#<--->
perl $toolDir/bin/fasta_to_tab.pl "$InputMultiFasta" > "${InputMultiFasta%.fasta}".tab;

#is gene expressed from minus transcript? if yes, take reverse transcript of multi-fasta file.
if [ "$strand" = "-" ]
then
  #if gene is from reverse strand:
  while read ll;
  do
    sample_name=`echo $ll|cut -d " " -f1`;
    seq=`echo $ll|cut -d " " -f2`;
    rcseq=`echo $seq|rev | tr "ATGCNatgcn" "TACGNTACGN"`;
    echo -e "$sample_name\t$rcseq";
  done < "${InputMultiFasta%.fasta}".tab > "${InputMultiFasta%.fasta}".rc.tab;

  awk '{print ">"$1"\n"$2 }' "${InputMultiFasta%.fasta}".rc.tab > "${InputMultiFasta%.fasta}".rc.fasta

#if no, keep multi-fasta as it is. just convert to upper case
else

  while read ll;
  do
    sample_name=`echo $ll|cut -d " " -f1`;
    seq=`echo $ll|cut -d " " -f2`;
    rcseq=`echo $seq| tr "ATGCNatgcn" "ATGCNATGCN"`;
    echo -e "$sample_name\t$rcseq";
  done < "${InputMultiFasta%.fasta}".tab > "${InputMultiFasta%.fasta}".c.tab;

  awk '{print ">"$1"\n"$2 }' "${InputMultiFasta%.fasta}".c.tab > "${InputMultiFasta%.fasta}".c.fasta

fi




if [ "$splitTypes" = "1" ]
then

  grep -w Type1 $toolDir/workspace/data/pop_info.txt | awk 'NR==FNR{a[$1]=$0;next}{print a[$1]}' "$geneOfvariant".aln.c.tab - | awk '{print ">"$1"\n"$2 }' > "$geneOfvariant".aln.c.Type1.fasta;
  Type1multifasta="$geneOfvariant".aln.c.Type1.fasta;
  grep -w Type2 $toolDir/workspace/data/pop_info.txt | awk 'NR==FNR{a[$1]=$0;next}{print a[$1]}' "$geneOfvariant".aln.c.tab - | awk '{print ">"$1"\n"$2 }' > "$geneOfvariant".aln.c.Type2.fasta;
  Type2multifasta="$geneOfvariant".aln.c.Type2.fasta;

  genemultifasta="$geneOfvariant".aln.fasta
else
  echo "done"
fi


rm "${InputMultiFasta%.fasta}".*.tab
#<--->
