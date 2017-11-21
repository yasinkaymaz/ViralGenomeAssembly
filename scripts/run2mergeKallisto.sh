#!/bin/bash

#form=isoforms
form=genes

toolDir='/home/yk42w/codes/ViralGenomeAssembly'

#Merge RSEM output files:
utils=/n/home13/yasinkaymaz/codes/utils
for type in 1 2;
do

  #EBV genome type
  if [ "$type" = "1" ]
  then
  echo "Genome is Type I"
  Type='type1'
  EBVgenome="NC_007605"
  VfatDIR="$toolDir/resources/Vfat"
  RefDIR="$toolDir/resources/Annotation/Type1/Vfat/EBV"
  AnnotationDir="$toolDir/resources/Annotation/Type1/"
  EBVINDEX="$toolDir/resources/RNAseq/kallisto/Type1/EBV_transcripts"
  else
  echo "Genome is Type II"
  Type='type2'
  EBVgenome="NC_009334"
  VfatDIR="$toolDir/resources/Vfat/Type2"
  RefDIR="$toolDir/resources/Annotation/Type2/Vfat/EBV"
  AnnotationDir="$toolDir/resources/Annotation/Type2/"
  EBVINDEX="$toolDir/resources/RNAseq/kallisto/Type2/EBV_transcripts"
  fi

  for value in TPM;
  do
$toolDir/utils/merge_Kallisto-outputs_single_table_${value}.pl \
kallisto.eBL01_"$type"_PolyA-RNA_EBV/eBL01_"$type"_PolyA-RNA \
kallisto.eBL01_"$type"_TotalRNA_EBV/eBL01_"$type"_TotalRNA \
kallisto.eBL01P_"$type"_TotalRNA_EBV/eBL01P_"$type"_TotalRNA \
kallisto.eBL02_"$type"_PolyA-RNA_EBV/eBL02_"$type"_PolyA-RNA \
kallisto.eBL02_"$type"_TotalRNA_EBV/eBL02_"$type"_TotalRNA \
kallisto.eBL02P_"$type"_TotalRNA_EBV/eBL02P_"$type"_TotalRNA \
kallisto.eBL03_"$type"_PolyA-RNA_EBV/eBL03_"$type"_PolyA-RNA \
kallisto.eBL03_"$type"_TotalRNA_EBV/eBL03_"$type"_TotalRNA \
kallisto.eBL03P_"$type"_TotalRNA_EBV/eBL03P_"$type"_TotalRNA \
kallisto.eBL04_"$type"_PolyA-RNA_EBV/eBL04_"$type"_PolyA-RNA \
kallisto.eBL04_"$type"_TotalRNA_EBV/eBL04_"$type"_TotalRNA \
kallisto.eBL04P_"$type"_TotalRNA_EBV/eBL04P_"$type"_TotalRNA \
kallisto.eBL05_"$type"_PolyA-RNA_EBV/eBL05_"$type"_PolyA-RNA \
kallisto.eBL05_"$type"_TotalRNA_EBV/eBL05_"$type"_TotalRNA \
kallisto.eBL06_"$type"_PolyA-RNA_EBV/eBL06_"$type"_PolyA-RNA \
kallisto.eBL06_"$type"_TotalRNA_EBV/eBL06_"$type"_TotalRNA \
kallisto.eBL07_"$type"_PolyA-RNA_EBV/eBL07_"$type"_PolyA-RNA \
kallisto.eBL07_"$type"_TotalRNA_EBV/eBL07_"$type"_TotalRNA \
kallisto.eBL08_"$type"_PolyA-RNA_EBV/eBL08_"$type"_PolyA-RNA \
kallisto.eBL08_"$type"_TotalRNA_EBV/eBL08_"$type"_TotalRNA \
kallisto.eBL09_"$type"_PolyA-RNA_EBV/eBL09_"$type"_PolyA-RNA \
kallisto.eBL09_"$type"_TotalRNA_EBV/eBL09_"$type"_TotalRNA \
kallisto.eBL10_"$type"_PolyA-RNA_EBV/eBL10_"$type"_PolyA-RNA \
kallisto.eBL10_"$type"_TotalRNA_EBV/eBL10_"$type"_TotalRNA \
kallisto.eBL11_"$type"_PolyA-RNA_EBV/eBL11_"$type"_PolyA-RNA \
kallisto.eBL11_"$type"_TotalRNA_EBV/eBL11_"$type"_TotalRNA \
kallisto.eBL12_"$type"_PolyA-RNA_EBV/eBL12_"$type"_PolyA-RNA \
kallisto.eBL12_"$type"_TotalRNA_EBV/eBL12_"$type"_TotalRNA \
kallisto.eBL13_"$type"_PolyA-RNA_EBV/eBL13_"$type"_PolyA-RNA \
kallisto.eBL13_"$type"_TotalRNA_EBV/eBL13_"$type"_TotalRNA \
kallisto.eBL14_"$type"_PolyA-RNA_EBV/eBL14_"$type"_PolyA-RNA \
kallisto.eBL14_"$type"_TotalRNA_EBV/eBL14_"$type"_TotalRNA \
kallisto.eBL17_"$type"_TotalRNA_EBV/eBL17_"$type"_TotalRNA \
kallisto.eBL19_"$type"_PolyA-RNA_EBV/eBL19_"$type"_PolyA-RNA \
kallisto.eBL19_"$type"_TotalRNA_EBV/eBL19_"$type"_TotalRNA \
kallisto.eBL20_"$type"_PolyA-RNA_EBV/eBL20_"$type"_PolyA-RNA \
kallisto.eBL20_"$type"_TotalRNA_EBV/eBL20_"$type"_TotalRNA \
kallisto.eBL21_"$type"_PolyA-RNA_EBV/eBL21_"$type"_PolyA-RNA \
kallisto.eBL21_"$type"_TotalRNA_EBV/eBL21_"$type"_TotalRNA \
kallisto.eBL22_"$type"_PolyA-RNA_EBV/eBL22_"$type"_PolyA-RNA \
kallisto.eBL22_"$type"_TotalRNA_EBV/eBL22_"$type"_TotalRNA \
kallisto.eBL23_"$type"_PolyA-RNA_EBV/eBL23_"$type"_PolyA-RNA \
kallisto.eBL23_"$type"_TotalRNA_EBV/eBL23_"$type"_TotalRNA \
kallisto.eBL24_"$type"_PolyA-RNA_EBV/eBL24_"$type"_PolyA-RNA \
kallisto.eBL24_"$type"_TotalRNA_EBV/eBL24_"$type"_TotalRNA \
kallisto.eBL25_"$type"_PolyA-RNA_EBV/eBL25_"$type"_PolyA-RNA \
kallisto.eBL25_"$type"_TotalRNA_EBV/eBL25_"$type"_TotalRNA \
kallisto.eBL26_"$type"_PolyA-RNA_EBV/eBL26_"$type"_PolyA-RNA \
kallisto.eBL26_"$type"_TotalRNA_EBV/eBL26_"$type"_TotalRNA \
kallisto.eBL27_"$type"_PolyA-RNA_EBV/eBL27_"$type"_PolyA-RNA \
kallisto.eBL27_"$type"_TotalRNA_EBV/eBL27_"$type"_TotalRNA \
kallisto.eBL28_"$type"_PolyA-RNA_EBV/eBL28_"$type"_PolyA-RNA \
kallisto.eBL28_"$type"_TotalRNA_EBV/eBL28_"$type"_TotalRNA \
kallisto.eBL29_"$type"_PolyA-RNA_EBV/eBL29_"$type"_PolyA-RNA \
kallisto.eBL29_"$type"_TotalRNA_EBV/eBL29_"$type"_TotalRNA \
kallisto.eBL30_"$type"_PolyA-RNA_EBV/eBL30_"$type"_PolyA-RNA \
kallisto.eBL30_"$type"_TotalRNA_EBV/eBL30_"$type"_TotalRNA \
kallisto.eBL31_"$type"_PolyA-RNA_EBV/eBL31_"$type"_PolyA-RNA \
kallisto.eBL31_"$type"_TotalRNA_EBV/eBL31_"$type"_TotalRNA \
kallisto.eBL32_"$type"_PolyA-RNA_EBV/eBL32_"$type"_PolyA-RNA \
kallisto.eBL32_"$type"_TotalRNA_EBV/eBL32_"$type"_TotalRNA \
kallisto.eBL33_"$type"_PolyA-RNA_EBV/eBL33_"$type"_PolyA-RNA \
kallisto.eBL33_"$type"_TotalRNA_EBV/eBL33_"$type"_TotalRNA \
kallisto.eBL34_"$type"_PolyA-RNA_EBV/eBL34_"$type"_PolyA-RNA \
kallisto.eBL34_"$type"_TotalRNA_EBV/eBL34_"$type"_TotalRNA \
kallisto.eBL35_"$type"_PolyA-RNA_EBV/eBL35_"$type"_PolyA-RNA \
kallisto.eBL35_"$type"_TotalRNA_EBV/eBL35_"$type"_TotalRNA \
kallisto.eBL36_"$type"_PolyA-RNA_EBV/eBL36_"$type"_PolyA-RNA \
kallisto.eBL36_"$type"_TotalRNA_EBV/eBL36_"$type"_TotalRNA \
kallisto.eBL37_"$type"_PolyA-RNA_EBV/eBL37_"$type"_PolyA-RNA \
kallisto.eBL37_"$type"_TotalRNA_EBV/eBL37_"$type"_TotalRNA \
kallisto.eBL38_"$type"_PolyA-RNA_EBV/eBL38_"$type"_PolyA-RNA \
kallisto.eBL38_"$type"_TotalRNA_EBV/eBL38_"$type"_TotalRNA \
kallisto.eBL39_"$type"_PolyA-RNA_EBV/eBL39_"$type"_PolyA-RNA \
kallisto.eBL39_"$type"_TotalRNA_EBV/eBL39_"$type"_TotalRNA \
kallisto.eBL42_"$type"_PolyA-RNA_EBV/eBL42_"$type"_PolyA-RNA \
kallisto.eBL42_"$type"_TotalRNA_EBV/eBL42_"$type"_TotalRNA \
kallisto.eBL44_"$type"_PolyA-RNA_EBV/eBL44_"$type"_PolyA-RNA \
kallisto.eBL45_"$type"_TotalRNA_EBV/eBL45_"$type"_TotalRNA \
kallisto.eBL46_"$type"_TotalRNA_EBV/eBL46_"$type"_TotalRNA \
kallisto.eBL47_"$type"_PolyA-RNA_EBV/eBL47_"$type"_PolyA-RNA \
kallisto.eBL47_"$type"_TotalRNA_EBV/eBL47_"$type"_TotalRNA \
kallisto.eBL48_"$type"_PolyA-RNA_EBV/eBL48_"$type"_PolyA-RNA \
kallisto.eBL48_"$type"_TotalRNA_EBV/eBL48_"$type"_TotalRNA \
> merged.kallisto.expression_"$Type"_${value}.txt

awk 'FNR==NR{a[substr($1,1,15)]=$1;next}{print $0"\t", a[$1]}' $toolDir/resources/Annotation/$Type/isoforms2genes.tab merged.kallisto.expression_"$Type"_"${value}".txt > merged.kallisto.expression_"$Type"_"${value}"_geneIds.txt
Rscript $toolDir/utils/MergeTranscriptsToGenes.R merged.kallisto.expression_"$Type"_"${value}"_geneIds.txt
#cat /n/home13/yasinkaymaz/LabSpace/results/test/merged.kallisto.expression_${value}.txt|sed 's/scRNAseq_paired_//g'|sed 's/.genes.results//g'|cut -d "_" -f2- > /n/home13/yasinkaymaz/LabSpace/results/test/Merged.kallisto.Expression_${value}.txt
echo "...";
done
done
