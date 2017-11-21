#!/bin/bash

toolDir='/home/yk42w/codes/RNAseq_tools'
export PATH=$toolDir/bin/linux/kallisto_linux-v0.43.1:$PATH

SAMPLE_NAME=$1
#stranded or nonstranded
LibraryType=$2
#number of threads to use.
nt=$3
Fastq1=$4
Fastq2=$5
#organism is ether "human" or "EBV"
Organism=$6
#optional EBV type. Required only when EBV is used. Either 1 or 2.
type=$7

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


if [ -z "$Fastq2" ] ; then
  Pairness='--single -l 200 -s 20'
else
  Pairness=''
fi

if [ Organism == "human" ];
then
  index_file=''
else
  index_file=$EBVINDEX
fi

kallisto quant -i $index_file \
-o kallisto."$SAMPLE_NAME" \
--rf-stranded \
$Pairness \
-t $nt \
--bias \
$Fastq1 $Fastq2
