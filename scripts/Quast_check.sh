#!/bin/bash

InputDir=$1
#type is either 1 or 2
type=$2
nt=2

QuastDir='/project/umw_jeffrey_bailey/share/bin_sync/quast-4.4/'
toolDir='/home/yk42w/codes/ViralGenomeAssembly'
INDEXGENOMESDIR="$toolDir/resources/Bowtie2Index"

#EBV genome type
if [ "$type" = "1" ]
then
echo "Genome is Type I"
Type='type1'
EBVgenome="NC_007605"
VfatDIR="$toolDir/resources/Vfat"
RefDIR="$toolDir/resources/Annotation/Type1/Vfat/EBV"
AnnotationDir="$toolDir/resources/Annotation/Type1/"
else
echo "Genome is Type II"
Type='type2'
EBVgenome="NC_009334"
VfatDIR="$toolDir/resources/Vfat/Type2"
RefDIR="$toolDir/resources/Annotation/Type2/Vfat/EBV"
AnnotationDir="$toolDir/resources/Annotation/Type2/"
fi

module load python/2.7.5
module load gcc/4.8.1
module load jdk/1.7.0_25
module load boost/1.57.0
module load MUMmer/3.23

$QuastDir/quast.py $InputDir/*.fa \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
-G $AnnotationDir/$EBVgenome.gff \
-m 150 \
--threads $nt \
--scaffolds \
-o $InputDir/Quast_Optimization_$Type

#--scaffolds



#bsub -q long -n 4 -R rusage[mem=12000] -W 04:00 ~/codes/EBVseq/Quast_check.sh ./ 1
