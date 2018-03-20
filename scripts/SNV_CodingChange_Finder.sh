#!/bin/bash
dir=`pwd`

toolDir='/home/yk42w/codes/ViralGenomeAssembly'
type=1
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

GATKdir="/project/umw_jeffrey_bailey/share/bin_sync/"

module load java/1.7.0_25
module load IGVTools/2.3.31

InputVCF=$1

if [ "$CreateAlternateFasta" = "1" ]
then

igvtools index $InputVCF

java -Xmx16g -XX:ParallelGCThreads=$nt -jar \
$GATKdir/GenomeAnalysisTK.jar \
-R $toolDir/resources/Bowtie2Index/$EBVgenome.fa \
-T FastaAlternateReferenceMaker \
--variant $InputVCF \
-o AlternativeGenome_with_SNVs.fasta

else
echo "Skipping Create Alternate Fasta"
fi
