#!/bin/bash

dir=`pwd`
InputAssembly=$1
InputBAM=$2
SAMPLE_NAME=$3
type=$4
toolDir='~/codes/ViralGenomeAssembly'
nt=10

Reads2Assembly=1
VarCall=1

. /n/home13/yasinkaymaz/miniconda3/etc/profile.d/conda.sh

conda activate py27
source new-modules.sh
module load R/3.4.2-fasrc01
module load gcc/7.1.0-fasrc01
module purge
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


if [ "$Reads2Assembly" = "1" ]
then

#sort bam file based on read name
picard \
FixMateInformation \
INPUT="$InputBAM" \
OUTPUT="${InputBAM%.bam}"_sorted.bam \
SORT_ORDER=queryname \
VALIDATION_STRINGENCY=LENIENT

#Dump paired fastq files from bam file
bamToFastq -i "${InputBAM%.bam}"_sorted.bam \
-fq "${InputBAM%.bam}"_1.fastq \
-fq2 "${InputBAM%.bam}"_2.fastq

#create bowtie2 index
bowtie2-build $InputAssembly ${InputAssembly%.fa}_bowtieIndex
samtools faidx $InputAssembly

#Map reads back to assembly
bowtie2 -p $nt \
-x ${InputAssembly%.fa}_bowtieIndex \
-1 "${InputBAM%.bam}"_1.fastq \
-2 "${InputBAM%.bam}"_2.fastq \
-S "${InputBAM%.bam}"_readsBack2assembly.sam

#Filter MQ<20 reads
samtools view -h -b -q 30 \
"${InputBAM%.bam}"_readsBack2assembly.sam \
-o "${InputBAM%.bam}"_readsBack2assembly_MQ30.bam

#Sort sam file and output as bam
picard \
SortSam \
INPUT="${InputBAM%.bam}"_readsBack2assembly_MQ30.bam \
OUTPUT="${InputBAM%.bam}"_readsBack2assembly_MQ30_sorted.bam \
SORT_ORDER=coordinate

picard \
BuildBamIndex \
INPUT="${InputBAM%.bam}"_readsBack2assembly_MQ30_sorted.bam

else
echo "no pileup"
fi

if [ "$VarCall" = "1" ]
then
  samtools mpileup -uf $InputAssembly "${InputBAM%.bam}"_readsBack2assembly_MQ30_sorted.bam |\
  bcftools call -c \
  --output-type v \
  --variants-only \
  --output "${InputBAM%.bam}".Assembly.variants.vcf
  #
  cat "${InputBAM%.bam}".Assembly.variants.vcf| \
  grep -v "#"|awk '($6 > 100)'|\
  awk 'NF { info=$8; gsub(/.*;DP4=|;MQ=.*/, "", info); split(info, a, /,/); print $0 "\t" a[1]"\t"a[2]"\t"a[3]"\t"a[4]}' |\
  awk '{if($NF+$(NF-1) != 0) print $1"_"$2"_"$4"_"$5 "\t" ($(NF)+$(NF-1)+$(NF-2)+$(NF-3)) "\t" ($NF+$(NF-1)) / ($(NF)+$(NF-1)+$(NF-2)+$(NF-3))}'|\
  awk '($2 > 99)'|\
  sed 's/,/_/g' > "${InputBAM%.bam}".AssemblyMinorfreq.txt
  #
  Rscript $toolDir/bin/MAF_plot.R $SAMPLE_NAME "${InputBAM%.bam}".AssemblyMinorfreq.txt

else
echo "no pileup"
fi

if [ -f "${InputBAM%.bam}".AssemblyMinorfreq.txt ]
then

rm "${InputBAM%.bam}"_readsBack2assembly_MQ30.bam \
"${InputBAM%.bam}"_readsBack2assembly.sam \
${InputAssembly%.fa}_bowtieIndex* \
"${InputBAM%.bam}"_*.fastq \
"${InputBAM%.bam}"_sorted.bam

else
echo "no pileup"
fi
