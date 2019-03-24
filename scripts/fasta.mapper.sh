#!/bin/bash
#module load bowtie2/2.3.2-fasrc02

type=$1
InputSam=$2
SAMPLE_NAME=$3

#cat BL643_eBLFNAtumorDNA_wga_velvet_all_type2.sam|grep GGGGACTTTATGTGACCCTTGGGCCT |awk '{if(  ($3 == "*") ) print ">"$1"\n"$10 }' |less

awk '{if(($3 == "*")) print ">"$1"\n"$10 }' $InputSam > "${InputSam%.sam}".fa

InputGenomeFasta="${InputSam%.sam}".fa

toolDir='~/codes/ViralGenomeAssembly'
PICARDPATH="$toolDir/bin/linux/picard-2.14.1/build/libs"
Bowtie2PATH="$toolDir/bin/linux/bowtie2-2.3.3.1-linux-x86_64"
nt=2

export PATH="${Bowtie2PATH}":"${toolDir}":"${toolDir}/bin/":"${PICARDPATH}":${PATH}

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

#python $toolDir/bin/Assembly_Splitter_from_Ns.py $InputGenomeFasta 1 "${InputGenomeFasta%.fa}"_contigs.fa.tmp

#seqnum=`grep ">" "${InputGenomeFasta%.fa}"_contigs.fa.tmp |wc -l`;
seqnum=`grep ">" $InputGenomeFasta |wc -l`;

for ((i=1;i<=$seqnum;i++));
  do
    let x=$i*2;
    head -$x $InputGenomeFasta | tail -2 | perl fasta_chop.pl;
done > "${InputGenomeFasta%.fa}"_contigs_chopped.fa

bowtie2 -p $nt -f -x $toolDir/resources/Bowtie2Index/$EBVgenome \
-U "${InputGenomeFasta%.fa}"_contigs_chopped.fa \
-S "${InputGenomeFasta%.fa}"_mapped2_"$Type".sam

#awk '($5 != 0)' "${InputFasta%.fa}"_mapped2_"$Type".sam > "${InputFasta%.fa}"_mapped2_"$Type"_noMQ0.sam

#Sort sam file and output as bam
java -Xmx10g -XX:ParallelGCThreads=$nt -jar \
picard.jar \
SortSam \
INPUT="${InputGenomeFasta%.fa}"_mapped2_"$Type".sam \
OUTPUT="${InputGenomeFasta%.fa}"_mapped2_"$Type"_sorted.bam \
SORT_ORDER=coordinate
#Index bam file
java -Xmx10g -XX:ParallelGCThreads=$nt -jar \
picard.jar \
BuildBamIndex \
INPUT="${InputGenomeFasta%.fa}"_mapped2_"$Type"_sorted.bam


rm *tmp
