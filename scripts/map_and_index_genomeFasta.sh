
type=$1
InputGenomeFasta=$2
SAMPLE_NAME=$3

toolDir='/home/yk42w/codes/ViralGenomeAssembly'
PICARDPATH="$toolDir/bin/mac/picard-2.14.1/build/libs"
Bowtie2PATH="$toolDir/bin/mac/bowtie2-2.3.3.1-macos-x86_64"
nt=2

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

python $toolDir/bin/Assembly_Splitter_from_Ns.py $InputGenomeFasta 1 "${InputGenomeFasta%.fa}"_contigs.fa.tmp

seqnum=`grep ">" "${InputGenomeFasta%.fa}"_contigs.fa.tmp |wc -l`;

for ((i=1;i<=$seqnum;i++));
  do
    let x=$i*2;
    head -$x "${InputGenomeFasta%.fa}"_contigs.fa.tmp |\
    tail -2| $toolDir/bin/fasta_chop.pl;
  done > "${InputGenomeFasta%.fa}"_contigs_chopped.fa

$Bowtie2PATH/bowtie2 -p $nt -f -x $toolDir/resources/Bowtie2Index/$EBVgenome \
-U "${InputGenomeFasta%.fa}"_contigs_chopped.fa \
-S "${InputGenomeFasta%.fa}"_mapped2_"$Type".sam

#awk '($5 != 0)' "${InputFasta%.fa}"_mapped2_"$Type".sam > "${InputFasta%.fa}"_mapped2_"$Type"_noMQ0.sam

#Sort sam file and output as bam
java -Xmx10g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/picard.jar \
SortSam \
INPUT="${InputGenomeFasta%.fa}"_mapped2_"$Type".sam \
OUTPUT="${InputGenomeFasta%.fa}"_mapped2_"$Type"_sorted.bam \
SORT_ORDER=coordinate
#Index bam file
java -Xmx10g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/picard.jar \
BuildBamIndex \
INPUT="${InputGenomeFasta%.fa}"_mapped2_"$Type"_sorted.bam


rm *tmp
