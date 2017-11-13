#!/bin/bash

# Set the number of threads for tools that can be run in parallel.
#This is usually one digit below the number of cores this script is reserved to use (-n value above).


#Set the working directory and Sample name
dir=`pwd`

Batch=$1

if [ "$Batch" = "1" ]
then
SAMPLE_NAME=$2
Kmer_start=$3
Kmer_end=$4
nt=$5
option=$6
type=$7
else
#submit this way: bsub -q long -n 6 -R rusage[mem=20000] -R "span[hosts=1]" -R "select[tmp>1000]" -W 96:00 -e err.%J.txt -o out.%J.txt ~/codes/EBVseq/EBVGenomeTyper.sh 0
SAMPLE_NAME=`basename $dir`
Kmer_start=21
Kmer_end=151
nt=6
option='n50'

fi


echo "job_id: $LSB_JOBID"
echo "index_id: $LSB_JOBINDEX"
echo "job_name: $LSB_JOBNAME"

################################################
############# Main Control Panel: ##############
################################################

#Setting the directories for required files
WORKINGDIR=$dir

#ls -1 *R1_00*|sort|xargs -i zcat {} >> "$SAMPLE_NAME"_R1.fastq
#ls -1 *R2_00*|sort|xargs -i zcat {} >> "$SAMPLE_NAME"_R2.fastq

#gzip $WORKINGDIR/"$SAMPLE_NAME"_R*.fastq

fastq_1=$WORKINGDIR/"$SAMPLE_NAME"_R1.fastq.gz
fastq_2=$WORKINGDIR/"$SAMPLE_NAME"_R2.fastq.gz

##########################################
######### NO NEED TO MODIFY ! ############
############# Sub-routines ###############
##########################################


#Load required modules
module unload openssl/1.0.1g

module load ncbi_cxx/12_0_0
module load cutadapt/1.7.1
module load samtools/0.0.19
module load bowtie2/2-2.1.0
module load fastx_toolkit/0.0.14
module load java/1.7.0_25
module load fastqc/0.10.1
module load prinseq/0.20.4
module load gcc/4.8.1
module load bzip2/1.0.6
module load tabix/0.2.6
module load IGVTools/2.3.31
module load R/3.0.2
module load bedtools/2.17.0
module load vicuna/1.3
module load fastqc/0.10.1
module load perl/5.18.1

module load python/2.7.5

#Set the working environment
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/pkg/ncbi_cxx/12_0_0/lib
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/tagcleaner-standalone-0.16/:$PATH
export PERL5LIB=/project/umw_jeffrey_bailey/share/bin_sync/perl5lib/
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/Vfat:$PATH
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/picard-tools-1.105/:$PATH
export PERL5LIB=/project/umw_jeffrey_bailey/share/bin_sync/perl5lib/lib/perl5/
#export PATH=/project/umw_jeffrey_bailey/share/bin_sync/bamtools-master/bin/:$PATH
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/velvet_1.2.10:$PATH
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/VelvetOptimiser-2.2.5:$PATH
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/velvet_1.2.10/contrib/shuffleSequences_fasta:$PATH
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/MetaVelvet-1.2.02:$PATH
#Variables for required directories
toolDir='/project/umw_jeffrey_bailey/OTHERS/dnaSeq_Tools'

BCFtoolsDIR="/project/umw_jeffrey_bailey/share/bin_sync/samtools-0.1.19/bcftools"
PICARDPATH="/project/umw_jeffrey_bailey/share/bin_sync/picard-tools-1.105"
GATKdir="/project/umw_jeffrey_bailey/share/bin_sync/"
ShuffleDir='/project/umw_jeffrey_bailey/share/bin_sync/velvet_1.2.10/contrib/shuffleSequences_fasta'
INDEXGENOMESDIR="$toolDir/resources/Bowtie2Index"

ReadPreProc=0
ReadMapping2Ref=0
SimpleStats=1


if [ "$ReadPreProc" = "1" ]
then

#Do Fastqc check
#fastqc $fastq_1
#fastqc $fastq_2
################################################## TEMPORARY !!!
gunzip $WORKINGDIR/"$SAMPLE_NAME"_R*.fastq.gz

python /home/yk42w/codes/RandomSampler.py "$SAMPLE_NAME"_R1.fastq "$SAMPLE_NAME"_R2.fastq $WORKINGDIR/"$SAMPLE_NAME"_sub_1.fastq $WORKINGDIR/"$SAMPLE_NAME"_sub_2.fastq 10000000
fastq_1=$WORKINGDIR/"$SAMPLE_NAME"_sub_1.fastq
fastq_2=$WORKINGDIR/"$SAMPLE_NAME"_sub_2.fastq
gzip $WORKINGDIR/"$SAMPLE_NAME"_R*.fastq
################################################## #############
#Remove adaptors
cutadapt -a ATCTCGTATGCCGTCTTCTGCTTG \
-a GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCCGATGTTT \
-e 0.1 -O 5 \
-o "$SAMPLE_NAME"_AdpTrimmed_1.fastq \
-p "$SAMPLE_NAME"_AdpTrimmed_2.fastq \
$fastq_1 \
$fastq_2

fastq_1="$SAMPLE_NAME"_AdpTrimmed_1.fastq
fastq_2="$SAMPLE_NAME"_AdpTrimmed_2.fastq

#Do Fastqc check
#fastqc $fastq_1
#fastqc $fastq_2

#Trim low qual sequence off of reads
perl /share/pkg/prinseq/0.20.4/prinseq-lite.pl \
-fastq $fastq_1 -fastq2 $fastq_2 \
-out_good "$SAMPLE_NAME"_goodQual \
-out_bad "$SAMPLE_NAME"_lowQual \
-log \
-trim_qual_right 20

fastq_1="$SAMPLE_NAME"_goodQual_1.fastq
fastq_2="$SAMPLE_NAME"_goodQual_2.fastq

rm "$SAMPLE_NAME"_lowQual_*.fastq

cd $WORKINGDIR/
else
echo "Skipping Read Preprocessing step"
fi



typelist=(1 2)
for type in in "${typelist[@]}";
do
echo $type

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



#*********************************************************************************************************#
if [ "$ReadMapping2Ref" = "1" ]
then

fastq_1=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_1.fastq
fastq_2=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_2.fastq
# RUN ALIGNMENT TO REF GENOME WT
mkdir $WORKINGDIR/ReadsMap2Ref/

bowtie2 -p $nt -x $INDEXGENOMESDIR/$EBVgenome \
-q \
-1 $fastq_1 \
-2 $fastq_2 \
--rg-id 1000 --rg LB:HYBRIDCAPTURE --rg PL:ILLUMINA --rg SM:$SAMPLE_NAME \
-S $WORKINGDIR/ReadsMap2Ref/"$SAMPLE_NAME"_alignment_$EBVgenome.sam

#Sort sam file and output as bam
java -Xmx10g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/SortSam.jar \
INPUT=$WORKINGDIR/ReadsMap2Ref/"$SAMPLE_NAME"_alignment_$EBVgenome.sam \
OUTPUT=$WORKINGDIR/ReadsMap2Ref/"$SAMPLE_NAME"_alignment_$EBVgenome.sorted.bam \
SORT_ORDER=coordinate

#Index bam file
java -Xmx10g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/BuildBamIndex.jar \
INPUT=$WORKINGDIR/ReadsMap2Ref/"$SAMPLE_NAME"_alignment_$EBVgenome.sorted.bam

if [ -f $WORKINGDIR/ReadsMap2Ref/"$SAMPLE_NAME"_alignment_$EBVgenome.sorted.bam ]
then
echo "shuffled files is present..."
rm $WORKINGDIR/ReadsMap2Ref/"$SAMPLE_NAME"_alignment_$EBVgenome.sam
else
echo "sorted read2ref bam file not found"
fi


cd $WORKINGDIR/
else
echo "Skipping read Mapping"
fi


if [ "$SimpleStats" = "1" ]
then

cd $WORKINGDIR/ReadsMap2Ref/

java -Xmx16g -XX:ParallelGCThreads=$nt -jar \
$GATKdir/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
-o $WORKINGDIR/ReadsMap2Ref/"$SAMPLE_NAME"_alignment_$EBVgenome.sorted_DepthofCoverage \
-I $WORKINGDIR/ReadsMap2Ref/"$SAMPLE_NAME"_alignment_$EBVgenome.sorted.bam

#STATS
else
echo "Skipping Stats"
fi


done

Rscript ~/codes/EBVseq/EBVGenomeTyper_PlotCov.R $SAMPLE_NAME

