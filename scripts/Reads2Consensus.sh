#!/bin/bash

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
#submit this way: bsub -q long -n 6 -R rusage[mem=20000] -R "span[hosts=1]" -R "select[tmp>1000]" -W 96:00 -e err.%J.txt -o out.%J.txt ~/project/OTHERS/dnaSeq_Tools/scripts/EBV_SequenceAnalysisPipeline_Velvet_v2.sh 0 1
SAMPLE_NAME=`basename $dir`
Kmer_start=21
Kmer_end=151
nt=6
option='n50'
type=$2
fi


echo "job_id: $LSB_JOBID"
echo "index_id: $LSB_JOBINDEX"
echo "job_name: $LSB_JOBNAME"


################################################
############# Main Control Panel: ##############
################################################

#Setting the directories for required files
WORKINGDIR=$dir




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

GetConsensusFromReads=1


#*********************************************************************************************************#
if [ "$GetConsensusFromReads" = "1" ]
then
export PERL5LIB=/project/umw_jeffrey_bailey/share/bin_sync/trinityrnaseq_r20131110/PerlLib/
module load python/2.7.5_packages/biopython/1.68
#RegionName
RegionName=LMP1
#NC_009334:169179-170461
samtools view -b $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal.bam "$EBVgenome:169179-170461" > "$SAMPLE_NAME"_$RegionName.bam
samtools index "$SAMPLE_NAME"_$RegionName.bam

python /project/umw_jeffrey_bailey/OTHERS/dnaSeq_Tools/scripts/Reads2ConcensusGenome.py \
Reads2ConcensusGenome \
-n $EBVgenome \
--regions_to_mask $AnnotationDir/"$EBVgenome"_repeatMask.bed \
-r "$SAMPLE_NAME"_$RegionName.bam \
-f $INDEXGENOMESDIR/"$EBVgenome".fa \
-o "$SAMPLE_NAME"_EBV_"$Type"_$RegionName


else
echo "no read cons."
fi


