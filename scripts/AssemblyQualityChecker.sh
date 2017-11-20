#!/bin/bash

# Set the number of threads for tools that can be run in parallel.
#This is usually one digit below the number of cores this script is reserved to use (-n value above).
#Set the working directory and Sample name
dir=`pwd`
InputAssembly=$1
fastq_1=$2
fastq_2=$3
SAMPLE_NAME=$4

toolDir='/home/yk42w/codes/ViralGenomeAssembly'
PICARDPATH="$toolDir/bin/linux/picard-2.14.1/build/libs"
Bowtie2PATH="$toolDir/bin/linux/bowtie2-2.3.3.1-linux-x86_64"
INDEXGENOMESDIR="$toolDir/resources/Bowtie2Index"
nt=2
module load java/1.8.0_77

################################################
############# Main Control Panel: ##############
################################################
Reads2Assembly=0
FixAssemblyWithReads=1

##########################################
######### NO NEED TO MODIFY ! ############
############# Sub-routines ###############
##########################################
#Load required modules
module unload openssl/1.0.1g
module load ncbi_cxx/12_0_0
module load cutadapt/1.7.1
module load samtools/1.4.1
#module load bowtie2/2-2.1.0
module load fastx_toolkit/0.0.14
#module load java/1.7.0_25
module load fastqc/0.10.1
module load prinseq/0.20.4
module load gcc/4.8.1
module load bzip2/1.0.6
module load tabix/0.2.6
module load IGVTools/2.3.31
module load R/3.2.2
module load bedtools/2.17.0
module load vicuna/1.3
module load fastqc/0.10.1
module load perl/5.18.1
module load python/2.7.5

#Set the working environment
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/pkg/ncbi_cxx/12_0_0/lib
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/tagcleaner-standalone-0.16/:$PATH
export PERL5LIB=/project/umw_jeffrey_bailey/share/bin_sync/perl5lib/
#export PATH=/project/umw_jeffrey_bailey/share/bin_sync/Vfat:$PATH
#export PATH=/project/umw_jeffrey_bailey/share/bin_sync/picard-tools-1.105/:$PATH
export PERL5LIB=/project/umw_jeffrey_bailey/share/bin_sync/perl5lib/lib/perl5/
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/velvet_1.2.10:$PATH
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/VelvetOptimiser-2.2.5:$PATH
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/velvet_1.2.10/contrib/shuffleSequences_fasta:$PATH
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/MetaVelvet-1.2.02:$PATH
export PYTHONPATH=~/.local/lib/python2.7/site-packages/:$PATH

#Variables for required directories
BCFtoolsDIR="/project/umw_jeffrey_bailey/share/bin_sync/samtools-0.1.19/bcftools"
#PICARDPATH="/project/umw_jeffrey_bailey/share/bin_sync/picard-tools-1.105"
GATKdir="/project/umw_jeffrey_bailey/share/bin_sync/"
ShuffleDir='/project/umw_jeffrey_bailey/share/bin_sync/velvet_1.2.10/contrib/shuffleSequences_fasta'


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

#create bowtie2 index
$Bowtie2PATH/bowtie2-build $InputAssembly ${InputAssembly%.fa}_bowtieIndex
samtools faidx $InputAssembly

#Map reads back to assembly
$Bowtie2PATH/bowtie2 -p $nt \
-x ${InputAssembly%.fa}_bowtieIndex \
-1 $fastq_1 \
-2 $fastq_2 \
-S "$SAMPLE_NAME"_readsBack2assembly.sam

#Filter MQ<20 reads
samtools view -h -b -q 30 \
"$SAMPLE_NAME"_readsBack2assembly.sam \
-o "$SAMPLE_NAME"_readsBack2assembly_MQ30.bam

#Sort sam file and output as bam
java -Xmx10g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/picard.jar \
SortSam \
INPUT="$SAMPLE_NAME"_readsBack2assembly_MQ30.bam \
OUTPUT="$SAMPLE_NAME"_readsBack2assembly_MQ30_sorted.bam \
SORT_ORDER=coordinate


#Remove Duplicated reads
java -Xmx10g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/picard.jar \
MarkDuplicates \
I="$SAMPLE_NAME"_readsBack2assembly_MQ30_sorted.bam \
O="$SAMPLE_NAME"_readsBack2assembly_MQ30_DD.bam \
REMOVE_DUPLICATES=true \
METRICS_FILE=$SAMPLE_NAME._dedup_metrics.txt

#Fix pair information, sort, and index.
java -Xmx10g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/picard.jar \
FixMateInformation \
INPUT="$SAMPLE_NAME"_readsBack2assembly_MQ30_DD.bam \
OUTPUT="$SAMPLE_NAME"_readsBack2assembly_MQ30_DD_fixed.bam \
SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=LENIENT

#Index bam file
java -Xmx10g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/picard.jar \
BuildBamIndex \
INPUT="$SAMPLE_NAME"_readsBack2assembly_MQ30_DD_fixed.bam

samtools mpileup -uf $InputAssembly "$SAMPLE_NAME"_readsBack2assembly_MQ30_DD_fixed.bam > $SAMPLE_NAME.AssemblyRawcalls.bcf
/project/umw_jeffrey_bailey/share/bin_sync/samtools-0.1.19/bcftools/bcftools view -v $SAMPLE_NAME.AssemblyRawcalls.bcf > $SAMPLE_NAME.Assembly.variants.vcf

cat $SAMPLE_NAME.Assembly.variants.vcf| \
grep -v "#"|awk '($6 > 100)'|awk 'NF { info=$8; gsub(/.*;DP4=|;MQ=.*/, "", info); split(info, a, /,/); print $0 "\t" a[1]"\t"a[2]"\t"a[3]"\t"a[4]}' |\
awk '{if($NF+$(NF-1) != 0) print $1"_"$2"_"$4"_"$5 "\t" ($(NF)+$(NF-1)+$(NF-2)+$(NF-3)) "\t" ($NF+$(NF-1)) / ($(NF)+$(NF-1)+$(NF-2)+$(NF-3))}'|\
awk '($2 > 99)'|sed 's/,/_/g' > $SAMPLE_NAME.AssemblyMinorfreq.txt

Rscript $toolDir/bin/MAF_plot.R $SAMPLE_NAME $SAMPLE_NAME.AssemblyMinorfreq.txt

else
echo "no pileup"
fi

#*********************************************************************************************************#
if [ "$FixAssemblyWithReads" = "1" ]
then
module load python/2.7.5_packages/biopython/1.68
#fastaName=`head -1 ../JB_ContigMerge_v3/BL614_FNAPairedPlasma_PCRsWGA_EBV_type1_genome.fa|sed 's/>//g'`
#sed 's/'$EBVgenome'/'$fastaName'/g' ~/codes/ViralGenomeAssembly/resources/Annotation/Type1/NC_007605_repeatMask.bed > tmp.repeatfile.tmp

python $toolDir/bin/FixAssembly_With_Reads.py \
Reads2ConcensusGenome \
-n $EBVgenome \
--regions_to_mask $AnnotationDir/"$EBVgenome"_repeatMask.bed \
-r "$SAMPLE_NAME"_$RegionName.bam \
-f $INDEXGENOMESDIR/"$EBVgenome".fa \
-o "$SAMPLE_NAME"_EBV_"$Type"_$RegionName


else
echo "no read cons."
fi
