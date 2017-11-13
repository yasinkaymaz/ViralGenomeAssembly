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
#submit this way: bsub -q long -n 6 -R rusage[mem=20000] -R "span[hosts=1]" -R "select[tmp>1000]" -W 96:00 -e err.%J.txt -o out.%J.txt ~/project/OTHERS/dnaSeq_Tools/scripts/EBV_SequenceAnalysisPipeline.sh 0 1
SAMPLE_NAME=`basename $dir`
Kmer_start=21
Kmer_end=151
nt=10
option='n50'
type=$2
fi

echo "job_id: $LSB_JOBID"
echo "index_id: $LSB_JOBINDEX"
echo "job_name: $LSB_JOBNAME"

################################################
############# Main Control Panel: ##############
################################################
Reads2Assembly=1

#Setting the directories for required files
WORKINGDIR=$dir

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
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/Vfat:$PATH
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/picard-tools-1.105/:$PATH
export PERL5LIB=/project/umw_jeffrey_bailey/share/bin_sync/perl5lib/lib/perl5/
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/velvet_1.2.10:$PATH
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/VelvetOptimiser-2.2.5:$PATH
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/velvet_1.2.10/contrib/shuffleSequences_fasta:$PATH
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/MetaVelvet-1.2.02:$PATH
export PYTHONPATH=~/.local/lib/python2.7/site-packages/:$PATH

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


if [ "$Reads2Assembly" = "1" ]
then
rm -r $WORKINGDIR/Velvet_assembly/readsback2contigs
mkdir $WORKINGDIR/Velvet_assembly/readsback2contigs
cd $WORKINGDIR/Velvet_assembly/readsback2contigs

InputAssembly=`grep -w "$SAMPLE_NAME" $WORKINGDIR/../genomes4.txt| cut -f5`
#InputAssembly=$WORKINGDIR/Velvet_assembly/$JBmergerVersion/"$SAMPLE_NAME"_EBV_"$Type"_genome.fa
DATAdir='/project/umw_jeffrey_bailey/yk42w/results/DNAseq/EBV/Capture_seq/CombinedBatch'

fastq_1=$DATAdir/$SAMPLE_NAME/"$SAMPLE_NAME"_Hg_filt_pairs_1.fastq.gz
fastq_2=$DATAdir/$SAMPLE_NAME/"$SAMPLE_NAME"_Hg_filt_pairs_2.fastq.gz

#fastq_1=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_1.fastq
#fastq_2=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_2.fastq

bowtie2-build $WORKINGDIR/../$InputAssembly $WORKINGDIR/../${InputAssembly%.fa}_bowtieIndex
samtools faidx $WORKINGDIR/../$InputAssembly

bowtie2 -p $nt \
-x $WORKINGDIR/../${InputAssembly%.fa}_bowtieIndex \
-1 $fastq_1 \
-2 $fastq_2 \
-S "$SAMPLE_NAME"_readsBack2assembly.sam

samtools view -bS "$SAMPLE_NAME"_readsBack2assembly.sam > "$SAMPLE_NAME"_readsBack2assembly.bam
samtools sort "$SAMPLE_NAME"_readsBack2assembly.bam "$SAMPLE_NAME"_readsBack2assembly_sorted
samtools index "$SAMPLE_NAME"_readsBack2assembly_sorted.bam
#Should we filter low MQ reads ?

samtools mpileup -uf $WORKINGDIR/../$InputAssembly "$SAMPLE_NAME"_readsBack2assembly_sorted.bam > AssemblyRawcalls.bcf
/project/umw_jeffrey_bailey/share/bin_sync/samtools-0.1.19/bcftools/bcftools view -v AssemblyRawcalls.bcf > Assembly.variants.vcf

cat Assembly.variants.vcf| \
grep -v "#"|awk '($6 > 100)'|awk 'NF { info=$8; gsub(/.*;DP4=|;MQ=.*/, "", info); split(info, a, /,/); print $0 "\t" a[1]"\t"a[2]"\t"a[3]"\t"a[4]}' |\
awk '{if($NF+$(NF-1) != 0) print $1"_"$2"_"$4"_"$5 "\t" ($(NF)+$(NF-1)+$(NF-2)+$(NF-3)) "\t" ($NF+$(NF-1)) / ($(NF)+$(NF-1)+$(NF-2)+$(NF-3))}'|\
awk '($2 > 99)'|sed 's/,/_/g' > AssemblyMinorfreq.txt

Rscript ~/codes/EBVseq/MAF_plot.R $SAMPLE_NAME AssemblyMinorfreq.txt

else
echo "no pileup"
fi
