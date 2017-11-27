#!/bin/bash


dir=`pwd`
InputAssembly=$1
SAMPLE_NAME=$2
type=$3

toolDir='/home/yk42w/codes/ViralGenomeAssembly'
PICARDPATH="$toolDir/bin/linux/picard-2.14.1/build/libs"
Bowtie2PATH="$toolDir/bin/linux/bowtie2-2.3.3.1-linux-x86_64"
INDEXGENOMESDIR="$toolDir/resources/Bowtie2Index"
nt=4
module load java/1.8.0_77


ExtractGeneSequences=1

module unload openssl/1.0.1g
module load samtools/1.4.1
#module load bowtie2/2-2.1.0
module load fastx_toolkit/0.0.14
#module load java/1.7.0_25
module load R/3.2.2
module load bedtools/2.17.0
module load perl/5.18.1
module load python/2.7.5

#Set the working environment
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/pkg/ncbi_cxx/12_0_0/lib
# export PATH=/project/umw_jeffrey_bailey/share/bin_sync/tagcleaner-standalone-0.16/:$PATH
export PERL5LIB=/project/umw_jeffrey_bailey/share/bin_sync/perl5lib/
#export PATH=/project/umw_jeffrey_bailey/share/bin_sync/Vfat:$PATH
#export PATH=/project/umw_jeffrey_bailey/share/bin_sync/picard-tools-1.105/:$PATH
export PERL5LIB=/project/umw_jeffrey_bailey/share/bin_sync/perl5lib/lib/perl5/
# export PATH=/project/umw_jeffrey_bailey/share/bin_sync/velvet_1.2.10:$PATH
# export PATH=/project/umw_jeffrey_bailey/share/bin_sync/VelvetOptimiser-2.2.5:$PATH
# export PATH=/project/umw_jeffrey_bailey/share/bin_sync/velvet_1.2.10/contrib/shuffleSequences_fasta:$PATH
# export PATH=/project/umw_jeffrey_bailey/share/bin_sync/MetaVelvet-1.2.02:$PATH
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


if [ "$ExtractGeneSequences" = "1" ]
then

perl $VfatDIR/annotate.pl \
-fa $InputAssembly \
-ref $toolDir/resources/Bowtie2Index/$EBVgenome.fa \
-genelist $RefDIR/EBV_Reference_genelist.txt \
-pepfolder $RefDIR/EBV_Peptides \
-o "$SAMPLE_NAME"_vfat_annotation

else
echo "No gene extraction"
fi
