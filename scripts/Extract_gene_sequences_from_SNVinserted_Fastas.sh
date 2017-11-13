#!/bin/bash
#BSUB -n 2
#BSUB -R rusage[mem=10000]
#BSUB -R "span[hosts=1]"
#BSUB -q long
#BSUB -W 24:00
#BSUB -o "/project/umw_jeffrey_bailey/OTHERS/std_out/%J.out"
#BSUB -e "/project/umw_jeffrey_bailey/OTHERS/std_err/%J.err"
#BSUB -R select[tmp>1000]

# Set the number of threads for tools that can be run in parallel.
#This is usually one digit below the number of cores this script is reserved to use (-n value above).
nt=2

#Set the working directory and Sample name
dir=`pwd`
SAMPLE_NAME=`basename $dir`
$LSB_JOBNAME=$SAMPLE_NAME
echo "job_id: $LSB_JOBID"
echo "index_id: $LSB_JOBINDEX"
echo "job_name: $LSB_JOBNAME"


################################################
############# Main Control Panel: ##############
################################################

#Setting the directories for required files
WORKINGDIR=$dir

#Choose the type of the EBV genome!
#type=1
type=2
ExtractGeneSequences=1
ExtractInfo=1
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

if [ "$ExtractGeneSequences" = "1" ]
then 

export PATH=/project/umw_jeffrey_bailey/share/bin_sync/wise2.2.0/wisecfg/:$PATH

if [ -f vfat_annotation_aligned.afa ]
then

VFATOut='vfat_annotation_afa'
perl $VfatDIR/annotate.pl \
-fa "$SAMPLE_NAME"_AlternativeGenome_with_gatk_filterPassed_SNVs.fasta \
-ref $INDEXGENOMESDIR/$EBVgenome.fa \
-genelist $RefDIR/EBV_Reference_genelist.txt \
-pepfolder $RefDIR/EBV_Peptides \
-align vfat_annotation_aligned.afa \
-o $VFATOut

else
VFATOut='vfat_annotation'
perl $VfatDIR/annotate.pl \
-fa "$SAMPLE_NAME"_AlternativeGenome_with_gatk_filterPassed_SNVs.fasta \
-ref $INDEXGENOMESDIR/$EBVgenome.fa \
-genelist $RefDIR/EBV_Reference_genelist.txt \
-pepfolder $RefDIR/EBV_Peptides \
-o $VFATOut

fi

else
echo "No gene extraction"
fi



if [ "$ExtractInfo" = "1" ]
then

export PATH=/project/umw_jeffrey_bailey/share/bin_sync/wise2.2.0/wisecfg/:$PATH

AlternateFastaFile="$SAMPLE_NAME"_AlternativeGenome_with_gatk_filterPassed_SNPs.fasta

#Choose type
if [ "$type" = "1" ]
then
VfatDIR='/project/umw_jeffrey_bailey/share/bin_sync/Vfat/'
RefDIR='/project/umw_jeffrey_bailey/share/EBV/Annotation/Type1/Vfat/EBV'
else
VfatDIR='/project/umw_jeffrey_bailey/share/bin_sync/Vfat/temp/'
RefDIR='/project/umw_jeffrey_bailey/share/EBV/Annotation/Type2/Vfat/EBV'
fi

#Vfat annotate
#perl $VfatDIR/annotate.pl \
#-fa $AlternateFastaFile \
#-ref $RefDIR/EBV_Reference.fasta \
#-genelist $RefDIR/EBV_Reference_genelist.txt \
#-pepfolder $RefDIR/EBV_Peptides/ \
#-align ./ExtractInfoOut_aligned.afa \
#-o ExtractInfoOut
#cd ExtractInfoOut


cd $VFATOut
#Parse GeneWise DNA output
for file in `ls -1|grep DNA`;\
do gene=${file%_genewiseDNA.txt}; \
sed 's,//,*,p' $file |\
sed -n '/].sp/,/*/p'|\
sed 's,*,,p'|\
sed ':a;N;$!ba;s/\n\n//g'|\
sed "s,>1,>$SAMPLE_NAME,p"|\
uniq|sed "s,sp,$gene.dna,p"|\
uniq|awk '!/^>/ { printf "%s", $0; n = "\n" }/^>/{ print n $0; n = "" }END{ printf "%s", n }'|\
sed ':a;N;$!ba;s/.dna\n/.dna\t/g'|\
awk '{print $1"\t"$2"\t"length($2)}'|\
sort -k3,3gr|head -1|\
awk '{print $1"\n"$2}' > $gene.dna.fa; \
done

cd $WORKINGDIR/

else
echo "Skipping Extract info!!!"
fi

WORKINGDIR=/home/yk42w/project/yk42w/results/DNAseq/EBV/Capture_seq/CombinedBatch/Hg_clean
for i in `cut -f1 samples.2 |grep -v WTSI|tail -n62`;
do
cd $i;
echo $i;
dir=`pwd`;
SAMPLE_NAME=$dir;
cd $dir/Velvet_assembly/GeneSequences/vfat_annotation;
for file in `ls -1|grep DNA`;\
do gene=${file%_genewiseDNA.txt}; \
echo $gene;
sed 's,//,*,p' $file |\
sed -n '/].sp/,/*/p'|\
sed 's,*,,p'|\
sed ':a;N;$!ba;s/\n\n//g'|\
sed "s,>1,>$SAMPLE_NAME,p"|\
uniq|sed "s,sp,$gene.dna,p"|\
uniq|awk '!/^>/ { printf "%s", $0; n = "\n" }/^>/{ print n $0; n = "" }END{ printf "%s", n }'|\
sed ':a;N;$!ba;s/.dna\n/.dna\t/g'|\
awk '{print $1"\t"$2"\t"length($2)}'|\
sort -k3,3gr|head -1|\
awk '{print $1"\n"$2}' > $gene.dna.fa; \
done
cd $WORKINGDIR;
done

