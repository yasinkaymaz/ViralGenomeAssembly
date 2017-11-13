#!/bin/bash

# Set the number of threads for tools that can be run in parallel.
#This is usually one digit below the number of cores this script is reserved to use (-n value above).


#Set the working directory and Sample name
dir=`pwd`

Batch=$1

if [ "$Batch" = "1" ]
then
SAMPLE_NAME=$2
nt=$3
type=$4
ploidy=1
else
#submit this way: bsub -q long -n 6 -R rusage[mem=20000] -R "span[hosts=1]" -R "select[tmp>1000]" -W 72:00 -e err.%J.txt -o out.%J.txt ~/project/OTHERS/dnaSeq_Tools/scripts/EBV_SequenceAnalysisPipeline_SNPAnalysis.sh 0 1
SAMPLE_NAME=`basename $dir`
nt=12
type=$2
ploidy=$3
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


PreProcessForGATK=0

# Recalibration of bases after realignment
Recalibration=0

#Calling variations
VariantCalling=0
#Parsing variations
MergeAndCall=1
MOI_analysis=0

SNVFasta2Genes=0
SNPdat=0

CleanIntermediateFiles=0



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



### Pre-process alignment file ###
#---------------------------------#
if [ "$PreProcessForGATK" = "1" ]
then

mkdir $WORKINGDIR/GATK
cd $WORKINGDIR/GATK

ln -s $WORKINGDIR/ReadsMap2Ref/"$SAMPLE_NAME"_alignment_$EBVgenome.sorted.bam ./

#MARK DUPLICATES
java -Xmx10g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/MarkDuplicates.jar \
INPUT=$WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.sorted.bam \
OUTPUT=$WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.sorted_DD.bam \
REMOVE_DUPLICATES=true \
METRICS_FILE=$WORKINGDIR/GATK/"$SAMPLE_NAME"_dedup_metrics.txt

#Run Picard Tools
java -Xmx10g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/AddOrReplaceReadGroups.jar \
I=$WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.sorted_DD.bam \
O=$WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.sorted_RGadded.bam \
SORT_ORDER=coordinate \
CREATE_INDEX=true \
RGPL=ILLUMINA \
RGID=$LSB_JOBID \
RGSM=$SAMPLE_NAME \
RGLB=Agilent \
RGPU=Run_Barcode \
VALIDATION_STRINGENCY=LENIENT

#Fix Mate pair
java -Xmx31g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/FixMateInformation.jar \
INPUT=$WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.sorted_RGadded.bam \
OUTPUT=$WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.sorted_RGadded_fixed.bam \
SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=LENIENT

#Index bam file
java -Xmx10g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/BuildBamIndex.jar \
INPUT=$WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.sorted_RGadded_fixed.bam

cd $WORKINGDIR/
else
echo "Skipping data pre-processing step"
fi

################# GATK ####################
#-----------------------------------------#
if [ "$Recalibration" = "1" ]
then
#Perform local realignment around indels to correct mapping-related artifacts.
cd $WORKINGDIR/GATK/

java -Xmx16g -XX:ParallelGCThreads=$nt -jar \
$GATKdir/GenomeAnalysisTK.jar \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
-I $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.sorted_RGadded_fixed.bam \
-T RealignerTargetCreator \
-nt $nt \
-o $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_target.intervals

java -Xmx16g -XX:ParallelGCThreads=$nt -jar \
$GATKdir/GenomeAnalysisTK.jar \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
-I $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.sorted_RGadded_fixed.bam \
-T IndelRealigner \
-targetIntervals $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_target.intervals \
-o $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_realigned.bam

#Analyze patterns of covariation in the sequence dataset.
java -Xmx16g -XX:ParallelGCThreads=$nt -jar $GATKdir/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
-I $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_realigned.bam \
-cov ReadGroupcovariate \
-cov QualityScoreCovariate \
-cov CycleCovariate \
-cov ContextCovariate \
--run_without_dbsnp_potentially_ruining_quality \
-o $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal_data.table

#Do a second pass to analyze covariation remaining after recalibration.
java -Xmx16g -XX:ParallelGCThreads=$nt -jar $GATKdir/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
-I $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_realigned.bam \
--run_without_dbsnp_potentially_ruining_quality \
-BQSR $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal_data.table \
-o $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_post_recal_data.table

#Generate before/after plots.
java -Xmx16g -XX:ParallelGCThreads=$nt -jar $GATKdir/GenomeAnalysisTK.jar \
-T AnalyzeCovariates \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
-before $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal_data.table \
-after $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_post_recal_data.table \
-plots $WORKINGDIR/GATK/"$SAMPLE_NAME"_recalibration_plots.pdf

#Apply the recalibration to your sequence data.
java -Xmx16g -XX:ParallelGCThreads=$nt -jar $GATKdir/GenomeAnalysisTK.jar \
-T PrintReads \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
-I $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_realigned.bam \
-nct $nt \
-BQSR $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal_data.table \
--read_filter MappingQualityZero \
-o $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal.bam

samtools index $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal.bam

java -Xmx16g -XX:ParallelGCThreads=$nt -jar $GATKdir/GenomeAnalysisTK.jar -T DepthOfCoverage -R $INDEXGENOMESDIR/$EBVgenome.fa -o $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal_DepthofCoverage -I $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal.bam -geneList $RefDIR/EBV_Reference_genelist.txt --minMappingQuality 20

cd $WORKINGDIR/
else
echo "Skipping Recalibration step"
fi

### Variant Discovery ###
#Variant calling
if [ "$VariantCalling" = "1" ]
then

cd $WORKINGDIR/GATK/

ploidy=2

java -Xmx16g -XX:ParallelGCThreads=$nt -jar \
$GATKdir/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
-I $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal.bam \
-nct $nt \
-ploidy $ploidy \
-glm BOTH \
-stand_call_conf 30.0 \
-stand_emit_conf 10.0 \
-o $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk.raw.vcf

#Skipping VQSR (Variant Quality Score Recalibration) step.
#Apply hard filters to the call set (cause the purpose of this pass is to build a truth/training set for next runs).

#Extract the SNPs from the call set
java -Xmx16g -XX:ParallelGCThreads=$nt -jar \
$GATKdir/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
-V $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk.raw.vcf \
-selectType SNP \
-o $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_snp_raw.vcf 

#Apply the filter to the SNP call set
java -Xmx16g -XX:ParallelGCThreads=$nt -jar \
$GATKdir/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
-V:VCF $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_snp_raw.vcf \
--filterExpression "QD < 2.0" \
--filterName "QDFilter" \
--filterExpression "MQ < 10.0" \
--filterName "MQFilter" \
--filterExpression "FS > 60.0" \
--filterName "FSFilter" \
--filterExpression "MQRankSum < -12.5" \
--filterName "MQRankSumFilter" \
--filterExpression "ReadPosRankSum < -8.0" \
--filterName "ReadPosRankSumFilter" \
--missingValuesInExpressionsShouldEvaluateAsFailing \
-o $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_snp_filtered.vcf
#DP filter can be included for whole genome sequences.

#Extract the Indels from the call set
java -Xmx16g -XX:ParallelGCThreads=$nt -jar \
$GATKdir/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
--variant $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk.raw.vcf \
-selectType INDEL \
-o $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_indel_raw.vcf 

#Apply the filter to the Indel call set
java -Xmx16g -XX:ParallelGCThreads=$nt -jar \
$GATKdir/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
-V:VCF $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_indel_raw.vcf \
--filterExpression "QD < 2.0" \
--filterName "QDFilter" \
--filterExpression "ReadPosRankSum < -20.0" \
--filterName "ReadPosRankSumFilter" \
--filterExpression "FS > 200.0" \
--filterName "FSFilter" \
--missingValuesInExpressionsShouldEvaluateAsFailing \
-o $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_indel_filtered.vcf
#!!InbreedingCoeff < -0.8	#Add this parameter when running more than 10 samples.
#DP filter can be included for whole genome sequences.

java -Xmx16g -XX:ParallelGCThreads=$nt -jar \
$GATKdir/GenomeAnalysisTK.jar \
-T CombineVariants \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
-V:VCF $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_indel_filtered.vcf \
-V:VCF $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_snp_filtered.vcf \
-o $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_V3.vcf

grep "#" "$SAMPLE_NAME"_alignment_"$EBVgenome".gatk_snp_filtered.vcf > "$SAMPLE_NAME"_alignment_"$EBVgenome".gatk_snp_filtered_PASSed.vcf
grep PASS "$SAMPLE_NAME"_alignment_"$EBVgenome".gatk_snp_filtered.vcf >> "$SAMPLE_NAME"_alignment_"$EBVgenome".gatk_snp_filtered_PASSed.vcf
gzip "$SAMPLE_NAME"_alignment_"$EBVgenome".gatk_snp_filtered_PASSed.vcf


cd $WORKINGDIR/
else
echo "Skipping Variant calling"
fi


if [ "$MOI_analysis" = "1" ]
then

mkdir $WORKINGDIR/estMOI
cd $WORKINGDIR/estMOI

ln -s $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal.ba* ./
ln -s $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_"$EBVgenome".gatk_snp_filtered.vcf* ./
ln -s $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_"$EBVgenome".gatk_snp_filtered_PASSed.vcf.gz ./

Outname=PASSed_estMOI_10Percent
minHap=10
$toolDir/scripts/estMOI.pl \
"$SAMPLE_NAME"_alignment_"$EBVgenome".gatk_recal.bam \
"$SAMPLE_NAME"_alignment_"$EBVgenome".gatk_snp_filtered_PASSed.vcf.gz \
$INDEXGENOMESDIR/$EBVgenome.fa \
-exclude $toolDir/resources/estMOI/"$EBVgenome"_estMOI_exclude.bed \
-out $Outname \
--minhap=$minHap

#sed 's, ,_,g'|sed 's,_NC_,NC,g'|grep -v _N|\
grep Hap $Outname.moi.3.10.500.$minHap.log |\
sed 's, ,_,g'|sed 's,_NC_,NC,g'|\
groupBy -g 2 -c 3,4,4 -o collapse,collapse,sum |\
sed 's,#_,,g'|\
expandCols -c 2,3|\
awk '{print $1"\t"$2"\t"$3"\t"($3*100/$4)"\t""1"}'|awk '($4>1.0)'|\
groupBy -g 1 -c 2,3,4,5 -o collapse,collapse,collapse,sum |\
sed 's,alignment,,g'|\
sed 's,gatk_recal.bam,,g'|\
sed 's,Hapotype,Haplotype,g' > $Outname.moi.3.10.500.$minHap.frequencies

##Column 5 is major haplo freq, 6 is minor.
bedtools expand -c 3 -i $Outname.moi.3.10.500.$minHap.frequencies |\
groupBy -g 1 -c 2,3,4,5 -o first,sum,first,first|\
awk '($3 > 99)'|awk '($5 < 3)'|sed 's/,/\t/g'|\
awk '{OFS="\t";if($NF == 1) print $1,$2,$2,$3,$4,$4,$5; else print $0}'|\
awk '{OFS="\t"; if($6 > $5) print $1,$2,$3,$4,$5,$7,$6,$8; else print $0}' > haplo.freq.txt
#bedtools expand -c 3 -i PASSed_estMOI_5Prct.moi.3.10.500.3.frequencies|groupBy -g 1 -c 2,3,4,5 -o first,sum,first,first > ../../haploFreqs/multi/freq.txt
echo $SAMPLE_NAME > $Outname.moi.3.10.500.$minHap.frequencies.hist
echo 'MOI	#ofHaplotypes	Frequencies' >> $Outname.moi.3.10.500.$minHap.frequencies.hist
cat $Outname.moi.3.10.500.$minHap.frequencies |\
cut -f5|sort -n|uniq -c|\
awk '{print "MOI:""\t"$2"\t"$1}'|\
groupBy -g 1 -c 2,3,3 -o collapse,collapse,sum|\
expandCols -c 2,3 |\
awk '{print $2"\t"$3"\t"($3*100/$4)}' >> $Outname.moi.3.10.500.$minHap.frequencies.hist

cd $WORKINGDIR/
else
echo "Skipping MOI analysis"
fi



if [ "$SNPdat" = "1" ]
then

mkdir $WORKINGDIR/GATK/SNPdat
cd $WORKINGDIR/GATK/SNPdat

ln -s $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_"$EBVgenome".gatk_snp_filtered_PASSed.vcf.gz ./
gunzip "$SAMPLE_NAME"_alignment_"$EBVgenome".gatk_snp_filtered_PASSed.vcf.gz

perl /project/umw_jeffrey_bailey/share/bin_sync/SNPdat_package_v1.0.5/SNPdat_v1.0.5.pl \
-i "$SAMPLE_NAME"_alignment_"$EBVgenome".gatk_snp_filtered_PASSed.vcf \
-g $toolDir/resources/Annotation/Type"$type"/custom_T"$type".gtf \
-f $toolDir/resources/Bowtie2Index/"$EBVgenome".fa \
-o "$SAMPLE_NAME"_SNPdat

else
echo "Skip SNPdat"
fi


if [ "$MergeAndCall" = "1" ]
then

module unload samtools/0.0.19
module load samtools/1.3

cd $WORKINGDIR/

samtools merge -fr "$SAMPLE_NAME"_$EBVgenome.merged.bam *recal.bam
samtools sort "$SAMPLE_NAME"_$EBVgenome.merged.bam -o "$SAMPLE_NAME"_$EBVgenome.merged_sorted.bam
samtools index "$SAMPLE_NAME"_$EBVgenome.merged_sorted.bam

java -Xmx16g -XX:ParallelGCThreads=$nt -jar \
$GATKdir/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
-I "$SAMPLE_NAME"_$EBVgenome.merged_sorted.bam \
-nct $nt \
-ploidy $ploidy \
-glm BOTH \
-stand_call_conf 30.0 \
-stand_emit_conf 10.0 \
-o "$SAMPLE_NAME"_$EBVgenome.merged_gatk.raw.vcf


#Extract the SNPs from the call set
java -Xmx16g -XX:ParallelGCThreads=$nt -jar \
$GATKdir/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
-V "$SAMPLE_NAME"_$EBVgenome.merged_gatk.raw.vcf \
-selectType SNP \
-o "$SAMPLE_NAME"_$EBVgenome.merged_gatk_snp_raw.vcf 

#Apply the filter to the SNP call set
java -Xmx16g -XX:ParallelGCThreads=$nt -jar \
$GATKdir/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
-V:VCF "$SAMPLE_NAME"_$EBVgenome.merged_gatk_snp_raw.vcf \
--filterExpression "QD < 2.0" \
--filterName "QDFilter" \
--filterExpression "MQ < 10.0" \
--filterName "MQFilter" \
--filterExpression "FS > 60.0" \
--filterName "FSFilter" \
--filterExpression "MQRankSum < -12.5" \
--filterName "MQRankSumFilter" \
--filterExpression "ReadPosRankSum < -8.0" \
--filterName "ReadPosRankSumFilter" \
--missingValuesInExpressionsShouldEvaluateAsFailing \
-o "$SAMPLE_NAME"_$EBVgenome.merged_gatk_snp_filtered.vcf
#DP filter can be included for whole genome sequences.

#Extract the Indels from the call set
java -Xmx16g -XX:ParallelGCThreads=$nt -jar \
$GATKdir/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
--variant "$SAMPLE_NAME"_$EBVgenome.merged_gatk.raw.vcf \
-selectType INDEL \
-o "$SAMPLE_NAME"_$EBVgenome.merged_gatk_indel_raw.vcf 

#Apply the filter to the Indel call set
java -Xmx16g -XX:ParallelGCThreads=$nt -jar \
$GATKdir/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
-V:VCF "$SAMPLE_NAME"_$EBVgenome.merged_gatk_indel_raw.vcf \
--filterExpression "QD < 2.0" \
--filterName "QDFilter" \
--filterExpression "ReadPosRankSum < -20.0" \
--filterName "ReadPosRankSumFilter" \
--filterExpression "FS > 200.0" \
--filterName "FSFilter" \
--missingValuesInExpressionsShouldEvaluateAsFailing \
-o "$SAMPLE_NAME"_$EBVgenome.merged_gatk_indel_filtered.vcf
#!!InbreedingCoeff < -0.8	#Add this parameter when running more than 10 samples.
#DP filter can be included for whole genome sequences.

java -Xmx16g -XX:ParallelGCThreads=$nt -jar \
$GATKdir/GenomeAnalysisTK.jar \
-T CombineVariants \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
-V:VCF "$SAMPLE_NAME"_$EBVgenome.merged_gatk_indel_filtered.vcf \
-V:VCF "$SAMPLE_NAME"_$EBVgenome.merged_gatk_snp_filtered.vcf \
-o "$SAMPLE_NAME"_$EBVgenome.merged_gatk_V3.vcf

grep "#" "$SAMPLE_NAME"_$EBVgenome.merged_gatk_V3.vcf > "$SAMPLE_NAME"_$EBVgenome.merged_gatk_PASSed.vcf
grep PASS "$SAMPLE_NAME"_$EBVgenome.merged_gatk_V3.vcf >> "$SAMPLE_NAME"_$EBVgenome.merged_gatk_PASSed.vcf
#gzip "$SAMPLE_NAME"_alignment_"$EBVgenome".gatk_snp_filtered_PASSed.vcf



else
echo "no merge and call"
fi


if [ "$CleanIntermediateFiles" = "1" ]
then


#Check if the output is generated correctly then, remove intermediates:
s=`ls -1 $WORKINGDIR/GATK/|grep "$EBVgenome".sorted_RGadded_fixed.bam |wc -l`
if [ $s = "1" ]
then
rm $WORKINGDIR/ReadsMap2Ref/*alignment_*.sam
rm $WORKINGDIR/GATK/*alignment_$EBVgenome.sorted_DD.ba*
rm $WORKINGDIR/GATK/*alignment_$EBVgenome.sorted_RGadded.ba*
rm $WORKINGDIR/GATK/*alignment_$EBVgenome.sorted_RGadded_fixed.ba*
rm $WORKINGDIR/GATK/*.gatk_target.intervals
rm $WORKINGDIR/GATK/*.gatk_recal_data.table
rm $WORKINGDIR/GATK/*.gatk_post_recal_data.table

else
echo "Proper output files not found"
fi
####CleanUp####

else
echo "NO need to clean"
fi


