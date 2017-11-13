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

#ls -1 *R1_00*|sort|xargs -i zcat {} >> "$SAMPLE_NAME"_R1.fastq
#ls -1 *R2_00*|sort|xargs -i zcat {} >> "$SAMPLE_NAME"_R2.fastq

#gzip $WORKINGDIR/"$SAMPLE_NAME"_R*.fastq

fastq_1=$WORKINGDIR/"$SAMPLE_NAME"_R1.fastq.gz
fastq_2=$WORKINGDIR/"$SAMPLE_NAME"_R2.fastq.gz

# READ Pre-PROCESS: This function will automatically
# 1) Trims adapter sequeces from both reads,
# 2) Excludes low complexity reads (poly As, Ts etc.),
# 3) Trims low quality based off of reads starting from 3'ends,
# 4) After each step of precessing, checks the quality/length distribution of reads.
### Please Turn this on if you need by changing '0' to '1' 
ReadPreProc=0
FilterHumanReads=0

#Denovo assembly with velvet
#Compiling velvet with specific settings
#compile=0
#velvet=0
#PostVelvet=0
Dedupreads=0
velvetOptimizer=0
velvet=0
metavelvet=0
#velvet2=0
PostvelvetOptimizer=0

#Post de-novo assembly process
ReadMapping2Ref=0
AlignStats=0
GetConsensusFromReads=0
BamDownSample=1
JB_Merger=1
ContigsReadPileup=0
AssemblyContigsReadPileup=0
Vfat_orient_merge=0
ICORN=0

Fix_FrameShifts=0
AssemblyContigsChecker=0
ExtractGeneSequences=0
ExtractInfo=0

SNPdat=0

AA_Change_Rates=0
SimpleStats=0

CircosPlotFiles=0
Vfat=0
# Remove intermediate processing files to save space. Keep this off for de-bugging purposes.
### Please Turn this on if you need by changing '0' to '1'
CleanIntermediateFiles=0

#array=(Lbp Lcon max tbp n50)
#array=(n50 Lcon Lbp)
array=(n50)

#array=(Lbp)
#LNbp = The total number of Ns in large contigs
#Lbp = The total number of base pairs in large contigs
#Lcon = The number of large contigs
#max = The length of the longest contig
#n50 = The n50
#ncon = The total number of contigs
#tbp = The total number of basepairs in contigs
#velvetsettings='AutoInsert_Min150_Kmer21t75';
velvetsettings='Velvet_51_151_kN50_cLbp';
# 21 to 53 and -cov_cutoff 5

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

if [ "$FilterHumanReads" = "1" ]
then

cd $WORKINGDIR/
fastq_1=$WORKINGDIR/"$SAMPLE_NAME"_goodQual_1.fastq
fastq_2=$WORKINGDIR/"$SAMPLE_NAME"_goodQual_2.fastq

#Align reads to Hg19 and keep unmapped reads.
BOWTIE2hg19DIR="/project/umw_jeffrey_bailey/share/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index"

bowtie2 -p $nt -x $BOWTIE2hg19DIR/genome \
--un-conc "$SAMPLE_NAME"_Hg_filt_pairs \
--un "$SAMPLE_NAME"_Hg_filt_single \
-1 $fastq_1 -2 $fastq_2 -S "$SAMPLE_NAME"_Hg_mapped.bam

#samtools view -bS "$SAMPLE_NAME"_Hg_mapped.sam > "$SAMPLE_NAME"_Hg_mapped.bam
#samtools sort "$SAMPLE_NAME"_Hg_mapped.bam "$SAMPLE_NAME"_Hg_mapped_sorted
#samtools index "$SAMPLE_NAME"_Hg_mapped_sorted.bam

#fastqc "$SAMPLE_NAME"_Hg_filt_pairs*

mv "$SAMPLE_NAME"_Hg_filt_pairs.1 "$SAMPLE_NAME"_Hg_filt_pairs_1.fastq
mv "$SAMPLE_NAME"_Hg_filt_pairs.2 "$SAMPLE_NAME"_Hg_filt_pairs_2.fastq

rm $WORKINGDIR/"$SAMPLE_NAME"_goodQual_*.fastq
rm $WORKINGDIR/"$SAMPLE_NAME"_AdpTrimmed_*.fastq
rm $WORKINGDIR/"$SAMPLE_NAME"_Hg_mapped.bam
rm $WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_single*

else
echo "no human filtration..."
fi

#########################################################
######  Alternative De Novo assembly with VelVet  #######
#########################################################

#*********************************************************************************************************#
if [ "$velvetOptimizer" = "1" ]
then
#Alternatively, Finding the best kmer size
#/project/umw_jeffrey_bailey/OTHERS/dnaSeq_Tools/scripts/velvetk.pl --genome $INDEXGENOMESDIR/$EBVgenome.fa EBVBL_40_wga_Hg_filt_pairs_1.fastq EBVBL_40_wga_Hg_filt_pairs_2.fastq

echo "optimizing assembly for $option";


mkdir $WORKINGDIR/Velvet_$velvetsettings
cd $WORKINGDIR/Velvet_$velvetsettings

#gunzip $WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_*.fastq.gz

fastq_1=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_1.fastq
fastq_2=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_2.fastq


if [ "$Dedupreads" = 1 ]
then
echo $fastq_1 > infiles
echo $fastq_2 >> infiles
/home/yk42w/project/share/bin_sync/FastUniq/source/fastuniq -i infiles -t q -o $WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_Deduped_1.fastq -p $WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_Deduped_2.fastq
fastq_1=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_Deduped_1.fastq
fastq_2=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_Deduped_2.fastq
else
echo "Not dedup reads"
fi

fastq_1=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_Deduped_1.fastq
fastq_2=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_Deduped_2.fastq

READCOUNT=`awk '{if(NR%4==1)print}' $fastq_1|wc -l`;
if (( READCOUNT > 10000000 ));
then
#subsample Reads.
python /home/yk42w/codes/RandomSampler.py $fastq_1 $fastq_2 $WORKINGDIR/"$SAMPLE_NAME"_sub_1.fastq $WORKINGDIR/"$SAMPLE_NAME"_sub_2.fastq 10000000
fastq_1=$WORKINGDIR/"$SAMPLE_NAME"_sub_1.fastq
fastq_2=$WORKINGDIR/"$SAMPLE_NAME"_sub_2.fastq
else
echo "Not ultra-deep sequenced."
fi

if [ -f $WORKINGDIR/Velvet_$velvetsettings/"$SAMPLE_NAME"_shuffled.fastq ]
then
echo "shuffled files is present..."
ln -s $WORKINGDIR/Velvet_$velvetsettings/"$SAMPLE_NAME"_shuffled.fastq ./
else
perl $ShuffleDir/shuffleSequences_fastq.pl $fastq_1 $fastq_2 "$SAMPLE_NAME"_shuffled.fastq
fi


VelvetHString='-create_binary -shortPaired -fastq '$SAMPLE_NAME'_shuffled.fastq';
VelvetGString='-ins_length 300 -min_contig_lgth 200'
#VelvetGString='-cov_cutoff 5 -ins_length 400 -min_contig_lgth 200'
#VelvetGString='-min_contig_lgth 200'
#VelvetGString='-min_contig_lgth 150'

VelvetOptimiser.pl -s $Kmer_start -e $Kmer_end \
-f "$VelvetHString" \
-o "$VelvetGString" \
-t $nt \
-k $option \
-c "Lbp" \
-d ""$option"_outDir";

cd "$option"_outDir/;

perl $toolDir/scripts/contigs_stats.pl -t Velvet \
contigs.fa -plot > contigs_stats.txt;

rm CnyUnifiedSeq*
rm *Graph*

bowtie2 -p $nt -f -x $INDEXGENOMESDIR/$EBVgenome -U contigs.fa -S "$SAMPLE_NAME"_velvet_"$option".sam;
samtools view -bS "$SAMPLE_NAME"_velvet_"$option".sam > "$SAMPLE_NAME"_velvet_"$option".bam;
samtools sort "$SAMPLE_NAME"_velvet_"$option".bam "$SAMPLE_NAME"_velvet_"$option"_sorted;
samtools index "$SAMPLE_NAME"_velvet_"$option"_sorted.bam;

cd $WORKINGDIR

else
echo "Skipping velvetoptimiser";
fi
#*********************************************************************************************************#

if [ "$velvet" = "1" ]
then

mkdir $WORKINGDIR/Velvet_assembly
cd $WORKINGDIR/Velvet_assembly

#gunzip $WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_*.fastq.gz

fastq_1=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_1.fastq
fastq_2=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_2.fastq

if [ "$Dedupreads" = 1 ]
then
echo $fastq_1 > infiles
echo $fastq_2 >> infiles
/home/yk42w/project/share/bin_sync/FastUniq/source/fastuniq -i infiles -t q -o $WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_Deduped_1.fastq -p $WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_Deduped_2.fastq
fastq_1=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_Deduped_1.fastq
fastq_2=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_Deduped_2.fastq
else
echo "Not dedup reads"
fi

fastq_1=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_Deduped_1.fastq
fastq_2=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_Deduped_2.fastq

#READCOUNT=`awk '{if(NR%4==1)print}' $fastq_1|wc -l`;
#if (( READCOUNT > 10000000 ));
#then
#subsample Reads.
#python /home/yk42w/codes/RandomSampler.py $fastq_1 $fastq_2 $WORKINGDIR/"$SAMPLE_NAME"_sub_1.fastq $WORKINGDIR/"$SAMPLE_NAME"_sub_2.fastq 3000000
fastq_1=$WORKINGDIR/"$SAMPLE_NAME"_sub_1.fastq
fastq_2=$WORKINGDIR/"$SAMPLE_NAME"_sub_2.fastq
#else
#echo "Not ultra-deep sequenced."
#fi

export OMP_THREAD_LIMIT=$nt
export OMP_NUM_THREADS=$nt
velveth $WORKINGDIR/Velvet_assembly/ $Kmer_start,$Kmer_end,2 \
-create_binary \
-fastq \
-shortPaired \
-separate $fastq_1 $fastq_2

let end=$Kmer_end-2;

#Run velvet graph for each kmer length separately. 
for((n=$Kmer_start; n<=$end; n=n+2));
do
velvetg $WORKINGDIR/Velvet_assembly/_"$n" \
-cov_cutoff 0 \
-exp_cov auto \
-min_contig_lgth 200 \
-ins_length 300 \
-clean yes \
-unused_reads no \
-min_pair_count 3 \
-scaffolding no

done

rm _*/CnyUnifiedSeq*
rm _*/*Graph*

cat _*/contigs.fa > all.contigs.fa

bowtie2 -p $nt -f -x $INDEXGENOMESDIR/$EBVgenome -U all.contigs.fa -S "$SAMPLE_NAME"_velvet_all_"$Type".sam;
samtools view -bS "$SAMPLE_NAME"_velvet_all_"$Type".sam > "$SAMPLE_NAME"_velvet_all_"$Type".bam;
samtools sort "$SAMPLE_NAME"_velvet_all_"$Type".bam "$SAMPLE_NAME"_velvet_all_"$Type"_sorted;
samtools index "$SAMPLE_NAME"_velvet_all_"$Type"_sorted.bam;


#/home/yk42w/project/share/bin_sync/trinityrnaseq_r20131110/util/fasta_to_tab.pl all.contigs.fa |sed 's/NODE/>NODE/g'|sed 's/\t/\n/g' > all.oneline.contigs.fa
else
echo "Skipping velvet without optimiser!!"
fi


#Filter out unmapped contigs
if [ "$PostvelvetOptimizer" = "1" ]
then

source /project/umw_jeffrey_bailey/share/bin_sync/PAGIT/PAGIT/sourceme.pagit

#Run Post velvet contig assembly processing to all kmer results

cd $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/

contigsFasta=contigs.fa

#Aligning assembled contigs to the reference EBV genomes
bowtie2 -p $nt -f -x $INDEXGENOMESDIR/$EBVgenome -U $contigsFasta -S "$SAMPLE_NAME"_velvet_"$Type".sam
samtools view -bS "$SAMPLE_NAME"_velvet_"$Type".sam > "$SAMPLE_NAME"_velvet_"$Type".bam
samtools sort "$SAMPLE_NAME"_velvet_"$Type".bam "$SAMPLE_NAME"_velvet_"$Type"_sorted
samtools index "$SAMPLE_NAME"_velvet_"$Type"_sorted.bam

#ABACAS
mkdir $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/runABACAS_First
cd $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/runABACAS_First
ln -s $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/$contigsFasta ./
ln -s $INDEXGENOMESDIR/$EBVgenome.fa ./

perl /project/umw_jeffrey_bailey/share/bin_sync/PAGIT/PAGIT/ABACAS/abacas.pl \
-r "$EBVgenome".fa \
-q $contigsFasta \
-p nucmer -b -t -c \
-o "$SAMPLE_NAME"_abacas

cat "$SAMPLE_NAME"_abacas.fasta "$SAMPLE_NAME"_abacas.contigsInbin.fas > "$SAMPLE_NAME"_abacas_mapAndUnmap.fasta

#IMAGE
mkdir $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/runIMAGE_First
cd $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/runIMAGE_First

ln -s $WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_1.fastq ./
ln -s $WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_2.fastq ./

ln -s $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/runABACAS_First/"$SAMPLE_NAME"_abacas_mapAndUnmap.fasta ./

perl /project/umw_jeffrey_bailey/share/bin_sync/PAGIT/PAGIT/IMAGE/image.pl \
-scaffolds "$SAMPLE_NAME"_abacas_mapAndUnmap.fasta \
-prefix "$SAMPLE_NAME"_Hg_filt_pairs \
-iteration 1 \
-all_iteration 10 \
-dir_prefix ite \
-kmer 97 \
-vel_ins_len 400

#Second ABACAS
mkdir $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/runABACAS_Second
cd $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/runABACAS_Second
ln -s $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/runIMAGE_First/ite10/new.fa ./
ln -s $INDEXGENOMESDIR/$EBVgenome.fa ./

contigsFasta=new.fa

perl /project/umw_jeffrey_bailey/share/bin_sync/PAGIT/PAGIT/ABACAS/abacas.pl \
-r "$EBVgenome".fa \
-q $contigsFasta \
-p nucmer -b -t -c \
-o "$SAMPLE_NAME"_abacas

mkdir $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/runIMAGE_last
cd $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/runIMAGE_last

ln -s $WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_1.fastq ./
ln -s $WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_2.fastq ./

ln -s $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/runABACAS_Second/"$SAMPLE_NAME"_abacas.fasta ./

scaffoldsFasta="$SAMPLE_NAME"_abacas.fasta
perl /project/umw_jeffrey_bailey/share/bin_sync/PAGIT/PAGIT/IMAGE/image.pl \
-scaffolds $scaffoldsFasta \
-prefix "$SAMPLE_NAME"_Hg_filt_pairs \
-iteration 1 -all_iteration 2 -dir_prefix ite -kmer 97 -vel_ins_len 400

cd $WORKINGDIR/

else
echo "Skipping Post velvet optimizer"
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

#*********************************************************************************************************#
if [ "$AlignStats" = "1" ]
then

cd $WORKINGDIR/ReadsMap2Ref/

#Collect Insert size 
java -Xmx32g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/CollectInsertSizeMetrics.jar \
INPUT=$WORKINGDIR/ReadsMap2Ref/"$SAMPLE_NAME"_alignment_$EBVgenome.sorted.bam \
OUTPUT=$WORKINGDIR/ReadsMap2Ref/"$SAMPLE_NAME"_$EBVgenome.AlignedReads_InsertSizeMetrics.txt \
HISTOGRAM_FILE=$WORKINGDIR/ReadsMap2Ref/"$SAMPLE_NAME"_$EBVgenome.AlignedReads_InsertSizeMetrics_histogram.pdf \
R=$INDEXGENOMESDIR/$EBVgenome.fa \
ASSUME_SORTED=true

#GC bias stats collection
java -Xmx32g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/CollectGcBiasMetrics.jar \
I=$WORKINGDIR/ReadsMap2Ref/"$SAMPLE_NAME"_alignment_$EBVgenome.sorted.bam \
R=$INDEXGENOMESDIR/$EBVgenome.fa \
O=GCmetrics.txt \
CHART=GC_Bias_Metrics_Plot_for_$EBVgenome.AlignedReads.pdf

cd $WORKINGDIR/

else
echo "Skipping align stats"
fi



#*********************************************************************************************************#
if [ "$GetConsensusFromReads" = "1" ]
then
export PERL5LIB=/project/umw_jeffrey_bailey/share/bin_sync/trinityrnaseq_r20131110/PerlLib/
module load python/2.7.5_packages/biopython/1.68

samtools index "$SAMPLE_NAME"_BART.bam

python /project/umw_jeffrey_bailey/OTHERS/dnaSeq_Tools/scripts/Reads2ConcensusGenome.py \
Reads2ConcensusGenome \
-n $EBVgenome \
--regions_to_mask $AnnotationDir/"$EBVgenome"_repeatMask.bed \
-r "$SAMPLE_NAME"_BART.bam \
-f $INDEXGENOMESDIR/"$EBVgenome".fa \
-o "$SAMPLE_NAME"_EBV_"$Type"

else
echo "no read cons."
fi



if [ "$JB_Merger" = "1" ]
then
export PERL5LIB=/project/umw_jeffrey_bailey/share/bin_sync/trinityrnaseq_r20131110/PerlLib/
module load python/2.7.5_packages/biopython/1.68
# Only version pysam 0.8.1 works... Installed in ~/.local/lib/python2.7/

mkdir $WORKINGDIR/Velvet_assembly/JB_ContigMerge_v3
cd $WORKINGDIR/Velvet_assembly/JB_ContigMerge_v3
#mkdir $WORKINGDIR/Velvet_assembly/JB_ContigMerge/test/
#cd $WORKINGDIR/Velvet_assembly/JB_ContigMerge/test/

#fastq_1=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_1.fastq
#fastq_2=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_2.fastq

ln -s $WORKINGDIR/Velvet_assembly/all.contigs.fa ./
/project/umw_jeffrey_bailey/share/bin_sync/trinityrnaseq_r20131110/util/fasta_remove_duplicates.pl all.contigs.fa |cut -f1 > all.contigs_noredun.fa
ln -s $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal.ba* ./

if [ "$BamDownSample" = "1" ]
then
java -Xmx32g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/DownsampleSam.jar \
INPUT="$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal.bam \
OUTPUT="$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal_down.bam \
PROBABILITY=0.05
#Sort sam file and output as bam
java -Xmx10g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/SortSam.jar \
INPUT="$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal_down.bam \
OUTPUT="$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal_down_sorted.bam \
SORT_ORDER=coordinate
#Index bam file
java -Xmx10g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/BuildBamIndex.jar \
INPUT="$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal_down_sorted.bam
#assign read bam file
ReadBamfile="$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal_down_sorted.bam
else
#assign read bam file
ReadBamfile="$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal.bam
echo "No downsampling"
fi

if [ -f $WORKINGDIR/Velvet_assembly/JB_ContigMerge/"$SAMPLE_NAME"_EBV_"$Type"_consensusSequence.fa ]
then
ln -s $WORKINGDIR/Velvet_assembly/JB_ContigMerge/"$SAMPLE_NAME"_EBV_"$Type"_consensusSequence.fa ./
else
python /project/umw_jeffrey_bailey/OTHERS/dnaSeq_Tools/scripts/Reads2ConcensusGenome.py \
Reads2ConcensusGenome \
-n $EBVgenome \
--regions_to_mask $AnnotationDir/"$EBVgenome"_repeatMask.bed \
-r $ReadBamfile \
-f $INDEXGENOMESDIR/"$EBVgenome".fa \
-o "$SAMPLE_NAME"_EBV_"$Type"
fi

python $toolDir/scripts/Assembly_Splitter_from_Ns.py "$SAMPLE_NAME"_EBV_"$Type"_consensusSequence.fa 1 "$SAMPLE_NAME"_EBV_"$Type"_consensusContigs.fa
cat "$SAMPLE_NAME"_EBV_"$Type"_consensusContigs.fa all.contigs_noredun.fa > "$SAMPLE_NAME"_EBV_"$Type"_NoRedundantContigs_and_ReadConsensus.fa

#INSERT CONSENSUS
#START HERE			!!!
#REMOVE REDUNDANT seq from all.contigs.fa

seqnum=`grep ">" "$SAMPLE_NAME"_EBV_"$Type"_NoRedundantContigs_and_ReadConsensus.fa |wc -l`;

for ((i=1;i<=$seqnum;i++));do let x=$i*2; head -$x "$SAMPLE_NAME"_EBV_"$Type"_NoRedundantContigs_and_ReadConsensus.fa|tail -2| ~/codes/EBVseq/fasta_chop.pl;done > "$SAMPLE_NAME"_EBV_"$Type"_NoRedundantContigs_and_ReadConsensus_chopped.fa

bowtie2 -p $nt -f -x $INDEXGENOMESDIR/$EBVgenome \
-U "$SAMPLE_NAME"_EBV_"$Type"_NoRedundantContigs_and_ReadConsensus_chopped.fa \
-S "$SAMPLE_NAME"_velvet_all_"$Type".sam

awk '($5 != 0)' "$SAMPLE_NAME"_velvet_all_"$Type".sam > "$SAMPLE_NAME"_velvet_all_"$Type"_noMQ0.sam

#Sort sam file and output as bam
java -Xmx10g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/SortSam.jar \
INPUT="$SAMPLE_NAME"_velvet_all_"$Type"_noMQ0.sam \
OUTPUT="$SAMPLE_NAME"_velvet_all_"$Type"_noMQ0_sorted.bam \
SORT_ORDER=coordinate
#Index bam file
java -Xmx10g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/BuildBamIndex.jar \
INPUT="$SAMPLE_NAME"_velvet_all_"$Type"_noMQ0_sorted.bam

#Merge overlapping contigs
python /project/umw_jeffrey_bailey/share/bin_sync/jbam_tools_v03.py \
ref_merge_contigs \
-n $EBVgenome \
-t 12 \
--regions_to_mask $AnnotationDir/"$EBVgenome"_repeatMask.bed \
-r $ReadBamfile \
-c "$SAMPLE_NAME"_velvet_all_"$Type"_noMQ0_sorted.bam \
-u 10 \
-f $INDEXGENOMESDIR/"$EBVgenome".fa \
-g 1 \
--min_unique_bp 60 \
-o "$SAMPLE_NAME"_EBV_"$Type"

#Contigs to Ref after merging
bowtie2 -p $nt -f -x $INDEXGENOMESDIR/$EBVgenome -U "$SAMPLE_NAME"_EBV_"$Type"_contig_new.fa -S "$SAMPLE_NAME"_contigs2Ref_JBmerged.sam
#Sort sam file and output as bam
java -Xmx10g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/SortSam.jar \
INPUT="$SAMPLE_NAME"_contigs2Ref_JBmerged.sam \
OUTPUT="$SAMPLE_NAME"_contigs2Ref_JBmerged_sorted.bam \
SORT_ORDER=coordinate
#Index bam file
java -Xmx10g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/BuildBamIndex.jar \
INPUT="$SAMPLE_NAME"_contigs2Ref_JBmerged_sorted.bam

else
echo "Skip JBmerger"
fi

if [ "$ContigsReadPileup" = "1" ]
then
rm -r $WORKINGDIR/Velvet_assembly/readsback2contigs
mkdir $WORKINGDIR/Velvet_assembly/readsback2contigs
cd $WORKINGDIR/Velvet_assembly/readsback2contigs

#ln -s $WORKINGDIR/Velvet_assembly/JB_ContigMerge_v2/"$SAMPLE_NAME"_EBV_"$Type"_contig_new.fa -S ./
ln -s $WORKINGDIR/Velvet_assembly/JB_ContigMerge_v3/"$SAMPLE_NAME"_EBV_"$Type"_contig_new.fa -S ./
#ln -s $WORKINGDIR/Velvet_assembly/JB_ContigMerge/test/"$SAMPLE_NAME"_EBV_"$Type"_contig_new.fa -S ./

ln -s $WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_1.fastq ./
ln -s $WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_2.fastq ./

fastq_1=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_1.fastq
fastq_2=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_2.fastq

bowtie2-build "$SAMPLE_NAME"_EBV_"$Type"_contig_new.fa "$SAMPLE_NAME"_EBV_"$Type"_contig_new_bowtieIndex
samtools faidx "$SAMPLE_NAME"_EBV_"$Type"_contig_new.fa

#Map reads back to contigs # This should be best match! so that reads will map to only one contig
bowtie2 -p $nt -x "$SAMPLE_NAME"_EBV_"$Type"_contig_new_bowtieIndex \
-1 $fastq_1 \
-2 $fastq_2 \
-S "$SAMPLE_NAME"_assembly_contigs_readsBack.sam

samtools view -bS "$SAMPLE_NAME"_assembly_contigs_readsBack.sam > "$SAMPLE_NAME"_assembly_contigs_readsBack.bam
samtools sort "$SAMPLE_NAME"_assembly_contigs_readsBack.bam "$SAMPLE_NAME"_assembly_contigs_readsBack_sorted
samtools index "$SAMPLE_NAME"_assembly_contigs_readsBack_sorted.bam

samtools mpileup -uf "$SAMPLE_NAME"_EBV_"$Type"_contig_new.fa "$SAMPLE_NAME"_assembly_contigs_readsBack_sorted.bam > rawcalls.bcf
/project/umw_jeffrey_bailey/share/bin_sync/samtools-0.1.19/bcftools/bcftools view -v rawcalls.bcf > contig.variants.vcf

cat contig.variants.vcf| \
grep -v "#"|awk '($6 > 100)'|awk 'NF { info=$8; gsub(/.*;DP4=|;MQ=.*/, "", info); split(info, a, /,/); print $0 "\t" a[1]"\t"a[2]"\t"a[3]"\t"a[4]}' |\
awk '{if($NF+$(NF-1) != 0) print $1"_"$2"_"$4"_"$5 "\t" ($(NF)+$(NF-1)+$(NF-2)+$(NF-3)) "\t" ($NF+$(NF-1)) / ($(NF)+$(NF-1)+$(NF-2)+$(NF-3))}'|\
awk '($2 > 99)'|sed 's/,/_/g' > minorfreq.txt

Rscript ~/codes/EBVseq/MAF_plot.R $SAMPLE_NAME minorfreq.txt


else
echo "no pileup"
fi



if [ "$AssemblyContigsReadPileup" = "1" ]
then
rm -r $WORKINGDIR/Velvet_assembly/readsback2contigs
mkdir $WORKINGDIR/Velvet_assembly/readsback2contigs
cd $WORKINGDIR/Velvet_assembly/readsback2contigs

InputAssembly=`grep -w "$SAMPLE_NAME" $WORKINGDIR/../MAFscreen.samples.txt| cut -f6`
InputBamForReads=`grep -w "$SAMPLE_NAME" $WORKINGDIR/../MAFscreen.samples.txt| cut -f7`

bowtie2-build $WORKINGDIR/../$InputAssembly $WORKINGDIR/../${InputAssembly%.fa}_bowtieIndex
samtools faidx $WORKINGDIR/../$InputAssembly

java -Xmx10g -XX:ParallelGCThreads=$nt -jar \
$PICARDPATH/SamToFastq.jar \
INPUT=$WORKINGDIR/../$InputBamForReads \
F=$WORKINGDIR/../${InputBamForReads%.bam}_1.fastq \
F2=$WORKINGDIR/../${InputBamForReads%.bam}_2.fastq \

bowtie2 -p $nt \
-x $WORKINGDIR/../${InputAssembly%.fa}_bowtieIndex \
-1 $WORKINGDIR/../${InputBamForReads%.bam}_1.fastq \
-2 $WORKINGDIR/../${InputBamForReads%.bam}_2.fastq \
-S "$SAMPLE_NAME"_assembly_contigs_readsBack.sam

samtools view -bS "$SAMPLE_NAME"_assembly_contigs_readsBack.sam > "$SAMPLE_NAME"_assembly_contigs_readsBack.bam
samtools sort "$SAMPLE_NAME"_assembly_contigs_readsBack.bam "$SAMPLE_NAME"_assembly_contigs_readsBack_sorted
samtools index "$SAMPLE_NAME"_assembly_contigs_readsBack_sorted.bam

samtools mpileup -uf $WORKINGDIR/../$InputAssembly "$SAMPLE_NAME"_assembly_contigs_readsBack_sorted.bam > rawcalls.bcf
/project/umw_jeffrey_bailey/share/bin_sync/samtools-0.1.19/bcftools/bcftools view -v rawcalls.bcf > contig.variants.vcf

cat contig.variants.vcf| \
grep -v "#"|awk '($6 > 100)'|awk 'NF { info=$8; gsub(/.*;DP4=|;MQ=.*/, "", info); split(info, a, /,/); print $0 "\t" a[1]"\t"a[2]"\t"a[3]"\t"a[4]}' |\
awk '{if($NF+$(NF-1) != 0) print $1"_"$2"_"$4"_"$5 "\t" ($(NF)+$(NF-1)+$(NF-2)+$(NF-3)) "\t" ($NF+$(NF-1)) / ($(NF)+$(NF-1)+$(NF-2)+$(NF-3))}'|\
awk '($2 > 99)'|sed 's/,/_/g' > minorfreq.txt

Rscript ~/codes/EBVseq/MAF_plot.R $SAMPLE_NAME minorfreq.txt
    
else
echo "no pileup"
fi




#*********************************************************************************************************#

#Run VFAT (Viral Finishing and Annotation Tool)
if [ "$Vfat_orient_merge" = "1" ]
then

export PATH=/project/umw_jeffrey_bailey/share/bin_sync/wise2.2.0/wisecfg/:$PATH
mkdir $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/runVfatMerger
cd $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/runVfatMerger

ln -s $WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_1.fastq ./
ln -s $WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_2.fastq ./
ln -s $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/runIMAGE_last/ite2/new.fa ./

VFATOutName='Vfat'
perl /project/umw_jeffrey_bailey/OTHERS/dnaSeq_Tools/resources/Vfat/orientContig.pl \
new.fa $INDEXGENOMESDIR/$EBVgenome.fa $VFATOutName

perl /project/umw_jeffrey_bailey/OTHERS/dnaSeq_Tools/resources/Vfat/contigMerger.pl \
"$VFATOutName"_orientedContigs \
$INDEXGENOMESDIR/$EBVgenome.fa \
"$VFATOutName"_merged

else
echo "not merging with Vfat"
fi
#*********************************************************************************************************#
#ICORN Iteration Number = it number -1
IIN=11

if [ "$ICORN" = "1" ]
then
source /project/umw_jeffrey_bailey/share/bin_sync/PAGIT/PAGIT/sourceme.pagit

mkdir $WORKINGDIR/Velvet_$velvetsettings/runICORN_VfatMerged/
cd $WORKINGDIR/Velvet_$velvetsettings/runICORN_VfatMerged/

ln -s $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/runVfatMerger/Vfat_merged_assembly.fa ./
ln -s $INDEXGENOMESDIR/$EBVgenome.fa ./

#ln -s $WORKINGDIR/Batch*.fastq.gz ./
#gunzip -d --force Batch*.fastq.gz

#ln -s $WORKINGDIR/ERR*.fastq.gz ./
#gunzip -d --force ERR*.fastq.gz

fastq_1=*1.fastq
fastq_2=*2.fastq

/project/umw_jeffrey_bailey/share/bin_sync/PAGIT/PAGIT/ICORN/icorn.start.sh \
Vfat_merged_assembly.fa 1 10 \
$fastq_1 $fastq_2 100,500 250

rm Batch*.fastq
rm ERR*.fastq

else
echo "Skipping Read mapping to assembled genome"
fi

#*********************************************************************************************************#
if [ "$Fix_FrameShifts" = "1" ]
then
Reads2Assembly=1

export PATH=/project/umw_jeffrey_bailey/share/bin_sync/wise2.2.0/wisecfg/:$PATH

mkdir $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/FixFrameshifts
cd $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/FixFrameshifts
rm "$SAMPLE_NAME"_abacas.fasta.$IIN

ln -s $WORKINGDIR/Velvet_$velvetsettings/runICORN_VfatMerged/Vfat_merged_assembly.fa.11 ./"$SAMPLE_NAME"_abacas.fasta.$IIN
assemblyFile="$SAMPLE_NAME"_abacas.fasta.$IIN

if [ "$Reads2Assembly" = "1" ]
then
fastq_1=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_1.fastq
fastq_2=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_2.fastq

bowtie2-build $assemblyFile ${assemblyFile%.fa}_bowtieIndex
bowtie2 -p $nt -x ${assemblyFile%.fa}_bowtieIndex \
-1 $fastq_1 \
-2 $fastq_2 \
-S ${assemblyFile%.fa}_readsBack.sam

samtools view -bS ${assemblyFile%.fa}_readsBack.sam > ${assemblyFile%.fa}_readsBack.bam
samtools sort ${assemblyFile%.fa}_readsBack.bam ${assemblyFile%.fa}_readsBack_sorted
samtools index ${assemblyFile%.fa}_readsBack_sorted.bam
else
echo "Skipping Read mapping to assembled genome"
fi

perl $VfatDIR/samToQlx.pl ${assemblyFile%.fa}_readsBack.sam $assemblyFile ${assemblyFile%.fa}_readsBack

#Fix Frame shifts
perl $VfatDIR/fixFrameshifts.pl \
-fa $assemblyFile \
-ref $INDEXGENOMESDIR/$EBVgenome.fa \
-genelist $RefDIR/EBV_Reference_genelist.txt \
-qlx ${assemblyFile%.fa}_readsBack.qlx \
-forcefix 0.20 \
-fixhomo \
-minhomosize 5 \
-o FixOut
#-fixhomo:	Forces correction of frameshifts in homopolymer regions. Will work without reads supplied, but requires at least 1 read support to fix if reads are supplied

#Annotate

rm ${assemblyFile%.fa}_readsBack.bam
rm ${assemblyFile%.fa}_readsBack.sam
rm ${assemblyFile%.fa}_readsBack.qlx

cd $WORKINGDIR

else
echo "Skipping VFAT"
fi

#*********************************************************************************************************#

if [ "$AssemblyContigsChecker" = "1" ]
then


mkdir $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/AssemblyCheck/
cd $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/AssemblyCheck/

ln -s $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/FixFrameshifts/FixOut_fixed_assembly.fa ./
ln -s $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/FixFrameshifts/FixOut_alignPair_fixed.afa ./

#Fix the edges of assembly
tip=`head -100 FixOut_alignPair_fixed.afa|\
sed ':a;N;$!ba;s/\n//g'|\
cut -d "A" -f1|\
tr -d -c '\-\n'|\
sed ':a;N;$!ba;s/\n//g'|\
awk '{print length}'`
let top=$tip+1;
grep -v ">" FixOut_fixed_assembly.fa|sed ':a;N;$!ba;s/\n//g'|cut -c$top- > top.tmp
grep -v ">" FixOut_fixed_assembly.fa|sed ':a;N;$!ba;s/\n//g'|cut -c-$tip > tip.tmp
echo ">""$SAMPLE_NAME" > FixOut_Edgefixed_assembly.fa

paste -d "" top.tmp tip.tmp >> FixOut_Edgefixed_assembly.fa
rm top.tmp tip.tmp

python $toolDir/scripts/Assembly_Splitter_from_Ns.py FixOut_Edgefixed_assembly.fa 1 "$SAMPLE_NAME"_assembly_contigs.fa;
#fastq_1=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_1.fastq
#fastq_2=$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_2.fastq

#bowtie2-build "$SAMPLE_NAME"_assembly_contigs.fa "$SAMPLE_NAME"_assembly_contigs_bowtieIndex

#Map reads back to contigs
#bowtie2 -p $nt -x "$SAMPLE_NAME"_assembly_contigs_bowtieIndex \
#-1 $fastq_1 \
#-2 $fastq_2 \
#-S "$SAMPLE_NAME"_assembly_contigs_readsBack.sam

#samtools view -bS "$SAMPLE_NAME"_assembly_contigs_readsBack.sam > "$SAMPLE_NAME"_assembly_contigs_readsBack.bam
#samtools sort "$SAMPLE_NAME"_assembly_contigs_readsBack.bam "$SAMPLE_NAME"_assembly_contigs_readsBack_sorted
#samtools index "$SAMPLE_NAME"_assembly_contigs_readsBack_sorted.bam

#Map contigs to Reference Genome
bowtie2 -p $nt -f -x $INDEXGENOMESDIR/$EBVgenome \
-U "$SAMPLE_NAME"_assembly_contigs.fa \
-S "$SAMPLE_NAME"_FinalAssembly_"$Type".sam

samtools view -bS "$SAMPLE_NAME"_FinalAssembly_"$Type".sam > "$SAMPLE_NAME"_FinalAssembly_"$Type".bam
samtools sort "$SAMPLE_NAME"_FinalAssembly_"$Type".bam "$SAMPLE_NAME"_FinalAssembly_"$Type"_sorted
samtools index "$SAMPLE_NAME"_FinalAssembly_"$Type"_sorted.bam

bedtools bamtobed -bed12 -i "$SAMPLE_NAME"_FinalAssembly_"$Type"_sorted.bam |\
cut -f1-3|awk '{print $0"\t""color=green"}' > "$SAMPLE_NAME"_FinalAssembly_contigs2Ref.bed12;

cd $WORKINGDIR

else
echo "Skipping Read mapping to assembled genome"
fi
#*********************************************************************************************************#

if [ "$ExtractGeneSequences" = "1" ]
then 

#mkdir $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/GeneSequences
#cd $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/GeneSequences
mkdir $WORKINGDIR/Velvet_assembly/GeneSequences
cd $WORKINGDIR/Velvet_assembly/GeneSequences

export PATH=/project/umw_jeffrey_bailey/share/bin_sync/wise2.2.0/wisecfg/:$PATH

InputAssembly=`grep -w "$SAMPLE_NAME" $WORKINGDIR/../MAFscreen.samples.txt| cut -f6`

#ln -s $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/AssemblyCheck/FixOut_Edgefixed_assembly.fa ./

perl $VfatDIR/annotate.pl \
-fa $WORKINGDIR/../$InputAssembly \
-ref $INDEXGENOMESDIR/$EBVgenome.fa \
-genelist $RefDIR/EBV_Reference_genelist.txt \
-pepfolder $RefDIR/EBV_Peptides \
-o vfat_annotation

else
echo "No gene extraction"
fi

#*********************************************************************************************************#

if [ "$ExtractInfo" = "1" ]
then

export PATH=/project/umw_jeffrey_bailey/share/bin_sync/wise2.2.0/wisecfg/:$PATH

cd $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/GeneSequences/vfat_annotation

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
#*********************************************************************************************************#

if [ "$SNPdat" = "1" ]
then
cd $WORKINGDIR/GATK/ploidy_10/

if [ "$type" = "1" ]
then
#FOR TYPE1
genome="NC_007605"
else
echo "Genome is Type II"
#FOR TYPE2; 
genome="NC_009334"
fi

perl /project/umw_jeffrey_bailey/share/bin_sync/SNPdat_package_v1.0.5/SNPdat_v1.0.5.pl \
-i "$SAMPLE_NAME"_alignment_"$genome".gatk_snp_filtered_PASSed.vcf \
-g /project/umw_jeffrey_bailey/share/EBV/Annotation/Type"$type"/custom_T"$type".gtf \
-f /project/umw_jeffrey_bailey/share/EBV/Genomes/Type"$type"/"$genome".fa \
-o "$SAMPLE_NAME"_SNPdat

else
echo "Skip SNPdat!"
fi

#*********************************************************************************************************#
if [ "$AA_Change_Rates" = "1" ]
then

#export PERL5LIB=/home/baileyj1/share/bin_sync/perl5lib/lib/perl5/
mkdir $WORKINGDIR/AAChange
cd $WORKINGDIR/AAChange
#export PATH=/home/baileyj1/share/bin_sync/samtools-0.1.19/:$PATH
samtools index $WORKINGDIR/ReadsMap2Ref/"$SAMPLE_NAME"_alignment_$EBVgenome.sorted.bam
perl /project/umw_jeffrey_bailey/share/bin_sync/josephhughes-btctools-03edf6b/btcutils.pl \
-bam $WORKINGDIR/ReadsMap2Ref/"$SAMPLE_NAME"_alignment_$EBVgenome.sorted.bam \
-ref $INDEXGENOMESDIR/$EBVgenome.fa \
-orfs $AnnotationDir/"$EBVgenome"_ORFs.txt \
-stub bct_"$SAMPLE_NAME"

else
echo "Skipping Vfat"
fi
#*********************************************************************************************************#
mergeBTC=0

if [ "$mergeBTC" = "1" ]
then
/project/umw_jeffrey_bailey/share/bin_sync/josephhughes-btctools-03edf6b/btcmerge_linux \
-files /project/umw_jeffrey_bailey/EBV/1583_04/GATK/bct_1583_04,/project/umw_jeffrey_bailey/EBV/BL534/GATK/bct_BL534,/project/umw_jeffrey_bailey/EBV/BL557/GATK/bct_BL557,/project/umw_jeffrey_bailey/EBV/BL573/GATK/bct_BL573,/project/umw_jeffrey_bailey/EBV/BL574/GATK/bct_BL574,/project/umw_jeffrey_bailey/EBV/BL577/GATK/bct_BL577,/project/umw_jeffrey_bailey/EBV/BL579/GATK/bct_BL579,/project/umw_jeffrey_bailey/EBV/BL606/GATK/bct_BL606,/project/umw_jeffrey_bailey/EBV/BL607/GATK/bct_BL607,/project/umw_jeffrey_bailey/EBV/BL609/GATK/bct_BL609,/project/umw_jeffrey_bailey/EBV/BL614/GATK/bct_BL614,/project/umw_jeffrey_bailey/EBV/BL620/GATK/bct_BL620,/project/umw_jeffrey_bailey/EBV/BL627/GATK/bct_BL627,/project/umw_jeffrey_bailey/EBV/BL628/GATK/bct_BL628,/project/umw_jeffrey_bailey/EBV/eBL099_plasma/GATK/bct_eBL099_plasma,/project/umw_jeffrey_bailey/EBV/Namalwa/GATK/bct_Namalwa \
-out BTCout_T1
/project/umw_jeffrey_bailey/share/bin_sync/josephhughes-btctools-03edf6b/btcmerge_linux \
-files /project/umw_jeffrey_bailey/EBV/1381_22/GATK/bct_1381_22,/project/umw_jeffrey_bailey/EBV/BL552/GATK/bct_BL552,/project/umw_jeffrey_bailey/EBV/BL564/GATK/bct_BL564,/project/umw_jeffrey_bailey/EBV/BL576/GATK/bct_BL576,/project/umw_jeffrey_bailey/EBV/BL578/GATK/bct_BL578,/project/umw_jeffrey_bailey/EBV/BL587/GATK/bct_BL587,/project/umw_jeffrey_bailey/EBV/BL629/GATK/bct_BL629,/project/umw_jeffrey_bailey/EBV/eBL103_plasma/GATK/bct_eBL103_plasma,/project/umw_jeffrey_bailey/EBV/eBL108_plasma/GATK/bct_eBL108_plasma \
-out BTCout_T2
else
echo "Skip btc merge!"
fi
#*********************************************************************************************************#
if [ "$CircosPlotFiles" = "1" ]
then
#variation files
gunzip -c $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_"$EBVgenome".gatk_snp_filtered_PASSed.vcf.gz|grep -v "#" |awk '{print "'$EBVgenome'""\t"$2"\t"$2+1"\t""1.0"}' > $WORKINGDIR/GATK/"$SAMPLE_NAME".circosnp.txt
gunzip -c $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_"$EBVgenome".gatk_indel_filtered_PASSed.vcf.gz|grep -v "#" |awk '{print "'$EBVgenome'"\t"$2"\t"$2+length($5)"\t""1.0"}'> $WORKINGDIR/GATK/"$SAMPLE_NAME".circosindel.txt
##Rub this code below for all files.
#cat *circosnp.txt| awk '{print $1"\t"$2"\t"$3"\t""snp""\t"$4"\t""+"}'|sort -k1,1 -k2,2g > tmp.bed
#genomeCoverageBed -d -i tmp.bed -g ebvgenome |awk ' {if(NR%100 >= 1)sum +=$3; else sum=0; print $0"\t"sum"\t"$1"\t"$2-99"\t"$2}'|awk '(NR%99 == 0){print $5"\t"$6"\t"$7"\t"$4}'|less
#file ebvgenome has this
#AJ507799        171823


#coverage files
igvtools count -w 5 $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal.bam $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal.tdf $INDEXGENOMESDIR/"$EBVgenome".fa

genomeCoverageBed -d -ibam $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal.bam -g $INDEXGENOMESDIR/"$EBVgenome".fa |\
awk '{OFS="\t"; print $1,$2,$2+1,log($3)/log(10)}' |\
awk '{if(($4*1) <= 1.0) x="color=vvlred";
	else if(1.0 < ($4*1) && ($4*1) <= 2.0 ) x="color=vlred";
	else if(2.0 < ($4*1) && ($4*1) <= 3.0 ) x="color=lred";
	else if(3.0 < ($4*1) && ($4*1) <= 4.0 ) x="color=red";
	else if(4.0 < ($4*1) && ($4*1) <= 5.0 ) x="color=dred";
	else if(5.0 < ($4*1) && ($4*1) <= 6.0 ) x="color=vdred";
	else x="no"; print $1"\t"$2"\t"$3"\t"$4"\t"x}' > $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal.bgcolor


else
echo "Not creating circos files"
fi

#*********************************************************************************************************#

if [ "$SimpleStats" = "1" ]
then

cd $WORKINGDIR/GATK
java -Xmx16g -XX:ParallelGCThreads=$nt -jar \
$GATKdir/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
-o $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal_DepthofCoverage \
-I $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal.bam \
-geneList $RefDIR/EBV_Reference_genelist.txt

#STATS
else
echo "Skipping Stats"
fi

#*********************************************************************************************************#
#Run VFAT (Viral Finishing and Annotation Tool)
if [ "$Vfat" = "1" ]
then

export PATH=/project/umw_jeffrey_bailey/share/bin_sync/wise2.2.0/wisecfg/:$PATH
cd $WORKINGDIR/Velvet_$velvetsettings/

ln -s $WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_1.fastq ./
ln -s $WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_2.fastq ./

perl /project/umw_jeffrey_bailey/share/bin_sync/Vfat/fqpair2fasta.pl \
$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_1.fastq \
$WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_2.fastq \
READS


mkdir $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/
mkdir $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/Vfat
cd $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/Vfat

ln -s $WORKINGDIR/Velvet_$velvetsettings/READS.* ./
ln -s $WORKINGDIR/Velvet_$velvetsettings/"$option"_outDir/runIMAGE_last/ite2/new.fa ./

contigFile=new.fa

perl $VfatDIR/vfat.pl \
-contigs $contigFile \
-readfa READS_1.fa \
-readq READS_1.qual \
-readfa2 READS_2.fa \
-readq2 READS_2.qual \
-ref $INDEXGENOMESDIR/$EBVgenome.fa \
-genelist $RefDIR/EBV_Reference_genelist.txt \
-pepfolder $RefDIR/EBV_Peptides \
-mincontlen 110 \
-virus EBV \
-details \
-o vfatOut

cd $WORKINGDIR/Velvet_$velvetsettings/


else
echo "Skipping Vfat"
fi

#*********************************************************************************************************#


if [ "$CleanIntermediateFiles" = "1" ]
then


#Check if the output is generated correctly then, remove intermediates:
s=`ls -1 $WORKINGDIR/GATK/|grep "$EBVgenome".sorted_RGadded_fixed.bam |wc -l`
if [ $s = "1" ]
then
rm $WORKINGDIR/ReadsMap2Ref/*alignment_*.sam
rm $WORKINGDIR/GATK/*alignment_$EBVgenome.sorted_DD.ba*
rm $WORKINGDIR/GATK/*alignment_$EBVgenome.sorted_RGadded.ba*
else
echo "Proper output files not found"
fi
####CleanUp####


rm -rf $WORKINGDIR/*_AdpTrimmed_*.fastq
rm -rf $WORKINGDIR/*goodComplex_*.fastq
rm -rf $WORKINGDIR/*_lowComplex_*.fastq
rm -fr $WORKINGDIR/*_lowQual_*.fastq

else
echo "NO need to clean"
fi
#*********************************************************************************************************#

#gzip $WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_pairs_*.fastq

