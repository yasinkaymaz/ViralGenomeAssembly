#!/bin/bash
#BSUB -n 6
#BSUB -R rusage[mem=8000]
#BSUB -R "span[hosts=1]"
#BSUB -q long
#BSUB -W 36:00
#BSUB -o "/project/umw_jeffrey_bailey/yk42w/std_out/%J.out"
#BSUB -e "/project/umw_jeffrey_bailey/yk42w/std_err/%J.err"

#####BSUB -m blades

#export PATH=/home/yk42w/biotools/annovar/:$PATH
#export PATH=/home/yk42w/pipeline/tabix-0.2.6/:$PATH
#export PATH=/home/yk42w/pipeline/vcftools_0.1.9/bin/:$PATH
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/picard-tools-1.105/:$PATH
#export PATH=/home/yk42w/pipeline/R-3.0.1/bin/:$PATH
PICARDPATH="/project/umw_jeffrey_bailey/share/bin_sync/picard-tools-1.105/"
#load modules
module load java/1.7.0_25
module load samtools/0.0.19
module load tophat/2.0.12
module load bowtie2/2-2.1.0
module load R/2.14.2

HOME=/home/yk42w
### working directory

dir=`pwd`
b=`basename $dir`
SAMPLE_NAME=$b
$LSB_JOBNAME=$SAMPLE_NAME

echo "job_id: $LSB_JOBID"
echo "index_id: $LSB_JOBINDEX"
echo "job_name: $LSB_JOBNAME"

###########################################
#Main Control Panel

TR=4

group="eBLs"
#group="sBLs"
#group="LCLs"

Stranded=1
Diagnose=1
DenovoAss=0
Ebv1Tophat=0
Ebv2Tophat=0
EBV1Rsem=0
EBV2Rsem=0
Map2DenovoTx=0
ReadMapping_tophat=1
ReadMapping_star=1
Hg19Tophat=0
FusionSearch=0
RSeQC=0
DataPreProcess=0
Reorder=0
RNASeQC=0
RepeatFilter=0
RsemRefCreate=0
Rsem=1

#### Sequencing data: ###
fastq_1=*1.fastq
fastq_2=*2.fastq
##########################################
RESULTSDIR="/project/umw_jeffrey_bailey/OTHERS/DENIZ/All_results/RNAseq"
#RESULTSDIR="/project/umw_jeffrey_bailey/yk42w/results/RNAseq/TotalRNA"
BOWTIE2ebvDIR="/project/umw_jeffrey_bailey/share/EBV/Bowtie2Index"
BOWTIE2hg19DIR="/project/umw_jeffrey_bailey/share/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index"
BOWTIEhg19DIR="/project/umw_jeffrey_bailey/share/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex"
GTFebvDIR="/project/umw_jeffrey_bailey/share/EBV/Annotation"
GTFhg19DIR="/project/umw_jeffrey_bailey/share/Homo_sapiens/UCSC/hg19/Annotation/Genes"
BOWTIEhg19RepeatDIR="/project/umw_jeffrey_bailey/share/Homo_sapiens/Repeats/Bowtie"
BOWTIE2hg19RepeatDIR="/project/umw_jeffrey_bailey/share/Homo_sapiens/Repeats/Bowtie2"
RSEMhg19FASTADIR="/project/umw_jeffrey_bailey/share/Homo_sapiens/UCSC/hg19/RSEM"
RSEMEBVFASTADIR="/project/umw_jeffrey_bailey/share/EBV/RSEM"
RSEMsubDirNM="/project/umw_jeffrey_bailey/yk42w/results/RNAseq/RSEM/RR/NM/PolyA"
EBVRSEMsubDir="/project/umw_jeffrey_bailey/yk42w/results/RNAseq/PolyA/EBV_Transcriptome/RSEM"
#RSEMsubDirNM="/project/umw_jeffrey_bailey/yk42w/results/RNAseq/RSEM/RR/NM"
#RSEMsubDir="/project/umw_jeffrey_bailey/yk42w/results/RNAseq/RSEM/RR/NMNR/"
RSEMsubDir="/project/umw_jeffrey_bailey/yk42w/results/RNAseq/RSEM/polyA/NMNR"
Gencode="/project/umw_jeffrey_bailey/share/Homo_sapiens/Gencode"

if [ "$Stranded" = "1" ]
then
library=fr-firststrand
insert=300
strand=FIRST_READ_TRANSCRIPTION_STRAND
Read_1_target='--nofw'
forwardProb=0
strandRule="-d '1+-,1-+,2++,2--'"
else
library=fr-unstranded
#library=fr-secondstrand
insert=300
#strand=NONE
#strand=SECOND_READ_TRANSCRIPTION_STRAND
Read_1_target=''
forwardProb=0.5
#forwardProb=1
strandRule=""
fi

if [ "$Diagnose" = "1" ]
then
module load fastqc/0.10.1
module load cutadapt/1.3

fastqc $fastq_1 --outdir RESULTSDIR/
fastqc $fastq_2 --outdir RESULTSDIR/

cutadapt -b ATCTCGTATGCCGTCTTCTGCTTG \
-e 0.1 -O 5 -o "$SAMPLE_NAME"_cutAdpted_1.fastq $fastq_1

cutadapt -b ATCTCGTATGCCGTCTTCTGCTTG \
-e 0.1 -O 5 -o "$SAMPLE_NAME"_cutAdpted_2.fastq $fastq_2

#fastqc "$SAMPLE_NAME"_cutAdpted_1.fastq
#fastqc "$SAMPLE_NAME"_cutAdpted_2.fastq

#gzip eBLRNAseqBatch3_"$b"_R1.fastq
#gzip eBLRNAseqBatch3_"$b"_R2.fastq

else
echo "Skipping Diagnostic step!!!"
fi

### Run Trinity for DeNovo Assembly
#_____________________________________________________________________________________________________________________________________#
if [ "$DenovoAss" == "1" ]
then
#cd $@
#cd /project/umw_jeffrey_bailey/yk42w/results/RNAseq/PolyA/EBV_Transcriptome/trinity_bams
cd /project/umw_jeffrey_bailey/OTHERS/DENIZ/human_testis_deniz/data/RSQ/

fastq_1=*_1.fastq
fastq_2=*_2.fastq

dir=`pwd`
SAMPLE_NAME=`basename $dir`

echo $SAMPLE_NAME
echo $fastq_1
echo $fastq_2

### Run Trinity ###

for r1 in *.1.fastq; do
    prefix=${r1%.1.fastq}
    r2=${prefix}.2.fastq
    echo "Trinity --max_memory 100G --seqType fq --left ${r1} --right ${r2} --CPU 16 --SS_lib_type RF --full_cleanup --output ${prefix}.trinity &> ${prefix}.trinity.log"
done

Trinity.pl --seqType fq --left $fastq_1 --right $fastq_2 \
--full_cleanup \
--JM 100G \
--CPU 8 \
--inchworm_cpu 8 \
--bflyHeapSpaceMax 40G \
--bflyCPU 8 \
--output /project/umw_jeffrey_bailey/OTHERS/DENIZ/trinity/

module load blat/35x1
blat /project/umw_jeffrey_bailey/share/EBV/Bowtie2Index/Type1/NC_007605.fa EBV_Type1_trinity_transcript_assemly.Trinity.fasta EBV_Type1_trinity_transcript_assemly.Trinity.psl
awk '{print $14"\t"$16"\t"$17"\t"$10"\t"$11"\t"$9"\t"$16"\t"$17"\t""255,0,0""\t"$18"\t"$19"\t"$20}' EBV_Type1_trinity_transcript_assemly.Trinity.psl > EBV_Type1_trinity_transcript_assemly.Trinity.bed



else
echo "Skipping de novo assembly!!!"
fi





if [ "$ReadMapping_tophat" = "1" ]
then
mkdir $RESULTSDIR/TopHat/tophatOut_Hg19_"$SAMPLE_NAME"

tophat2 -p $TR --output-dir $RESULTSDIR/TopHat/tophatOut_Hg19_"$SAMPLE_NAME" \
--library-type=$library -r $insert \
$BOWTIE2hg19DIR/$i \
$fastq_1 \
$fastq_2 \
-G $Gencode/v19/gencode.v19.annotation.gtf

cd $RESULTSDIR/TopHat/tophatOut_Hg19_"$SAMPLE_NAME"/
mv accepted_hits.bam "$SAMPLE_NAME"_alignment_"$i".bam
samtools sort "$SAMPLE_NAME"_alignment_"$i".bam "$SAMPLE_NAME"_alignment_"$i"_sorted
samtools index "$SAMPLE_NAME"_alignment_"$i"_sorted.bam

else
echo "Skipping read Mapping with TopHat!"
fi

if [ "$ReadMapping_star" = "1" ]
then

if [ -f *1.fastq.gz ];
then
   echo "Need to unzip first";
gunzip -d --force *.fastq.gz
else
   echo "No need to unzip.";
fi

RESULTSDIR=$RESULTSDIR/STAR/

fastq_1=*_1.fastq
fastq_2=*_2.fastq
mkdir $RESULTSDIR/Star_"$SAMPLE_NAME"
mkdir $RESULTSDIR/Star_"$SAMPLE_NAME"/1pass

#First Pass
STAR --genomeDir $STARhg19DIR --readFilesIn $fastq_1 $fastq_2 --runThreadN $TR \
--outFileNamePrefix $RESULTSDIR/Star_"$SAMPLE_NAME"/1pass/

#Second Pass
mkdir $RESULTSDIR/Star_"$SAMPLE_NAME"/2pass/
mkdir $RESULTSDIR/Star_"$SAMPLE_NAME"/2pass/hg19_2pass
TwoPassGenomeDir="$RESULTSDIR/Star_"$SAMPLE_NAME"/2pass/hg19_2pass"
#For the 2-pass STAR, a new index is then created using splice junction information contained in the file SJ.out.tab from the first pass
STAR --runMode genomeGenerate --genomeDir $TwoPassGenomeDir \
--genomeFastaFiles $BOWTIE2hg19DIR/$i.fa \
--sjdbFileChrStartEnd $RESULTSDIR/Star_"$SAMPLE_NAME"/1pass/SJ.out.tab \
--sjdbOverhang 75 \
--runThreadN $TR

STAR --genomeDir $TwoPassGenomeDir --readFilesIn $fastq_1 $fastq_2 --runThreadN $TR \
--outFileNamePrefix $RESULTSDIR/Star_"$SAMPLE_NAME"/2pass/
cd $RESULTSDIR/Star_"$SAMPLE_NAME"/2pass/
mv Aligned.out.sam "$SAMPLE_NAME"_alignment_"$i".bam

samtools sort "$SAMPLE_NAME"_alignment_"$i".bam "$SAMPLE_NAME"_alignment_"$i"_sorted
samtools index "$SAMPLE_NAME"_alignment_"$i"_sorted.bam

else
echo "Skipping read Mapping with STAR!"
fi







### Run TopHat-2 with Bowtie2 for EBV Type-1
#_____________________________________________________________________________________________________________________________________#
if [ "$Ebv1Tophat" = "1" ]
then
#fastq_1=*cutAdpted_1.fastq
#fastq_2=*cutAdpted_2.fastq

tophat2 -p $TR --output-dir tophatOut_EBVType-1_$b \
--library-type=$library -r $insert \
$BOWTIE2ebvDIR/wt_ebv \
$fastq_1 \
$fastq_2 \
-G $GTFebvDIR/EBV_WT.gb2.gff

cd tophatOut_EBVType-1_$b
rm unmapped*
mv accepted_hits.bam accepted_hits_EBVType-1_$b.bam
#Remove Duplicates
#First need to sort file:
java -jar $PICARDPATH/SortSam.jar \
INPUT=accepted_hits_EBVType-1_"$b".bam \
OUTPUT=accepted_hits_EBVType-1_"$b"_sorted.bam \
SORT_ORDER=coordinate
#To remove duplicates:
java -jar $PICARDPATH/MarkDuplicates.jar \
INPUT=accepted_hits_EBVType-1_"$b"_sorted.bam \
OUTPUT=accepted_hits_EBVType-1_"$b"_sorted_DupRmved.bam \
METRICS_FILE=RmvDups_outfile_Metrics_"$b" \
REMOVE_DUPLICATES=true

samtools index accepted_hits_EBVType-1_"$b"_sorted.bam
samtools index accepted_hits_EBVType-1_"$b"_sorted_DupRmved.bam

cd ../
mv tophatOut_EBVType-1_$b $RESULTSDIR/EBV_Transcriptome/
#cd RESULTSDIR/EBV_Transcriptome/tophatOut_EBVType-1_$b
#Run cufflinks
#cufflinks accepted_hits.bam ...

else
echo "Skipping read alignment to EBV Type-1"
fi
#_________________()__________________#

### Run TopHat-2 with Bowtie2 for EBV Type-2
#_____________________________________________________________________________________________________________________________________#
if [ "$Ebv2Tophat" = "1" ]
then
fastq_1=*cutAdpted_1.fastq
fastq_2=*cutAdpted_2.fastq

tophat2 -p $TR --output-dir tophatOut_EBVType-2_$b \
--library-type=$library -r $insert \
$BOWTIE2ebvDIR/AG876_ebv \
$fastq_1 \
$fastq_2 \
-G $GTFebvDIR/EBV_WT.gb2.gff

cd tophatOut_EBVType-2_$b
rm unmapped*

mv accepted_hits.bam accepted_hits_EBVType-2_$b.bam
#Remove Duplicates
#First need to sort file:
java -jar $PICARDPATH/SortSam.jar \
INPUT=accepted_hits_EBVType-2_"$b".bam \
OUTPUT=accepted_hits_EBVType-2_"$b"_sorted.bam \
SORT_ORDER=coordinate
#To remove duplicates:
java -jar $PICARDPATH/MarkDuplicates.jar \
INPUT=accepted_hits_EBVType-2_"$b"_sorted.bam \
OUTPUT=accepted_hits_EBVType-2_"$b"_sorted_DupRmved.bam \
METRICS_FILE=RmvDups_outfile_Metrics_"$b" \
REMOVE_DUPLICATES=true

samtools index accepted_hits_EBVType-2_"$b"_sorted.bam
samtools index accepted_hits_EBVType-2_"$b"_sorted_DupRmved.bam

cd ../
mv tophatOut_EBVType-2_$b $RESULTSDIR/EBV_Transcriptome/
#cd RESULTSDIR/EBV_Transcriptome/tophatOut_EBVType-1_$b
#Run cufflinks
#cufflinks accepted_hits.bam ...

else
echo "Skipping read alignment to EBV Type-2"
fi
#_________________()__________________#

### Map RNAseq reads to known viral transcripts using RSEM.
#______________________________________________________________________________________#
#Type1
if [ "$EBV1Rsem" = "1" ]
then
module unload R/2.14.2
#Load modules
module load bowtie/1.0.0
module load RSEM/1.2.11
module load R/3.0.1
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/jdk1.6.0_45/bin/:$PATH
fastq_1=*cutAdpted_1.fastq
fastq_2=*cutAdpted_2.fastq

mkdir $EBVRSEMsubDir/RsemType1_$SAMPLE_NAME/

cd $EBVRSEMsubDir/RsemType1_$SAMPLE_NAME/

rsem-calculate-expression -p $TR \
--paired-end \
$dir/$fastq_1 \
$dir/$fastq_2 \
$RSEMEBVFASTADIR/Type1/EBV_type1_genome_reference.rsem \
$SAMPLE_NAME.rsem \
--no-bam-output \
--forward-prob $forwardProb

#unload modules
module unload bowtie/1.0.0
module unload RSEM/1.2.11
module unload R/3.0.1

else
echo "Skipping Rsem estimations of viral expressions for type1!!!"
fi
#Type2
if [ "$EBV2Rsem" = "1" ]
then
module unload R/2.14.2
#Load modules
module load bowtie/1.0.0
module load RSEM/1.2.11
module load R/3.0.1
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/jdk1.6.0_45/bin/:$PATH
fastq_1=*cutAdpted_1.fastq
fastq_2=*cutAdpted_2.fastq

mkdir $EBVRSEMsubDir/RsemType2_$SAMPLE_NAME/

cd $EBVRSEMsubDir/RsemType2_$SAMPLE_NAME/

rsem-calculate-expression -p $TR \
--paired-end \
$dir/$fastq_1 \
$dir/$fastq_2 \
$RSEMEBVFASTADIR/Type2/EBV_type2_genome_reference.rsem \
$SAMPLE_NAME.rsem \
--no-bam-output \
--forward-prob $forwardProb

#unload modules
module unload bowtie/1.0.0
module unload RSEM/1.2.11
module unload R/3.0.1

else
echo "Skipping Rsem estimations of viral expressions for type2!!!"
fi
#_________________()__________________#

### Map RNAseq reads to DeNovo assembled transcript sequences
#____________________________________________________________________________________________________________________________________#
if [ "$Map2DenovoTx" = "1" ]
then
module unload bowtie2/2-2.1.0
module unload R/2.14.2
module unload R/3.1.0
module load R/3.0.1
module load bowtie/1.0.0
module load RSEM/1.2.11
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/trinityrnaseq_r20131110/util/RSEM_util/:$PATH
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/trinityrnaseq_r20131110/util/:$PATH
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/trinityrnaseq_r20131110/:$PATH

fastq_1=*cutAdpted_1.fastq
fastq_2=*cutAdpted_2.fastq
TrinityOut='/project/umw_jeffrey_bailey/yk42w/results/RNAseq/PolyA/EBV_Transcriptome/'
mkdir $RESULTSDIR/EBV_Transcriptome/Trinity_Rsem_$SAMPLE_NAME

run_RSEM_align_n_estimate.pl --transcripts $RESULTSDIR/EBV_Transcriptome/Expressed_transcripts.Trinity.fasta \
--seqType fq \
--left $fastq_1 \
--right $fastq_2 \
--SS_lib_type RF \
--prefix RSEMcountTrinity -- \
--no-bam-output 
mv RSEMcountTrinity* $RESULTSDIR/EBV_Transcriptome/Trinity_Rsem_$SAMPLE_NAME/

module unload bowtie/1.0.0
module unload RSEM/1.2.11
module load bowtie2/2-2.1.0
module load R/2.14.2
else
echo "Skipping DeNovo transcriptome expression quantification!!!"
fi
#_________________()__________________#

### Run TopHat-2 with Bowtie2 for Human
#_____________________________________________________________________________________________________________________________________#
if [ "$Hg19Tophat" = "1" ]
then
fastq_1=*cutAdpted_1.fastq
fastq_2=*cutAdpted_2.fastq

tophat2 -p $TR --output-dir tophatOut_Hg19_$b \
--library-type=$library -r $insert \
$BOWTIE2hg19DIR/genome \
$fastq_1 \
$fastq_2 \
-G $GTFhg19DIR/genes.gtf

cd tophatOut_Hg19_$b

mv accepted_hits.bam accepted_hits_Hg19_$b.bam
#Remove Duplicates
#First need to sort file:
java -jar $PICARDPATH/SortSam.jar \
INPUT=accepted_hits_Hg19_"$b".bam \
OUTPUT=accepted_hits_Hg19_"$b"_sorted.bam \
SORT_ORDER=coordinate
#To remove duplicates:
java -jar $PICARDPATH/MarkDuplicates.jar \
INPUT=accepted_hits_Hg19_"$b"_sorted.bam \
OUTPUT=accepted_hits_Hg19_"$b"_sorted_DupRmved.bam \
METRICS_FILE=RmvDups_outfile_Metrics_"$b" \
REMOVE_DUPLICATES=true

#Alignment Metrics
java -jar $PICARDPATH/CollectAlignmentSummaryMetrics.jar \
I=accepted_hits_Hg19_"$b"_sorted_DupRmved.bam \
O=alignmentMetrics_$b.out \
R=/project/umw_jeffrey_bailey/share/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome.fa

java -jar $PICARDPATH/CollectRnaSeqMetrics.jar \
REF_FLAT=/project/umw_jeffrey_bailey/share/Homo_sapiens/UCSC/hg19/Annotation/Genes/refFlat.txt.gz \
STRAND_SPECIFICITY=FIRST_READ_TRANSCRIPTION_STRAND \
CHART_OUTPUT="$b"_chart \
INPUT=accepted_hits_Hg19_"$b"_sorted_DupRmved.bam \
OUTPUT="$b"_TranscriptUniformity_Metrix_out \
MINIMUM_LENGTH=500 RRNA_FRAGMENT_PERCENTAGE=0.8 METRIC_ACCUMULATION_LEVEL=ALL_READS \
ASSUME_SORTED=true \
STOP_AFTER=0 \
VERBOSITY=INFO \
QUIET=false \
VALIDATION_STRINGENCY=STRICT \
COMPRESSION_LEVEL=5 \
MAX_RECORDS_IN_RAM=500000 \
CREATE_INDEX=false \
CREATE_MD5_FILE=false \

#Plot the RNAseq metrics output
Rscript /home/yk42w/codes/TranscriptUniformityMetrics.R "$b"_TranscriptUniformity_Metrix_out "$b"_TranscriptUniformity.pdf Sample_$b

cd ../
mv tophatOut_Hg19_$b $RESULTSDIR/Hg19_Transcriptome/


else
echo "Skipping read alignment to Hg19!!!"
fi
#Run cufflinks
#cufflinks -p 8 -o cufflinksOut_$b accepted_hits_"$b"_sorted.bam

#_________________()__________________#

if [ "$FusionSearch" = "1" ]
then
fastq_1=*cutAdpted_1.fastq
fastq_2=*cutAdpted_2.fastq

module load bowtie/1.0.0

tophat2 -p $TR --output-dir tophat_$b \
--fusion-search \
--keep-fasta-order \
--bowtie1 \
--no-coverage-search \
-r 0 --mate-std-dev 80 \
--max-intron-length 100000 \
--fusion-min-dist 100000 \
--fusion-anchor-length 13 \
--fusion-ignore-chromosomes chrM \
$BOWTIEhg19DIR/genome \
$fastq_1 \
$fastq_2 \
-G $GTFhg19DIR/genes.gtf

mv tophat_$b $RESULTSDIR/Hg19_Transcriptome/Fusion/$group

module unload bowtie/1.0.0
else
echo "Skipping Fusion transcript search!!!"
fi

module unload java/1.7.0_25
module unload samtools/0.0.19
module unload tophat/2.0.12
module unload bowtie2/2-2.1.0
module unload R/2.14.2



if [ "$RSeQC" = "1" ]
then
#RSeQC Metrics
module load python/2.7.5
module load samtools/0.0.19
module load R/2.14.2

export PYTHONPATH=/home/yk42w/.local/lib/python2.7/site-packages:$PYTHONPATH
export PATH=/home/yk42w/.local/bin:$PATH
cd $RESULTSDIR/tophatOut_Hg19_$b/

#samtools sort "$b"_alignment_genome.bam "$b"_alignment_genome_sorted
#samtools index "$b"_alignment_genome_sorted.bam
#infer_experiment.py -r $GTFhg19DIR/hg19_RefSeq.bed -i "$b"_alignment_genome_sorted.bam
#read_distribution.py -i "$b"_alignment_genome_sorted.bam -r $GTFhg19DIR/hg19_RefSeq.bed > "$b"_read_dist.txt
#split_bam.py -i "$b"_alignment_genome_sorted.bam -r /project/umw_jeffrey_bailey/share/Homo_sapiens/UCSC/hg19/hg19_rRNA.bed -o outputrRNA
#read_duplication.py -i "$b"_alignment_genome_sorted.bam -o outputDup
#clipping_profile.py -i "$b"_alignment_genome_sorted.bam -o outputClip
#junction_annotation.py -i "$b"_alignment_genome_sorted.bam -o outputJunction -r $GTFhg19DIR/hg19_RefSeq.bed
#junction_saturation.py -i "$b"_alignment_genome_sorted.bam -o outputJuncSaturation -r $GTFhg19DIR/hg19_RefSeq.bed
junction_saturation.py -i accepted_hits_Hg19_"$b"_sorted_DupRmved.bam -o outputJuncSaturation -r $GTFhg19DIR/hg19_RefSeq.bed
#read_GC.py -i "$b"_alignment_genome_sorted.bam -o outputGC
#read_NVC.py -i "$b"_alignment_genome_sorted.bam -o outputNVC
#read_quality.py -i "$b"_alignment_genome_sorted.bam -o outputreadQual
#RPKM_count.py -i "$b"_alignment_genome_sorted.bam -d $strandRule -r $GTFhg19DIR/hg19_RefSeq.bed -o outputFPKM
#RPKM_saturation.py -i "$b"_alignment_genome_sorted.bam -r $GTFhg19DIR/hg19_RefSeq.bed -o "$b"_outputFPKMsaturation
RPKM_saturation.py -i accepted_hits_Hg19_"$b"_sorted_DupRmved.bam -d '1+-,1-+,2++,2--' -r $GTFhg19DIR/hg19_RefSeq.bed -o outputFPKMsaturation

else
echo "Skipping RSeQC!!!"
fi

############ Pre-process alignment file #############
#___________________________________________________#
if [ "$DataPreProcess" = "1" ]
then
module load java/1.7.0_25

cd $RESULTSDIR/hgtophats/

#Sort bam file and output as sorted.bam
java -Xmx20g -jar \
$PICARDPATH/SortSam.jar \
INPUT="$SAMPLE_NAME"_alignment_genome.bam \
OUTPUT="$SAMPLE_NAME"_alignment_genome_sorted.bam \
VALIDATION_STRINGENCY=LENIENT \
SORT_ORDER=coordinate

#Run Picard Tools
java -Xmx20g -jar \
$PICARDPATH/AddOrReplaceReadGroups.jar \
I="$SAMPLE_NAME"_alignment_genome_sorted.bam \
O="$SAMPLE_NAME"_alignment_genome_RGadded.bam \
SORT_ORDER=coordinate \
CREATE_INDEX=true \
RGPL=ILLUMINA \
RGID=$LSB_JOBID \
RGSM=$SAMPLE_NAME \
RGLB=RNAseq \
RGPU=Run_Barcode \
VALIDATION_STRINGENCY=LENIENT

#MARK DUPLICATES
java -Xmx20g -jar \
$PICARDPATH/MarkDuplicates.jar \
INPUT="$SAMPLE_NAME"_alignment_genome_RGadded.bam \
OUTPUT="$SAMPLE_NAME"_alignment_genome_DupMarked.bam \
VALIDATION_STRINGENCY=LENIENT \
METRICS_FILE="$SAMPLE_NAME"_genome_dedup_metrics.txt

#Fix Mate pair
java -Xmx31g -jar \
$PICARDPATH/FixMateInformation.jar \
INPUT="$SAMPLE_NAME"_alignment_genome_DupMarked.bam \
OUTPUT="$SAMPLE_NAME"_alignment_genome_MateFixed.bam \
SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=LENIENT

#Sort bam file and output as sorted.bam
java -Xmx20g -jar \
$PICARDPATH/SortSam.jar \
INPUT="$SAMPLE_NAME"_alignment_genome_MateFixed.bam \
OUTPUT="$SAMPLE_NAME"_alignment_genome_sorted_final.bam \
VALIDATION_STRINGENCY=LENIENT \
SORT_ORDER=coordinate

else
echo "skip preprocessing!!!"
fi

if [ "$Reorder" = "1" ]
then
module load java/1.7.0_25

cd $RESULTSDIR/hgtophats/

#Reorder bam file in karyotypic order rather than Lexicographical order.
#java -Xmx10g -jar \
#$PICARDPATH/ReorderSam.jar \
#REFERENCE=$BOWTIE2hg19DIR/genome.fa \
#INPUT="$SAMPLE_NAME"_alignment_genome_sorted_final.bam \
#VALIDATION_STRINGENCY=LENIENT \
#OUTPUT="$SAMPLE_NAME"_alignment_alignment_genome_sorted_final_KaryoSorted.bam

#Index bam file
java -Xmx10g -jar \
$PICARDPATH/BuildBamIndex.jar \
VALIDATION_STRINGENCY=LENIENT \
INPUT="$SAMPLE_NAME"_alignment_alignment_genome_sorted_final_KaryoSorted.bam

else
echo "No Karyo sort!!"
fi



if [ "$RNASeQC" = "1" ]
then
module load java/1.7.0_25
module load RNA-SeQC/1.1.7
module load bwa/0.7.5a
JarPATH="/share/pkg/RNA-SeQC/1.1.7/"
cd $RESULTSDIR/hgtophats/

#RNAmetrics
java -Xmx20g -jar $JarPATH/RNA-SeQC_v1.1.7.jar \
-BWArRNA $Gencode/human_all_rRNA.fasta \
-gc $Gencode/gencode.v7.gc.txt \
-o RNA_SeQC \
-r $BOWTIE2hg19DIR/genome.fa \
-s /home/yk42w/codes/rnaseqc-2.txt \
-t $Gencode/gencode.v7.annotation.gtf \
-transcriptDetails



module unload java/1.7.0_25
module unload RNA-SeQC/1.1.7
module unload bwa/0.7.5a
else
echo "Skipping RNASeQC!!!"
fi



#Remove reads map to Repeat regions and rRNA.
#------------------------------------------------------
if [ "$RepeatFilter" = "1" ]
then
#module load bowtie/1.0.0
module load bowtie2/2-2.1.0
module load samtools/0.0.19
module load bedtools/2.17.0

fastq_1=*cutAdpted_1.fastq
fastq_2=*cutAdpted_2.fastq

mkdir $RSEMsubDir/Rsem_$SAMPLE_NAME
mkdir $RSEMsubDir/Rsem_$SAMPLE_NAME/repeatOut

bowtie2 -p $TR  --fr $Read_1_target -x $BOWTIE2hg19RepeatDIR/Hg19_Repeat_Masker_100ntLarger \
--un-conc $RSEMsubDir/Rsem_$SAMPLE_NAME/repeatOut/filt_pairs \
--un $RSEMsubDir/Rsem_$SAMPLE_NAME/repeatOut/filt_single \
-1 $fastq_1 -2 $fastq_2 -S $RSEMsubDir/Rsem_$SAMPLE_NAME/"$SAMPLE_NAME"_repeats.sam

#bowtie -S --fr $Read_1_target $BOWTIEhg19RepeatDIR/Hg19_Repeat_Masker_100ntLarger -1 $fastq_1 -2 $fastq_2 > \
#$RSEMsubDir/Rsem_$SAMPLE_NAME/"$SAMPLE_NAME"_repeats_alignment.sam

#cd $RSEMsubDir/Rsem_$SAMPLE_NAME/

#samtools view -hSf 4 "$SAMPLE_NAME"_repeats_alignment.sam > "$SAMPLE_NAME"_repeatFiltered_alignment.sam
#samtools view -bS "$b"_repeatFiltered_alignment.sam > "$b"_repeatFiltered_alignment.bam
#samtools sort -n "$b"_repeatFiltered_alignment.bam "$b"_repeatFiltered_alignment_sorted

#bedtools bamtofastq -i "$b"_repeatFiltered_alignment_sorted.bam -fq "$b"_repeatFiltered_alignment_sorted_1.fq -fq2 "$b"_repeatFiltered_alignment_sorted_2.fq

#convert bam to paired end fastq
#java -Xmx20g -jar \
#$PICARDPATH/SamToFastq.jar \
#I="$b"_repeatFiltered_alignment_sorted.bam \
#F="$b"_repeatFiltered_alignment_sorted_1.fq \
#F2="$b"_repeatFiltered_alignment_sorted_2.fq

module unload bowtie/1.0.0
module unload samtools/0.0.19
module unload bedtools/2.17.0
else
echo "Skipping Repeat filtering!!!"
fi
#-------------------------------------------------------

if [ "$RsemRefCreate" = "1" ]
then
module load bowtie2/2-2.1.0
module load RSEM/1.2.11
module load R/3.0.1

rsem-prepare-reference --gtf /home/yk42w/project/share/Homo_sapiens/Gencode/v19/gencode.v19.annotation.gtf \
--bowtie2 \
/project/umw_jeffrey_bailey/share/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome.fa \
gencode.v19.rsem

#For Long-noncoding Transcripts
rsem-prepare-reference --gtf /project/umw_jeffrey_bailey/share/Homo_sapiens/Gencode/v19/gencode.v19.long_noncoding_RNAs.gtf \
--bowtie2 \
/project/umw_jeffrey_bailey/share/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome.fa \
gencode.v19.NC.rsem


else
echo "No Rsem ref create!"
fi



#RSEM
#_____________________________________________________________________________________________________________________________________#
if [ "$Rsem" = "1" ]
then
#Load modules
module load bowtie2/2-2.1.0
module load RSEM/1.2.11
module load R/3.0.1
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/jdk1.6.0_45/bin/:$PATH

mkdir $RSEMsubDir/Rsem_$SAMPLE_NAME/

cd $RSEMsubDir/Rsem_$SAMPLE_NAME/

#Run Rsem first pass to get transcriptome alignments to be able to determine PCR dup reads locations.
rsem-calculate-expression -p $TR \
--paired-end \
$RSEMsubDir/Rsem_$SAMPLE_NAME/repeatOut/*.1 \
$RSEMsubDir/Rsem_$SAMPLE_NAME/repeatOut/*.2 \
$RSEMhg19FASTADIR/hg19genome_reference.rsem \
$SAMPLE_NAME.rsem \
--forward-prob $forwardProb

#To remove duplicates:
java -jar $PICARDPATH/MarkDuplicates.jar \
INPUT="$b".rsem.transcript.sorted.bam \
OUTPUT="$b".rsem.transcript.sorted_DupRmved.bam \
METRICS_FILE=RmvDups_outfile_Metrics_"$b" \
REMOVE_DUPLICATES=true

mkdir $RSEMsubDir/Rsem_$SAMPLE_NAME/Dedupped
mv "$b".rsem.transcript.sorted_DupRmved.bam $RSEMsubDir/Rsem_$SAMPLE_NAME/Dedupped

cd $RSEMsubDir/Rsem_$SAMPLE_NAME/Dedupped
#Run RSEM for gene and Isoform level quantification
rsem-calculate-expression -p $TR \
--paired-end \
--bam \
"$b".rsem.transcript.sorted_DupRmved.bam \
$RSEMhg19FASTADIR/hg19genome_reference.rsem \
$SAMPLE_NAME.rsem \
--no-bam-output \
--forward-prob $forwardProb

#CleanUp
#Check if the output is generated correctly then, remove intermediates:
s=`ls -1 $RSEMsubDir/Rsem_$SAMPLE_NAME/Dedupped/|grep isoforms.results|wc -l`
if [ $s = "1" ]
then
rm $RSEMsubDir/Rsem_$SAMPLE_NAME/Dedupped/*.bam
rm $RSEMsubDir/Rsem_$SAMPLE_NAME/*.bam
rm $RSEMsubDir/Rsem_$SAMPLE_NAME/*.sam
rm $RSEMsubDir/Rsem_$SAMPLE_NAME/*.fq
else
echo "RSEM output files not found!"
fi


else
echo "Skipping RSEM expression Quantification!!!"
fi
#Unload modules
module unload bowtie/1.0.0
module unload RSEM/1.2.11
module unload R/3.0.1
#_________________()__________________#
#######################################
