#!/bin/bash
#BSUB -n 10
#BSUB -R rusage[mem=15000]
#BSUB -R "span[hosts=1]"
#BSUB -q long
#BSUB -W 10:00
#BSUB -o "/project/umw_jeffrey_bailey/OTHERS/std_out/%J.out"
#BSUB -e "/project/umw_jeffrey_bailey/OTHERS/std_err/%J.err"
#BSUB -R select[tmp>1000]

# Set the number of threads for tools that can be run in parallel.
#This is usually one digit below the number of cores this script is reserved to use (-n value above).
nt=10

#Set the working directory and Sample name
dir=`pwd`
SAMPLE_NAME=`basename $dir`
LSB_JOBNAME=$SAMPLE_NAME
echo "job_id: $LSB_JOBID"
echo "index_id: $LSB_JOBINDEX"
echo "job_name: $LSB_JOBNAME"


################################################
############# Main Control Panel: ##############
################################################

#Setting the directories for required files
WORKINGDIR=$dir

#Choose the type of the EBV genome!
type=1
#type=2

#fastq_1=$WORKINGDIR/*1.fastq.gz
#fastq_2=$WORKINGDIR/*2.fastq.gz

#fastq_1=$WORKINGDIR/*R1_001.fastq.gz
#fastq_2=$WORKINGDIR/*R2_001.fastq.gz

# READ Pre-PROCESS: This function will automatically
# 1) Trims adapter sequeces from both reads,
# 2) Excludes low complexity reads (poly As, Ts etc.),
# 3) Trims low quality based off of reads starting from 3'ends,
# 4) After each step of precessing, checks the quality/length distribution of reads.
### Please Turn this on if you need by changing '0' to '1' 
ReadPreProc=0
FilterHumanReads=0
# READ Mapping to Reference genomic Sequence:
#ReadMapping2Ref=0

#Denovo assembly with velvet
#Compiling velvet with specific settings
#compile=0
#velvet=0
#PostVelvet=0

velvetOptimizer=1
PostvelvetOptimizer=0

#Post de-novo assembly process
JB_Merger=0
Vfat=0
FixAnnotate=0

SimpleStats=0
AssemblyContigsChecker=0
ICORN=0
CircosPlotFiles=0

# Remove intermediate processing files to save space. Keep this off for de-bugging purposes.
### Please Turn this on if you need by changing '0' to '1'
CleanIntermediateFiles=0

#array=(Lbp Lcon max tbp n50)
array=(n50)
#LNbp = The total number of Ns in large contigs
#Lbp = The total number of base pairs in large contigs
#Lcon = The number of large contigs
#max = The length of the longest contig
#n50 = The n50
#ncon = The total number of contigs
#tbp = The total number of basepairs in contigs


##########################################
######### NO NEED TO MODIFY ! ############
############# Sub-routines ###############
##########################################


#Load required modules
module load ncbi_cxx/12_0_0
module load cutadapt/1.7.1
module load samtools/0.0.19
module load bowtie2/2-2.1.0
module load fastx_toolkit/0.0.14
module load java/1.7.0_25
module load fastqc/0.10.1
module load prinseq/0.20.4
module load gcc/4.8.1
module load python/2.7.5
module load bzip2/1.0.6
module load tabix/0.2.6
module load IGVTools/2.3.31
module load R/3.0.1
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
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/bamtools-master/bin/:$PATH
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


if [ "$ReadPreProc" = "1" ]
then

#Do Fastqc check
#fastqc $fastq_1
#fastqc $fastq_2

#Remove adaptors
#cutadapt -a ATCTCGTATGCCGTCTTCTGCTTG \
#-a GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCCGATGTTT \
#-e 0.1 -O 5 \
#-o "$SAMPLE_NAME"_AdpTrimmed_1.fastq \
#-p "$SAMPLE_NAME"_AdpTrimmed_2.fastq \
#$fastq_1 \
#$fastq_2

fastq_1="$SAMPLE_NAME"_AdpTrimmed_1.fastq

#Do Fastqc check
#fastqc $fastq_1
#fastqc $fastq_2

#Trim low qual sequence off of reads
perl /share/pkg/prinseq/0.20.4/prinseq-lite.pl \
-fastq $fastq_1 \
-out_good "$SAMPLE_NAME"_goodQual \
-out_bad "$SAMPLE_NAME"_lowQual \
-log \
-trim_qual_right 20

fastq_1=$WORKINGDIR/"$SAMPLE_NAME"_goodQual.fastq

#Do Fastqc check
#fastqc $fastq_1


cd $WORKINGDIR/
else
echo "Skipping Read Preprocessing step"
fi


if [ "$FilterHumanReads" = "1" ]
then

cd $WORKINGDIR/
fastq_1=$WORKINGDIR/"$SAMPLE_NAME"_goodQual.fastq


#Align reads to Hg19 and keep unmapped reads.
BOWTIE2hg19DIR="/project/umw_jeffrey_bailey/share/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index"

bowtie2 -p $nt -x $BOWTIE2hg19DIR/genome \
--un "$SAMPLE_NAME"_Hg_filt_single \
-U $fastq_1 -S "$SAMPLE_NAME"_Hg_mapped.bam

#samtools view -bS "$SAMPLE_NAME"_Hg_mapped.sam > "$SAMPLE_NAME"_Hg_mapped.bam
#samtools sort "$SAMPLE_NAME"_Hg_mapped.bam "$SAMPLE_NAME"_Hg_mapped_sorted
#samtools index "$SAMPLE_NAME"_Hg_mapped_sorted.bam

#fastqc "$SAMPLE_NAME"_Hg_filt_pairs*

mv "$SAMPLE_NAME"_Hg_filt_single "$SAMPLE_NAME"_Hg_filt_single.fastq

else
echo "no human filtration..."
fi

#########################################################
######  Alternative De Novo assembly with VelVet  #######
#########################################################


if [ "$velvetOptimizer" = "1" ]
then
mkdir $WORKINGDIR/Velvet_v2
cd $WORKINGDIR/Velvet_v2

ln -s $WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_single.fastq ./

VelvetHString='-create_binary -short -fastq '$SAMPLE_NAME'_Hg_filt_single.fastq';
VelvetGString='-read_trkg yes -min_contig_lgth 100'
#VelvetGString='-cov_cutoff 2 -read_trkg yes -min_contig_lgth 100'
#VelvetGString='-cov_cutoff 5 -read_trkg yes -min_contig_lgth 100'

for option in "${array[@]}";
do
echo "optimizing assembly for $option";

VelvetOptimiser.pl -s 21 -e 75 \
-f "$VelvetHString" \
-o "$VelvetGString" \
-t $nt \
-k $option \
-c "Lbp" \
-d ""$option"_v2";

rm CnyUnifiedSeq*
rm *Graph*

cd "$option"_v2/;

perl $toolDir/scripts/contigs_stats.pl -t Velvet \
contigs.fa -plot > contigs_stats.txt;

bowtie2 -p $nt -f -x $INDEXGENOMESDIR/$EBVgenome -U contigs.fa -S "$SAMPLE_NAME"_velvet_"$option".sam;
samtools view -bS "$SAMPLE_NAME"_velvet_"$option".sam > "$SAMPLE_NAME"_velvet_"$option".bam;
samtools sort "$SAMPLE_NAME"_velvet_"$option".bam "$SAMPLE_NAME"_velvet_"$option"_sorted;
samtools index "$SAMPLE_NAME"_velvet_"$option"_sorted.bam;

bedtools bamtobed -bed12 -i "$SAMPLE_NAME"_velvet_"$option"_sorted.bam > "$SAMPLE_NAME"_velvet_"$option"_sorted.bed
cd ../;

done

else
echo "Skipping velvetoptimiser";
fi


#Filter out unmapped contigs
if [ "$PostvelvetOptimizer" = "1" ]
then

source /project/umw_jeffrey_bailey/share/bin_sync/PAGIT/PAGIT/sourceme.pagit

#Run Post velvet contig assembly processing to all kmer results

for option in "${array[@]}";
do

cd $WORKINGDIR/Velvet/"$option"_outDir/

contigsFasta=contigs.fa

#Aligning assembled contigs to the reference EBV genomes
bowtie2 -p $nt -f -x $INDEXGENOMESDIR/$EBVgenome -U $contigsFasta -S "$SAMPLE_NAME"_velvet_"$Type".sam
samtools view -bS "$SAMPLE_NAME"_velvet_"$Type".sam > "$SAMPLE_NAME"_velvet_"$Type".bam
samtools sort "$SAMPLE_NAME"_velvet_"$Type".bam "$SAMPLE_NAME"_velvet_"$Type"_sorted
samtools index "$SAMPLE_NAME"_velvet_"$Type"_sorted.bam

#ABACAS
mkdir $WORKINGDIR/Velvet/"$option"_outDir/runABACAS_First
cd $WORKINGDIR/Velvet/"$option"_outDir/runABACAS_First
ln -s $WORKINGDIR/Velvet/"$option"_outDir/$contigsFasta ./
ln -s $INDEXGENOMESDIR/$EBVgenome.fa ./

perl /project/umw_jeffrey_bailey/share/bin_sync/PAGIT/PAGIT/ABACAS/abacas.pl \
-r "$EBVgenome".fa \
-q $contigsFasta \
-p nucmer -b -t -c \
-o "$SAMPLE_NAME"_abacas

cat "$SAMPLE_NAME"_abacas.fasta "$SAMPLE_NAME"_abacas.contigsInbin.fas > "$SAMPLE_NAME"_abacas_mapAndUnmap.fasta

#IMAGE
mkdir $WORKINGDIR/Velvet/"$option"_outDir/runIMAGE_First
cd $WORKINGDIR/Velvet/"$option"_outDir/runIMAGE_First

ln -s $WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_single.fastq ./

ln -s $WORKINGDIR/Velvet/"$option"_outDir/runABACAS_First/"$SAMPLE_NAME"_abacas_mapAndUnmap.fasta ./

perl /project/umw_jeffrey_bailey/share/bin_sync/PAGIT/PAGIT/IMAGE/image.pl \
-scaffolds "$SAMPLE_NAME"_abacas_mapAndUnmap.fasta \
-prefix $WORKINGDIR/"$SAMPLE_NAME"_Hg_filt_single \
-iteration 1 \
-all_iteration 10 \
-dir_prefix ite \
-kmer 30

#Second ABACAS
mkdir $WORKINGDIR/Velvet/"$option"_outDir/runABACAS_Second
cd $WORKINGDIR/Velvet/"$option"_outDir/runABACAS_Second
ln -s $WORKINGDIR/Velvet/"$option"_outDir/runIMAGE_First/ite10/new.fa ./
ln -s $INDEXGENOMESDIR/$EBVgenome.fa ./

contigsFasta=new.fa

perl /project/umw_jeffrey_bailey/share/bin_sync/PAGIT/PAGIT/ABACAS/abacas.pl \
-r "$EBVgenome".fa \
-q $contigsFasta \
-p nucmer -b -t -c \
-o "$SAMPLE_NAME"_abacas

mkdir $WORKINGDIR/Velvet/"$option"_outDir/runIMAGE_last
cd $WORKINGDIR/Velvet/"$option"_outDir/runIMAGE_last

ln -s $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recall_1.fastq ./
ln -s $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recall_2.fastq ./

ln -s $WORKINGDIR/Velvet/"$option"_outDir/runABACAS_Second/"$SAMPLE_NAME"_abacas.fasta ./

scaffoldsFasta="$SAMPLE_NAME"_abacas.fasta
perl /project/umw_jeffrey_bailey/share/bin_sync/PAGIT/PAGIT/IMAGE/image.pl \
-scaffolds $scaffoldsFasta \
-prefix "$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recal \
-iteration 1 -all_iteration 2 -dir_prefix ite -kmer 30

cd $WORKINGDIR/

done

else
echo "Skipping Post velvet optimizer"
fi


#Run VFAT (Viral Finishing and Annotation Tool)
if [ "$Vfat" = "1" ]
then

export PATH=/project/umw_jeffrey_bailey/share/bin_sync/wise2.2.0/wisecfg/:$PATH
cd $WORKINGDIR/Velvet/

ln -s $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recall_1.fastq ./
ln -s $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recall_2.fastq ./

perl /project/umw_jeffrey_bailey/share/bin_sync/Vfat/fqpair2fasta.pl \
"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recall_1.fastq \
"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recall_2.fastq \
READS


for option in "${array[@]}";
do

mkdir $WORKINGDIR/Velvet/"$option"_outDir/
mkdir $WORKINGDIR/Velvet/"$option"_outDir/Vfat
cd $WORKINGDIR/Velvet/"$option"_outDir/Vfat

ln -s $WORKINGDIR/Velvet/READS.* ./
ln -s $WORKINGDIR/Velvet/"$option"_outDir/runIMAGE_last/ite2/new.fa ./

contigFile=new.fa

perl $VfatDIR/vfat.pl \
-contigs $contigFile \
-readfa READS_1.fa \
-readq READS_1.qual \
-readfa2 READS_2.fa \
-readq2 READS_2.qual \
-ref $RefDIR/EBV_Reference.fasta \
-genelist $RefDIR/EBV_Reference_genelist.txt \
-pepfolder $RefDIR/EBV_Peptides \
-mincontlen 110 \
-virus EBV \
-details \
-o vfatOut

cd $WORKINGDIR/Velvet/

done

else
echo "Skipping Vfat"
fi


if [ "$FixAnnotate" = "1" ]
then
Reads2Assembly=1

export PATH=/project/umw_jeffrey_bailey/share/bin_sync/wise2.2.0/wisecfg/:$PATH


for option in "${array[@]}";
do

mkdir $WORKINGDIR/Velvet/"$option"_outDir/Vfat/FixAnnotate
cd $WORKINGDIR/Velvet/"$option"_outDir/Vfat/FixAnnotate

assemblyFile=$WORKINGDIR/Velvet/"$option"_outDir/Vfat/vfatOut_fixed_assembly.fa

if [ "$Reads2Assembly" = "1" ]
then

fastq_1=$WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recall_1.fastq
fastq_2=$WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recall_2.fastq

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
-ref $RefDIR/EBV_Reference.fasta \
-genelist $RefDIR/EBV_Reference_genelist.txt \
-qlx ${assemblyFile%.fa}_readsBack.qlx \
-forcefix 0.20 \
-o FixOut
#-fixhomo	Forces correction of frameshifts in homopolymer regions. Will work without reads supplied, but requires at least 1 read support to fix if reads are supplied

#Annotate
cd -
done


else
echo "Skipping VFAT"
fi



if [ "$AssemblyContigsChecker" = "1" ]
then


for option in "${array[@]}";
do

mkdir $WORKINGDIR/Velvet/"$option"_outDir/Vfat/FixAnnotate/AssemblyCheck/
cd $WORKINGDIR/Velvet/"$option"_outDir/Vfat/FixAnnotate/AssemblyCheck/

ln -s $WORKINGDIR/Velvet/"$option"_outDir/Vfat/FixAnnotate/FixOut_fixed_assembly.fa ./
ln -s $WORKINGDIR/Velvet/"$option"_outDir/Vfat/FixAnnotate/FixOut_alignPair_fixed.afa ./

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

fastq_1=$WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recall_1.fastq
fastq_2=$WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recall_2.fastq

bowtie2-build "$SAMPLE_NAME"_assembly_contigs.fa "$SAMPLE_NAME"_assembly_contigs_bowtieIndex

#Map reads back to contigs
bowtie2 -p $nt -x "$SAMPLE_NAME"_assembly_contigs_bowtieIndex \
-1 $fastq_1 \
-2 $fastq_2 \
-S "$SAMPLE_NAME"_assembly_contigs_readsBack.sam

samtools view -bS "$SAMPLE_NAME"_assembly_contigs_readsBack.sam > "$SAMPLE_NAME"_assembly_contigs_readsBack.bam
samtools sort "$SAMPLE_NAME"_assembly_contigs_readsBack.bam "$SAMPLE_NAME"_assembly_contigs_readsBack_sorted
samtools index "$SAMPLE_NAME"_assembly_contigs_readsBack_sorted.bam

#Map contigs to Reference Genome
bowtie2 -p $nt -f -x $INDEXGENOMESDIR/$EBVgenome \
-U "$SAMPLE_NAME"_assembly_contigs.fa \
-S "$SAMPLE_NAME"_FinalAssembly_"$Type".sam

samtools view -bS "$SAMPLE_NAME"_FinalAssembly_"$Type".sam > "$SAMPLE_NAME"_FinalAssembly_"$Type".bam
samtools sort "$SAMPLE_NAME"_FinalAssembly_"$Type".bam "$SAMPLE_NAME"_FinalAssembly_"$Type"_sorted
samtools index "$SAMPLE_NAME"_FinalAssembly_"$Type"_sorted.bam

bedtools bamtobed -bed12 -i "$SAMPLE_NAME"_FinalAssembly_"$Type"_sorted.bam |\
cut -f1-3|awk '{print $0"\t""color=green"}' > "$SAMPLE_NAME"_FinalAssembly_contigs2Ref.bed12;

#Run VFAT for final
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/wise2.2.0/wisecfg/:$PATH

contigFile="$SAMPLE_NAME"_assembly_contigs.fa
ln -s $WORKINGDIR/Velvet/READS.* ./

perl $VfatDIR/vfat.pl \
-contigs $contigFile \
-readfa READS_1.fa \
-readq READS_1.qual \
-readfa2 READS_2.fa \
-readq2 READS_2.qual \
-ref $RefDIR/EBV_Reference.fasta \
-genelist $RefDIR/EBV_Reference_genelist.txt \
-pepfolder $RefDIR/EBV_Peptides \
-mincontlen 110 \
-virus EBV \
-details \
-o vfatOut

cd -
done

else
echo "Skipping Read mapping to assembled genome"
fi


if [ "$ICORN" = "1" ]
then
source /project/umw_jeffrey_bailey/share/bin_sync/PAGIT/PAGIT/sourceme.pagit

mkdir $WORKINGDIR/Vicuna/runICORN/
cd $WORKINGDIR/Vicuna/runICORN/

ln -s $WORKINGDIR/Vicuna/AssemblyCheck/"$SAMPLE_NAME"_assembly_contigs.fa ./

ln -s $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recall_1.fastq ./
ln -s $WORKINGDIR/GATK/"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recall_2.fastq ./

/project/umw_jeffrey_bailey/share/bin_sync/PAGIT/PAGIT/ICORN/icorn.start.sh \
"$SAMPLE_NAME"_assembly_contigs.fa 1 5 \
"$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recall_1.fastq "$SAMPLE_NAME"_alignment_$EBVgenome.gatk_recall_2.fastq 100,500 250

cd -

else
echo "Skipping Read mapping to assembled genome"
fi



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


if [ "$SNVFasta2Genes" = "1" ]
then

export PATH=/project/umw_jeffrey_bailey/share/bin_sync/wise2.2.0/wisecfg/:$PATH

mkdir $WORKINGDIR/SNVFasta2Genes
cd $WORKINGDIR/SNVFasta2Genes

#Vfat annotate
perl $VfatDIR/annotate.pl \
-fa $WORKINGDIR/SNV2fasta/"$SAMPLE_NAME"_AlternativeGenome_with_gatk_filterPassed_SNVs.fasta \
-ref $RefDIR/EBV_Reference.fasta \
-genelist $RefDIR/EBV_Reference_genelist.txt \
-pepfolder $RefDIR/EBV_Peptides/ \
-align $WORKINGDIR/SNVFasta2Genes/GeneSequences_aligned.afa \
-o GeneSequences

cd GeneSequences

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
echo "Skipping SNVFasta to Genes"
fi


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




