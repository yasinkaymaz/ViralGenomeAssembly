#!/bin/bash

#BSUB -n 3
#BSUB -R rusage[mem=24000]
#BSUB -R "span[hosts=1]"
#BSUB -q long
#BSUB -W 36:00
#BSUB -o "/home/co54w/std_out/%J.out"
#BSUB -e "/home/co54w/std_err/%J.err"
#BSUB -R select[tmp>1000]

module load vicuna/1.3
module load ncbi_cxx/12_0_0
module load cutadapt/1.3
module load samtools/0.0.19
module load bowtie2/2-2.1.0
module load fastx_toolkit/0.0.14
module load java/1.7.0_25
module load fastqc/0.10.1
module load prinseq/0.20.4
module load perl/5.18.1
module load velvet/1.2.10
module load bedtools/2.22.0

ebvT1DIR="/project/umw_jeffrey_bailey/OTHERS/EBVseq/Ref_Genome/Type1"
BOWTIE2hg19DIR="/project/umw_jeffrey_bailey/OTHERS/EBVseq/hg19/Bowtie2Index"
BOWTIE2ebvT1="/project/umw_jeffrey_bailey/OTHERS/EBVseq/Bowtie2Index/Type1"

TR=2
dir=`pwd`
SAMPLE_NAME=`basename $dir`

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/pkg/ncbi_cxx/12_0_0/lib
export PERL5LIB=/project/umw_jeffrey_bailey/share/bin_sync/perl5lib/

gunzip *fastq.gz

fastq_1=*1.fastq
fastq_2=*2.fastq

#Remove Illumina adaptors sequences
cutadapt -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCC \
-b ATCTCGTATGCCGTCTTCTGCTTG \
-e 0.1 -O 5 -o "$SAMPLE_NAME"_AdpTrimmed_1.fastq $fastq_1

cutadapt -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCC \
-b ATCTCGTATGCCGTCTTCTGCTTG \
-e 0.1 -O 5 -o "$SAMPLE_NAME"_AdpTrimmed_2.fastq $fastq_2

fastq_1=*_AdpTrimmed_1.fastq
fastq_2=*_AdpTrimmed_2.fastq

#Do Fastqc check
fastqc $fastq_1
fastqc $fastq_2

#Remove low complexity reads
perl /share/pkg/prinseq/0.20.4/prinseq-lite.pl \
-fastq $fastq_1 -fastq2 $fastq_2 \
-out_good "$SAMPLE_NAME"_AdpTrimmed_goodComplex \
-out_bad "$SAMPLE_NAME"_AdpTrimmed_lowComplex \
-log \
-lc_method dust \
-lc_threshold 15

#Trim low qual sequence off of reads
#perl /share/pkg/prinseq/0.20.4/prinseq-lite.pl \
#-fastq $fastq_1 -fastq2 $fastq_2 \
#-out_good "$SAMPLE_NAME"_tagCleaned_goodQual \
#-out_bad "$SAMPLE_NAME"_tagCleaned_lowQual \
#-log \
#-trim_qual_right 20

fastq_1=*AdpTrimmed_goodComplex_1.fastq
fastq_2=*AdpTrimmed_goodComplex_2.fastq

#Do Fastqc check
fastqc $fastq_1
fastqc $fastq_2

#Align reads to Hg19 and keep unmapped reads.
bowtie2 -p $TR -x $BOWTIE2hg19DIR/genome \
--un-conc "$SAMPLE_NAME"_Hg_filt_pairs \
--un "$SAMPLE_NAME"_Hg_filt_single \
-1 $fastq_1 -2 $fastq_2 -S "$SAMPLE_NAME"_Hg_mapped.sam
/share/pkg/samtools/0.0.19/samtools view -bS "$SAMPLE_NAME"_Hg_mapped.sam > "$SAMPLE_NAME"_Hg_mapped.bam
/share/pkg/samtools/0.0.19/samtools sort "$SAMPLE_NAME"_Hg_mapped.bam "$SAMPLE_NAME"_Hg_mapped_sorted
/share/pkg/samtools/0.0.19/samtools index "$SAMPLE_NAME"_Hg_mapped_sorted.bam

fastqc "$SAMPLE_NAME"_Hg_filt_pairs*

mv "$SAMPLE_NAME"_Hg_filt_pairs.1 vic/"$SAMPLE_NAME"_Hg_filt_pairs_1.fastq
mv "$SAMPLE_NAME"_Hg_filt_pairs.2 vic/"$SAMPLE_NAME"_Hg_filt_pairs_2.fastq

cd vic/

#Denovo genome assembly using VICUNA. Make sure you have the vicuna_config.txt in the current dir and it is executable.
OMP_NUM_THREADS=32 /project/umw_jeffrey_bailey/share/bin_sync/vicuna/VICUNA_v1.3/bin/vicuna-omp-v1.0 \
./vicuna_config.txt

#Aligning vicuna assembled contigs to the reference EBV genomes to visualize on IGV
#Type1 
bowtie2 -f -x $BOWTIE2ebvT1/NC_007605 -U contig.fasta -S "$SAMPLE_NAME"_vic_type1.sam
/share/pkg/samtools/0.0.19/samtools view -bS "$SAMPLE_NAME"_vic_type1.sam > "$SAMPLE_NAME"_vic_type1.bam
/share/pkg/samtools/0.0.19/samtools sort "$SAMPLE_NAME"_vic_type1.bam "$SAMPLE_NAME"_vic_type1_sorted
/share/pkg/samtools/0.0.19/samtools index "$SAMPLE_NAME"_vic_type1_sorted.bam

#Also aligning contigs with low frequent variantions the reference for visualization on IGV
mkdir vic_lfv
mv contig.lfv.fasta vic_lfv/
cd vic_lfv/
#Aligning assembled contigs with low frequent variations to the reference
bowtie2 -f -x $BOWTIE2ebvT1/NC_007605 -U contig.lfv.fasta -S "$SAMPLE_NAME"_lfv_t1.sam
/share/pkg/samtools/0.0.19/samtools view -bS "$SAMPLE_NAME"_lfv_t1.sam > "$SAMPLE_NAME"_lfv_t1.bam
/share/pkg/samtools/0.0.19/samtools sort "$SAMPLE_NAME"_lfv_t1.bam "$SAMPLE_NAME"_lfv_t1_sorted
/share/pkg/samtools/0.0.19/samtools index "$SAMPLE_NAME"_lfv_t1_sorted.bam

#Aligning the adaptor trimmed fastq reads to the reference
cd ../
mkdir read_alignments
cp "$SAMPLE_NAME"_Hg_filt_pairs_1.fastq read_alignments/ 
cp "$SAMPLE_NAME"_Hg_filt_pairs_2.fastq read_alignments/
cd read_alignments/
ln -s /project/umw_jeffrey_bailey/co54w/DNAseq/EBV_EricsProbes/"$SAMPLE_NAME"/vic/"$SAMPLE_NAME"_Hg_filt_pairs_1.fastq .
ln -s /project/umw_jeffrey_bailey/co54w/DNAseq/EBV_EricsProbes/"$SAMPLE_NAME"/vic/"$SAMPLE_NAME"_Hg_filt_pairs_2.fastq .

#Running read alignment to EBV ref.
fastq_1="$SAMPLE_NAME"_Hg_filt_pairs_1.fastq
fastq_2="$SAMPLE_NAME"_Hg_filt_pairs_2.fastq

#Read Alignment to Type 1
bowtie2 -p $TR -x $BOWTIE2ebvT1/NC_007605 \
-q \
-1 $fastq_1 \
-2 $fastq_2 \
--rg-id 1000 --rg LB:POLYA --rg PL:ILLUMINA --rg SM:$SAMPLE_NAME \
-S "$SAMPLE_NAME"_reads_t1.sam
/share/pkg/samtools/0.0.19/samtools view -bS "$SAMPLE_NAME"_reads_t1.sam > "$SAMPLE_NAME"_reads_t1.bam
/share/pkg/samtools/0.0.19/samtools sort "$SAMPLE_NAME"_reads_t1.bam "$SAMPLE_NAME"_reads_t1_sorted
/share/pkg/samtools/0.0.19/samtools index "$SAMPLE_NAME"_reads_t1_sorted.bam


cd ../../
#Running PAGIT
source /project/umw_jeffrey_bailey/share/bin_sync/PAGIT/PAGIT/sourceme.pagit

#Running ABACAS- Orientation of contigs
#EBV TYPE1
mkdir runABACAS_T1
cd runABACAS_T1
ln -s /project/umw_jeffrey_bailey/EBV/"$SAMPLE_NAME"/vic/contig.fasta .
ln -s /project/umw_jeffrey_bailey/share/EBV/Genomes/Type1/NC_007605.fa .

perl /project/umw_jeffrey_bailey/share/bin_sync/PAGIT/PAGIT/ABACAS/abacas.pl -r NC_007605.fa -q contig.fasta -p nucmer -b -t -c -o "$SAMPLE_NAME"_T1abacas

#Concatenating together the mapped sequences and unmapped sequences. This is important so that we don't loses highly variable regions
cat "$SAMPLE_NAME"_T1abacas.fasta "$SAMPLE_NAME"_T1abacas.contigsInbin.fas > "$SAMPLE_NAME"_T1abacas_mapAndUnmap.fasta
 

#Running IMAGE- Filling the gaps between contigs using the adaptor trimmed fastq reads
#EBV TYPE1
cd ../
mkdir runIMAGE_T1
cd runIMAGE_T1
ln -s /project/umw_jeffrey_bailey/EBV/"$SAMPLE_NAME"/"$SAMPLE_NAME"_tagCleaned_1.fastq .
ln -s /project/umw_jeffrey_bailey/EBV/"$SAMPLE_NAME"/"$SAMPLE_NAME"_tagCleaned_2.fastq .
ln -s /project/umw_jeffrey_bailey/EBV/"$SAMPLE_NAME"/runABACAS_T1/"$SAMPLE_NAME"_T1abacas_mapAndUnmap.fasta

perl /project/umw_jeffrey_bailey/share/bin_sync/PAGIT/PAGIT/IMAGE/image.pl -scaffolds "$SAMPLE_NAME"_T1abacas_mapAndUnmap.fasta -prefix "$SAMPLE_NAME"_tagCleaned -iteration 1 -all_iteration 5 -dir_prefix ite -kmer 30


cd ite5

bowtie2 -f -x $BOWTIE2ebvT1/NC_007605 -U new.fa -S "$SAMPLE_NAME"_T1_IMAGE.sam
/share/pkg/samtools/0.0.19/samtools view -bS "$SAMPLE_NAME"_T1_IMAGE.sam > "$SAMPLE_NAME"_T1_IMAGE.bam
/share/pkg/samtools/0.0.19/samtools sort "$SAMPLE_NAME"_T1_IMAGE.bam "$SAMPLE_NAME"_T1_IMAGE_sorted
/share/pkg/samtools/0.0.19/samtools index "$SAMPLE_NAME"_T1_IMAGE_sorted.bam

#Filling gaps between the supercontigs withs N's
perl /project/umw_jeffrey_bailey/share/bin_sync/PAGIT/PAGIT/IMAGE/contigs2scaffolds.pl new.fa new.read.placed 2 0 "$SAMPLE_NAME"_scaffold

#Determing Contig coverage for circos plot
genomeCoverageBed -bg -ibam "$SAMPLE_NAME"_T1_IMAGE_sorted.bam -g $BOWTIE2ebvT1/NC_007605 | cut -f1-3 > "$SAMPLE_NAME"_T1_CONTIG_circos.cov

#Running ICORN
#cd ../../
#mkdir runICORN
#cd runICORN
#ln -s /project/umw_jeffrey_bailey/co54w/DNAseq/EBV_EricsProbes/"$SAMPLE_NAME"/vic_t/outReads1.fq .
#ln -s /project/umw_jeffrey_bailey/co54w/DNAseq/EBV_EricsProbes/"$SAMPLE_NAME"/vic_t/outReads2.fq .
#ln -s /project/umw_jeffrey_bailey/co54w/DNAseq/EBV_EricsProbes/"$SAMPLE_NAME"/runIMAGE2/ite10/"$SAMPLE_NAME"_scaffolds.fa .
#dir=/project/umw_jeffrey_bailey/co54w/DNAseq/EBV_EricsProbes/"$SAMPLE_NAME"/runIMAGE2/ite10

#$ICORN_HOME/icorn.start.sh "$SAMPLE_NAME"_scaffolds.fa 1 3 outReads1.fq outReads2.fq 100,500 180

#Running RATT
#Annotating our query genome for EBV type1
cd ../../
mkdir runRATT
cd runRATT
mkdir embl
cd embl
cp -s /project/umw_jeffrey_bailey/share/EBV/Annotation/NC_007605.embl .
cd ..
ln -s /project/umw_jeffrey_bailey/EBV/"$SAMPLE_NAME"/runIMAGE/ite5/"$SAMPLE_NAME"_scaffold.fa
/project/umw_jeffrey_bailey/share/bin_sync/PAGIT/PAGIT/RATT/start.ratt.sh ./embl "$SAMPLE_NAME"_scaffold.fa "$SAMPLE_NAME" Strain.Repetitive > "$SAMPLE_NAME"_ratt.output.txt
