#!/bin/bash

#This is how to run this script: input list files have to end with .list   
#bsub -q long -n 6 -R rusage[mem=20000] -R "span[hosts=1]" -R "select[tmp>1000]" -W 240:00 -e err.%J.txt -o out.%J.txt ~/codes/EBVseq/EBV_Association_v1.sh List_of_Case_Vcfs.list List_of_Control_Vcfs.list EBVTYPE[1/2]		  

List_of_Case_Vcfs=$1
List_of_Control_Vcfs=$2
#Choose the type of the EBV genome!
type=$3
#type=2
nt=2

WORKINGDIR=`pwd`

mergeVars=1
PLINK=1

module load java/1.7.0_25
module load samtools/0.0.19
module load tophat/2.0.9
module load bowtie2/2-2.1.0
module load R/3.0.1
module load tabix/0.2.6
module load bedtools/2.22.0
module load plinkseq/0.08
export PATH=/project/umw_jeffrey_bailey/share/bin_sync/picard-tools-1.111/:$PATH
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
snpEffDIR="/project/umw_jeffrey_bailey/share/bin_sync/snpEff/snpEff"
#snpEffDIR="/project/umw_jeffrey_bailey/share/bin_sync/snpEff/"
BCFtoolsDIR="/project/umw_jeffrey_bailey/share/bin_sync/samtools-0.1.19/bcftools"
PICARDPATH="/project/umw_jeffrey_bailey/share/bin_sync/picard-tools-1.111"
BOWTIE2hg19DIR="/project/umw_jeffrey_bailey/share/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index"
GATKdir="/project/umw_jeffrey_bailey/share/bin_sync/"
GATKbundleDIR="/project/umw_jeffrey_bailey/share/Homo_sapiens/GATK_bundle"
#RESULTSDIR="/project/umw_jeffrey_bailey/yk42w/results/RNAseq/PolyA/Hg19_Transcriptome"
RESULTSDIR="/project/umw_jeffrey_bailey/yk42w/results/RNAseq/PolyA/Hg19_Transcriptome/Subsample"

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


#mkdir $WORKINGDIR/Plink
#cd  $WORKINGDIR/Plink

#-------------------------
if [ "$mergeVars" = "1" ]
then

java -Xmx16g -jar \
$GATKdir/GenomeAnalysisTK.jar \
-T CombineVariants \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
-V:Case $List_of_Case_Vcfs \
-nt $nt \
-genotypeMergeOptions UNIQUIFY \
-o Cases_combined_$EBVgenome.vcf

java -Xmx16g -jar \
$GATKdir/GenomeAnalysisTK.jar \
-T CombineVariants \
-R $INDEXGENOMESDIR/$EBVgenome.fa \
-V:Control $List_of_Control_Vcfs \
-nt $nt \
-genotypeMergeOptions UNIQUIFY \
-o Controls_combined_$EBVgenome.vcf

rm pool.bed

for i in `cat $List_of_Case_Vcfs $List_of_Control_Vcfs `;
do
name=`echo $i|sed 's/_a/\t/g'|cut -f1`;
echo $name;
awk '{if($3 > 1) print $1":"$2}' "$name"_genomecov.bed >> pool.bed;
done
totalsamples=`cat $List_of_Case_Vcfs $List_of_Control_Vcfs|wc -l`
echo "VAR" > covered.regions
sort pool.bed|uniq -c|awk 'x="'$totalsamples'"{if($1 == x ) print $2}' >> covered.regions

echo -e "##phe1,Integer,-9,"1/2=eBLs/Controls"" > pop.phe
echo -e "#ID\tphe1" >> pop.phe

grep "#CHROM" Cases_combined_$EBVgenome.vcf|sed 's/\t/\n/g'|grep Case|awk '{print $1"\t"1}' >> pop.phe
#grep "#CHROM" Cases_combined_$EBVgenome.vcf|sed 's/\t/\n/g'|tail -n+10|awk '{print $1"\t"1}' >> pop.phe
grep "#CHROM" Controls_combined_$EBVgenome.vcf|sed 's/\t/\n/g'|grep Control|awk '{print $1"\t"2}' >> pop.phe
#grep "#CHROM" Controls_combined_$EBVgenome.vcf|sed 's/\t/\n/g'|tail -n+10|awk '{print $1"\t"2}' >> pop.phe

else
echo "Skipping merge step!"
fi


if [ "$PLINK" = "1" ]
then

rm -r proj_res proj_out proj

pseq proj new-project
pseq proj load-vcf --vcf Cases_combined_$EBVgenome.vcf
pseq proj load-vcf --vcf Controls_combined_$EBVgenome.vcf
pseq proj var-summary
pseq proj tag-file --id 1 --name eBLs
pseq proj tag-file --id 2 --name Controls
pseq proj var-summary
#pseq proj v-view |head
#pseq proj v-view --samples --vmeta|head
pseq proj load-pheno --file pop.phe
pseq proj v-view --geno --phenotype phe1 --mask limit=1
pseq proj i-stats 

#pseq proj counts --phenotype phe1 --annotate refseq --mask reg=chr22 #important
pseq proj v-freq --mask hwe=0:1e-7 | gcol VAR HWE 
pseq proj v-dist --phenotype phe1
pseq proj counts --phenotype phe1 --mask case.control=8-10,0-1
pseq proj v-assoc --phenotype phe1 | gcol MINA MINU OBSA OBSU VAR OR P 
pseq proj v-assoc --phenotype phe1 --perm 10000 --vmeta > Plink.all.out.$EBVgenome.txt
awk 'FNR==NR{a[$1]=$0;next}{print $0"\t", a[$1]}' Plink.all.out.$EBVgenome.txt covered.regions|awk '{if(NF >2)print}'|gcol VAR REF ALT OR P > Plink.covered.out.$EBVgenome.txt

#pseq proj v-assoc --phenotype phe1 --perm 10000 --vmeta | gcol MINA MINU OBSA OBSU VAR OR P |awk '($7 < 0.05)' > Plink.out.Psig.$EBVgenome.txt
#pseq proj v-freq --phenotype phe1

else
echo "No plink"
fi

