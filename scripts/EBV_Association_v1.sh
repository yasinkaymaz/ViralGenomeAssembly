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


rm -r proj_res proj_out proj

pseq proj new-project
pseq proj load-vcf --vcf input.vcf
pseq proj tag-file --id 1 --name eBLs
pseq proj tag-file --id 2 --name Controls
pseq proj load-pheno --file pop.phe
pseq proj v-assoc --phenotype phe1 --perm 10000 --vmeta > Plink.all
awk 'FNR==NR{a[$1]=$0;next}{print $0"\t", a[$1]}' Plink.all|awk '{if(NF >2)print}'|gcol VAR REF ALT OR P > Plink.covered.out.$EBVgenome.txt

module load plinkseq/0.08
module load bedtools/2.22.0

header=`pseq proj v-assoc --phenotype phe1|head -1`
echo -e "chrom\tposition\tend\t$header" > assocP.filtered.bed
pseq proj v-assoc --phenotype phe1 --fix-null|grep -v VAR|awk '{OFS="\t";gsub("chr1:","",$1); print "NC_007605""\t"$1"\t"$1+1"\t"$0}'|bedtools subtract -a - -b filter.regions.bed >> assocP.filtered.bed


header=`pseq proj v-assoc --phenotype phe1 --perm 10|head -1`
echo -e "chrom\tposition\tend\t$header" > assocP.filtered.bed
pseq proj v-assoc --phenotype phe1 --fix-null --perm 10000000|grep -v VAR|awk '{OFS="\t";gsub("chr1:","",$1); print "NC_007605""\t"$1"\t"$1+1"\t"$0}'|bedtools subtract -a - -b filter.regions.bed >> assocP.filtered.bed


##phe1,Integer,-9,1/2=eBLs/Controls
#ID     phe1
1381_22_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
1454-90_eBLFNAtumorDNA_wga_genome_fixed_assembly_Fixed  1
1583_04_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
1612-09_eBLFNAtumorDNA_wga_genome_fixed_assembly_Fixed  1
BL006_NoFNA-Plasma_PCRsWGA_genome_fixed_assembly_Fixed  1
BL008_NoFNAPaired-Plasma_PCRsWGA_genome_fixed_assembly_Fixed  1
BL033_NoFNAPaired-Plasma_PCRsWGA_genome_fixed_assembly_Fixed  1
BL035_NoFNAPaired-Plasma_PCRsWGA_genome_fixed_assembly_Fixed  1
BL044_NoFNAPaired-Plasma_PCRsWGA_genome_fixed_assembly_Fixed  1
BL099_NoFNAPaired-Plasma_genome_fixed_assembly_Fixed  1
BL103_NoFNAPaired-Plasma_genome_fixed_assembly_Fixed  1
BL108_NoFNAPaired-Plasma_genome_fixed_assembly_Fixed  1
BL120_NoFNAPaired-Plasma_genome_fixed_assembly_Fixed  1
BL204_NoFNAPaired-Plasma_genome_fixed_assembly_Fixed  1
BL234_NoFNAPaired-Plasma_PCRsWGA_genome_fixed_assembly_Fixed  1
BL247_NoFNAPaired-Plasma_PCRsWGA_genome_fixed_assembly_Fixed  1
BL251_NoFNAPaired-Plasma_genome_fixed_assembly_Fixed  1
BL392_NoFNAPaired-Plasma_genome_fixed_assembly_Fixed  1
BL534_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
BL539_eBLFNAtumorDNA_wga_genome_fixed_assembly_Fixed  1
BL541_eBLFNAtumorDNA_wga_genome_fixed_assembly_Fixed  1
BL546_eBLFNAtumorDNA_PCRsWGA_genome_fixed_assembly_Fixed  1
BL552_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
BL556_eBLFNAtumorDNA_wga_genome_fixed_assembly_Fixed  1
BL557_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
BL558_eBLFNAtumorDNA_wga_genome_fixed_assembly_Fixed  1
BL560_eBLFNAtumorDNA_wga_genome_fixed_assembly_Fixed  1
BL564_eBLFNAtumorDNA_PCRsWGA_genome_fixed_assembly_Fixed  1
BL565_eBLFNAtumorDNA_PCRsWGA_genome_fixed_assembly_Fixed  1
BL573_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
BL574_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
BL576_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
BL577_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
BL578_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
BL579_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
BL587_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
BL590_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
BL606_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
BL607_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
BL609_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
BL611_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
BL614_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
BL620_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
BL624_eBLFNAtumorDNA_PCRsWGA_genome_fixed_assembly_Fixed  1
BL627_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
BL628_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
BL629_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
BL635_eBLFNAtumorDNA_PCRsWGA_genome_fixed_assembly_Fixed  1
BL642_eBLFNAtumorDNA_wga_genome_fixed_assembly_Fixed  1
BL643_eBLFNAtumorDNA_wga_genome_fixed_assembly_Fixed  1
BL644_eBLFNAtumorDNA_wga_genome_fixed_assembly_Fixed  1
BL645_eBLFNAtumorDNA_genome_fixed_assembly_Fixed  1
BL646_eBLFNAtumorDNA_wga_genome_fixed_assembly_Fixed  1
BL647_eBLFNAtumorDNA_wga_genome_fixed_assembly_Fixed  1
BL648_eBLFNAtumorDNA_wga_genome_fixed_assembly_Fixed  1
BL717_CellLine_genome_fixed_assembly_Fixed  1
BL719_CellLine_genome_fixed_assembly_Fixed  1
BL720_CellLine_genome_fixed_assembly_Fixed  1
EC-C005-31V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
EC-C008-50V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
EC-C013-40V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
EC-C018-39V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
EC-C023-30V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
EC-C032-30V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
EC-C042-39V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
EC-C086-54V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
EC-C087-72V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
EC-C091-38V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
EC-C095-41V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
EC-C101-30V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
EC-C102-29V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
EC-C105-28V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
MC-C006-51V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
MC-C012-51V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
MC-C013-57V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
MC-C024-72V_CP_genome_fixed_assembly_Fixed  2
MC-C025-63V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
MC-C028-41V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
MC-C030-43V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
MC-C031-50V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
MC-C065-66V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
MC-C104-41V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
MC-C108-40V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
MC-C113-52V_CP_genome_fixed_assembly_Fixed  2
MC-C119-60V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
MC-C121-53V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
MC-C122-46V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
MC-C137-36V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
MC-C142-58V_NonBLkid_PCRsWGA_genome_fixed_assembly_Fixed  2
11003406_10_TTR_CP_genome_fixed_assembly_Fixed  2
1-3-0430-07-9_CP_genome_fixed_assembly_Fixed  2
1-4-0554-06-9_CP_genome_fixed_assembly_Fixed  2
1-4-0589-06-9_CP_genome_fixed_assembly_Fixed  2
1-4-0621-14-11_CP_genome_fixed_assembly_Fixed 2
1-6-0663-10-11_CP_genome_fixed_assembly_Fixed 2
1-6-0679-06-11_CP_genome_fixed_assembly_Fixed 2
16069207_9_TTR_CP_genome_fixed_assembly_Fixed 2
1-6-6007-08-11_CP_genome_fixed_assembly_Fixed 2



100795  100808  BBPL1
143511  143532  RPMS1intron



1454-90_eBLFNAtumorDNA_wga.Case 1
1583_04_eBLFNAtumorDNA.Case2    1
1612-09_eBLFNAtumorDNA_wga.Case3        1
BL006_NoFNA-Plasma_PCRsWGA.Case4        1
BL033_NoFNAPaired-Plasma_PCRsWGA.Case5  1
BL035_NoFNAPaired-Plasma_PCRsWGA.Case6  1
BL044_NoFNAPaired-Plasma_PCRsWGA.Case7  1
BL099_NoFNAPaired-Plasma.Case8  1
BL120_NoFNAPaired-Plasma.Case9  1
BL204_NoFNAPaired-Plasma.Case10 1
BL234_NoFNAPaired-Plasma_PCRsWGA.Case11 1
BL247_NoFNAPaired-Plasma_PCRsWGA.Case12 1
BL392_NoFNAPaired-Plasma.Case13 1
BL534_eBLFNAtumorDNA.Case14     1
BL546_eBLFNAtumorDNA_PCRsWGA.Case15     1
BL556_eBLFNAtumorDNA_wga.Case16 1
BL557_eBLFNAtumorDNA.Case17     1
BL558_eBLFNAtumorDNA_wga.Case18 1
BL560_eBLFNAtumorDNA_wga.Case19 1
BL573_eBLFNAtumorDNA.Case20     1
BL574_eBLFNAtumorDNA.Case21     1
BL577_eBLFNAtumorDNA.Case22     1
BL579_eBLFNAtumorDNA.Case23     1
BL590_eBLFNAtumorDNA.Case24     1
BL606_eBLFNAtumorDNA.Case25     1
BL607_eBLFNAtumorDNA.Case26     1
BL609_eBLFNAtumorDNA.Case27     1
BL614_eBLFNAtumorDNA.Case28     1
BL620_eBLFNAtumorDNA.Case29     1
BL624_eBLFNAtumorDNA_PCRsWGA.Case30     1
BL627_eBLFNAtumorDNA.Case31     1
BL628_eBLFNAtumorDNA.Case32     1
BL635_eBLFNAtumorDNA_PCRsWGA.Case33     1
BL642_eBLFNAtumorDNA_wga.Case34 1
BL644_eBLFNAtumorDNA_wga.Case35 1
BL645_eBLFNAtumorDNA.Case36     1
BL646_eBLFNAtumorDNA_wga.Case37 1
BL647_eBLFNAtumorDNA_wga.Case38 1
BL648_eBLFNAtumorDNA_wga.Case39 1
BL717_CellLine.Case40   1
BL720_CellLine.Case41   1
EC-C005-31V_NonBLkid_PCRsWGA.Control    2
EC-C013-40V_NonBLkid_PCRsWGA.Control2   2
EC-C042-39V_NonBLkid_PCRsWGA.Control3   2
EC-C086-54V_NonBLkid_PCRsWGA.Control4   2
EC-C087-72V_NonBLkid_PCRsWGA.Control5   2
EC-C091-38V_NonBLkid_PCRsWGA.Control6   2
EC-C102-29V_NonBLkid_PCRsWGA.Control7   2
MC-C006-51V_NonBLkid_PCRsWGA.Control8   2
MC-C013-57V_NonBLkid_PCRsWGA.Control9   2
MC-C030-43V_NonBLkid_PCRsWGA.Control10  2
MC-C031-50V_NonBLkid_PCRsWGA.Control11  2
MC-C108-40V_NonBLkid_PCRsWGA.Control12  2
MC-C119-60V_NonBLkid_PCRsWGA.Control13  2
MC-C142-58V_NonBLkid_PCRsWGA.Control14  2
