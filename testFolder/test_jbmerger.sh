

type=$1
ReadBamfile=$2
ContigsBamFile=$3
SAMPLE_NAME=$4

toolDir='/home/yk42w/codes/ViralGenomeAssembly'
module unload openssl/1.0.1g
module load samtools/0.0.19
module load R/3.2.2
module load python/2.7.5
module load python/2.7.5_packages/biopython/1.68

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


#Merge overlapping contigs
python $toolDir/bin/jbam_tools_v03.py \
ref_merge_contigs \
-n $EBVgenome \
-t 12 \
--regions_to_mask $AnnotationDir/"$EBVgenome"_repeatMask.bed \
-r $ReadBamfile \
-c $ContigsBamFile \
-u 10 \
-f $toolDir/resources/Bowtie2Index/"$EBVgenome".fa \
-g 1 \
--min_unique_bp 60 \
-o "$SAMPLE_NAME"_EBV_"$Type"_JBv3merged

#Map assembly contigs and index bam file.
$toolDir/scripts/map_and_index_genomeFasta.sh $type "$SAMPLE_NAME"_EBV_"$Type"_JBv3merged_genome.fa $SAMPLE_NAME
