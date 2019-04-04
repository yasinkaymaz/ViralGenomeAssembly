#!/bin/bash


#take patch sequences
#seqName  str end seq

#take MSA fasta in tab format
# fasta_to_tab.pl
#
# for i in length(seq):
#   if i is in patch seq coordinates:
#     replace the base with patch[i]
#   else
#     base=genome[i]
#
toolDir='/Users/yasinkaymaz/Dropbox/codes/ViralGenomeAssembly'

type=1
nt=2

if [ "$type" = "1" ]
then
echo "Genome is Type I"
Type='type1'
EBVgenome="NC_007605"
VfatDIR="$toolDir/resources/Vfat"
RefDIR="$toolDir/resources/Annotation/Type1/Vfat/EBV"
AnnotationDir="$toolDir/resources/Annotation/Type1/"
RefGenesDir="$toolDir/resources/Annotation/Type1/"
else
echo "Genome is Type II"
Type='type2'
EBVgenome="NC_009334"
VfatDIR="$toolDir/resources/Vfat/Type2"
RefDIR="$toolDir/resources/Annotation/Type2/Vfat/EBV"
AnnotationDir="$toolDir/resources/Annotation/Type2/"
fi

genomeMSAfasta=$1
patchfasta=$2

perl $toolDir/bin/fasta_to_tab.pl "$genomeMSAfasta" > "${genomeMSAfasta%.fasta}".tab;
perl $toolDir/bin/fasta_to_tab.pl "$patchfasta" > "${patchfasta%.fasta}".tab;




for genome in eBL-Tumor-0012 eBL-Plasma-0049 eBL-Tumor-0033 Namalwa_new Jijoye_new;
do
  python ~/codes/ViralGenomeAssembly/bin/SimilarityPlotter.py $genome all.samples.aln.patched.extended.c.fasta
done
