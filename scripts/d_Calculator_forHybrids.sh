#!/bin/bash
dir=`pwd`

#toolDir='/home/yk42w/codes/ViralGenomeAssembly'
toolDir='/Users/yasinkaymaz/Dropbox/codes/ViralGenomeAssembly'

type=1
nt=2
DistanceCalc=1
#EBV genome type
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

use-java () {
  export JAVA_HOME=`/usr/libexec/java_home -v 1.$1`
}


GATKdir="$toolDir/bin/mac"
IGVtoolsdir="$toolDir/bin/mac"
#module load java/1.7.0_25
#module load IGVTools/2.3.31

#This input genome MSA fasta file has to have Ref genome.
GenomeMSAfasta=$1

if [ "$DistanceCalc" = "1" ]
then
  #do some stuff
  #Create an array of genomic positions of repetitive and non-coding regions.
  #module load bedtools/2.17.0
  rm *pairwisedistances*
  rm *pairwiseMeandNdS*
  #for each coding gene
  for geneOfvariant in $(cut -f4 $AnnotationDir/EBV_Reference_genelist_Genenames_stranded.txt|sort|uniq);
  do
    echo $geneOfvariant;
    grep -w $geneOfvariant $AnnotationDir/EBV_Reference_genelist_Genenames_stranded.txt > "${geneOfvariant}".bed;
    strand=`cut -f5 "$geneOfvariant".bed|uniq`;
    echo -e "Strand is $strand";

    #Extract gene sequence from Ref sequence.
    python $toolDir/bin/MSA_gene_extractor.py $GenomeMSAfasta 1 "$geneOfvariant".bed;

    #<--->
    perl $toolDir/bin/fasta_to_tab.pl "$geneOfvariant".bed.aln.fasta > "$geneOfvariant".aln.tab;

    #is gene expressed from minus transcript? if yes, take reverse transcript of multi-fasta file.
    if [ "$strand" = "-" ]
    then
      #if gene is from reverse strand:
      while read ll;
      do
        sample_name=`echo $ll|cut -d " " -f1`;
        seq=`echo $ll|cut -d " " -f2`;
        rcseq=`echo $seq|rev | tr "ATGCNatgcn" "TACGNTACGN"`;
        echo -e "$sample_name\t$rcseq";
      done < "$geneOfvariant".aln.tab |grep -v "NC_007605" > "$geneOfvariant".aln.rc.tab;

      awk '{print ">"$1"\n"$2 }' "$geneOfvariant".aln.rc.tab > "$geneOfvariant".aln.rc.fasta


  #    grep -w Type1 $toolDir/workspace/data/pop_info.txt | awk 'NR==FNR{a[$1]=$0;next}{print a[$1]}' "$geneOfvariant".aln.rc.tab - | awk '{print ">"$1"\n"$2 }' > "$geneOfvariant".aln.rc.Type1.fasta;
  #    Type1multifasta="$geneOfvariant".aln.rc.Type1.fasta;
  #    grep -w Type2 $toolDir/workspace/data/pop_info.txt | awk 'NR==FNR{a[$1]=$0;next}{print a[$1]}' "$geneOfvariant".aln.rc.tab - | awk '{print ">"$1"\n"$2 }' > "$geneOfvariant".aln.rc.Type2.fasta;
  #    Type2multifasta="$geneOfvariant".aln.rc.Type2.fasta;

      genemultifasta="$geneOfvariant".aln.rc.fasta

    #if no, keep multi-fasta as it is. just convert to upper case
    else

      while read ll;
      do
        sample_name=`echo $ll|cut -d " " -f1`;
        seq=`echo $ll|cut -d " " -f2`;
        rcseq=`echo $seq| tr "ATGCNatgcn" "ATGCNATGCN"`;
        echo -e "$sample_name\t$rcseq";
      done < "$geneOfvariant".aln.tab |grep -v "NC_007605" > "$geneOfvariant".aln.c.tab

      awk '{print ">"$1"\n"$2 }' "$geneOfvariant".aln.c.tab > "$geneOfvariant".aln.fasta

  #    grep -w Type1 $toolDir/workspace/data/pop_info.txt | awk 'NR==FNR{a[$1]=$0;next}{print a[$1]}' "$geneOfvariant".aln.c.tab - | awk '{print ">"$1"\n"$2 }' > "$geneOfvariant".aln.c.Type1.fasta;
  #    Type1multifasta="$geneOfvariant".aln.c.Type1.fasta;
  #    grep -w Type2 $toolDir/workspace/data/pop_info.txt | awk 'NR==FNR{a[$1]=$0;next}{print a[$1]}' "$geneOfvariant".aln.c.tab - | awk '{print ">"$1"\n"$2 }' > "$geneOfvariant".aln.c.Type2.fasta;
  #    Type2multifasta="$geneOfvariant".aln.c.Type2.fasta;

      genemultifasta="$geneOfvariant".aln.fasta

    fi
    #<--->

  #  bash $toolDir/scripts/DNA2Protein.sh $genemultifasta

    #Calculate d-Kimura two parameter for All genomes
    python $toolDir/bin/MSA_distanceCalc.py MeanPairwiseKimuraDist -af $genemultifasta -rn $geneOfvariant >> pairwisedistances.txt
#    python $toolDir/bin/MSA_distanceCalc.py MeandNdS -af $genemultifasta -rn $geneOfvariant >> pairwiseMeandNdS.txt
    ##Calculate d-Kimura two parameter for Type1 genomes
    # python $toolDir/bin/MSA_distanceCalc.py MeanPairwiseKimuraDist -af $Type1multifasta -rn $geneOfvariant >> Type1.pairwisedistances.txt
    # python $toolDir/bin/MSA_distanceCalc.py MeandNdS -af $Type1multifasta -rn $geneOfvariant >> Type1.pairwiseMeandNdS.txt
    # #Calculate d-Kimura two parameter for Type2 genomes
    # python $toolDir/bin/MSA_distanceCalc.py MeanPairwiseKimuraDist -af $Type2multifasta -rn $geneOfvariant >> Type2.pairwisedistances.txt
    # python $toolDir/bin/MSA_distanceCalc.py MeandNdS -af $Type2multifasta -rn $geneOfvariant >> Type2.pairwiseMeandNdS.txt

    #Clean up intermediate files.
#    rm "$geneOfvariant".bed "$geneOfvariant".aln.tab "$geneOfvariant".bed.aln.fasta $genemultifasta
    #rm "$geneOfvariant".bed "$geneOfvariant".aln.tab "$geneOfvariant".bed.aln.fasta $genemultifasta $Type1multifasta $Type2multifasta
#    rm "$geneOfvariant".aln.c.tab "$geneOfvariant".aln.rc.tab
  done

  awk 'NR==FNR{a[$2]=$0;next}{print a[$1]"\t"$0}' $toolDir/workspace/data/EBV_gene_order.txt pairwisedistances.txt|cut -f1,3,5-|sort -k1,1n > pairwisedistances.ordered.txt
  # awk 'NR==FNR{a[$2]=$0;next}{print a[$1]"\t"$0}' $toolDir/workspace/data/EBV_gene_order.txt Type1.pairwisedistances.txt|cut -f1,3,5-|sort -k1,1n > Type1.pairwisedistances.ordered.txt
  # awk 'NR==FNR{a[$2]=$0;next}{print a[$1]"\t"$0}' $toolDir/workspace/data/EBV_gene_order.txt Type2.pairwisedistances.txt|cut -f1,3,5-|sort -k1,1n > Type2.pairwisedistances.ordered.txt

#  awk 'NR==FNR{a[$2]=$0;next}{print a[$1]"\t"$0}' $toolDir/workspace/data/EBV_gene_order.txt pairwiseMeandNdS.txt|cut -f1,3,5-|sort -k1,1n > pairwiseMeandNdS.ordered.txt
  # awk 'NR==FNR{a[$2]=$0;next}{print a[$1]"\t"$0}' $toolDir/workspace/data/EBV_gene_order.txt Type1.pairwiseMeandNdS.txt|cut -f1,3,5-|sort -k1,1n > Type1.pairwiseMeandNdS.ordered.txt
  # awk 'NR==FNR{a[$2]=$0;next}{print a[$1]"\t"$0}' $toolDir/workspace/data/EBV_gene_order.txt Type2.pairwiseMeandNdS.txt|cut -f1,3,5-|sort -k1,1n > Type2.pairwiseMeandNdS.ordered.txt

else
  echo "skip gene extracting..."
fi
