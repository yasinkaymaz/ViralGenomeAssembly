#!/bin/bash
dir=`pwd`

#toolDir='/home/yk42w/codes/ViralGenomeAssembly'
toolDir='/Users/yasinkaymaz/Dropbox/codes/ViralGenomeAssembly'

type=1
nt=2
VariantFilterer=1
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

InputVCF=$1

if [ "$VariantFilterer" = "1" ]
then
  #do some stuff
  #Create an array of genomic positions of repetitive and non-coding regions.
  #module load bedtools/2.17.0
  bedtools intersect -header -wa -a $InputVCF -b $AnnotationDir/EBV.genes.repeatMask-regions-substracted3.bed > "${InputVCF%.vcf}"_filtered.vcf
  rm "${InputVCF%.vcf}"_filtered-OnlyNonSynVariants.vcf "${InputVCF%.vcf}"_filtered-OnlyNonSynVariants_Effects.txt
  #for each line of variants
  while read line
  do
    #if line starts with #, don't do anything other than passing line to output
    if [[ $line =~ "#" ]];
      then
      echo "$line" >> "${InputVCF%.vcf}"_filtered-OnlyNonSynVariants.vcf
      :
    elif [[ $line =~ "NC_007605" ]];
      then
    #get genomic postion of variant


      #Insert variant into gene ref seq.
      grep "#" $InputVCF > varline.vcf
      echo -e "$line" >> varline.vcf
      $IGVtoolsdir/IGVTools/igvtools index varline.vcf

      use-java 7

      java -Xmx16g -XX:ParallelGCThreads=$nt -jar \
      $GATKdir/GenomeAnalysisTK.jar \
      -R $toolDir/resources/Bowtie2Index/"${EBVgenome}".fa \
      -T FastaAlternateReferenceMaker \
      --variant varline.vcf \
      -o AlternativeGenome_with_variant.fa

      #concatinate both sequences above into multi-fasta file.
      cat $toolDir/resources/Bowtie2Index/"${EBVgenome}".fa AlternativeGenome_with_variant.fa > Ref.Variant.multi.fasta


      #which gene is the variant position in?
      #bedtools intersect -wb -a varline.vcf -b $AnnotationDir/EBV.genes.repeatMask-regions-substracted3.bed |awk '{print $(NF-2)}'
      geneOfvariant=`bedtools intersect -wb -a varline.vcf -b $AnnotationDir/EBV.genes.repeatMask-regions-substracted3.bed |awk '{print $(NF-2)}'`
      echo -e "Variant is in gene: $geneOfvariant"
      grep -w $geneOfvariant $RefGenesDir/EBV_Reference_genelist_Genenames_stranded.txt > "$geneOfvariant".bed;
      strand=`cut -f5 "$geneOfvariant".bed|uniq`;
      echo -e "Strand is $strand";

      #Extract gene sequence from Ref sequence.
      python $toolDir/bin/MSA_gene_extractor.py Ref.Variant.multi.fasta 1 "$geneOfvariant".bed;
      perl $toolDir/bin/fasta_to_tab.pl "$geneOfvariant".bed.aln.fasta > "$geneOfvariant".aln.tab;

      #is gene expressed from minus transcript? if yes, take reverse transcript of multi-fasta file.
      if [ "$strand" = "-" ]
      then
        #if gene is from reverse strand:
        while read ll;
        do
          sample_name=`echo $ll|cut -d " " -f1`;
          seq=`echo $ll|cut -d " " -f2`;
          rcseq=`echo $seq|rev | tr "ATGCNatgcn" "TACGNtacgn"`;
          echo -e "$sample_name $rcseq";
        done < "$geneOfvariant".aln.tab > "$geneOfvariant".aln.rc.tab;
        testalignmentFile="$geneOfvariant".aln.rc.tab
      #if no, keep multi-fasta as it is.
      else
        #else
        testalignmentFile="$geneOfvariant".aln.tab
      fi

      #Run SNAP on multi-fasta.
      perl $toolDir/bin/SNAP.pl $testalignmentFile;
      #Check codon.* file if it has "synon" or "nonsynon".
      VarType=`grep syn codons.*| grep -v "#"|awk '{print $2}'`
      echo -e "Variant Type is $VarType."

      #if "nonsynon", keep the variant line. Otherwise, stop. Next variant line.
      if [ "$VarType" = "nonsynon" ]
      then
        First5cols=`echo "$line"|cut -f1-5`
        VarEffect=`grep syn codons.*| grep -v "#"`
        echo "$line" >> "${InputVCF%.vcf}"_filtered-OnlyNonSynVariants.vcf
        echo -e "$First5cols\t$VarEffect\t$geneOfvariant\t$strand" >> "${InputVCF%.vcf}"_filtered-OnlyNonSynVariants_Effects.txt
      else
        echo "Variant is synonymous, passing"
      fi

      #Clean up intermediate files.
      rm summary.* codons.* background.* "$geneOfvariant".bed "$geneOfvariant".aln.tab "$geneOfvariant".bed.aln.fasta
      rm AlternativeGenome_with_variant.fa Ref.Variant.multi.fasta varline.vcf.idx

    fi
  done < "${InputVCF%.vcf}"_filtered.vcf

  rm "${InputVCF%.vcf}"_filtered.vcf

else
  echo "Not filtering any variants"
fi
