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

#Take input dna multi fasta file. this file hast to end with .fasta
InputEpitopesFile=$1
DNAsDir='/Users/yasinkaymaz/Documents/EBV/ourbatch/CorrectFastas/backup'
rm Epitopes_Table.tab

while read line;
do
  if [[ $line =~ "#" ]];
    then
    echo "Started..."
    :
  else
    EpitopeType=`echo $line|cut -d " " -f1`;
    gene=`echo $line|cut -d " " -f2`;
    antigen=`echo $line|cut -d " " -f3`;
    start=`echo $line|cut -d " " -f4`;
    end=`echo $line|cut -d " " -f5`;
    CommonSeq=`echo $line|cut -d " " -f6`;
    echo "$antigen"_"$start"_"$end";

    weblogo -A protein -c chemistry -l $start -u $end -D fasta -F pdf \
    --fineprint "$antigen"_p."$start"-"$end"_"$EpitopeType"_"$CommonSeq" \
    --title-fontsize 6 --resolution 400 --yaxis 5 \
    -f $DNAsDir/my_ICed_"$gene".aln.filtered.protein.fasta -o "$antigen"_"$start"-"$end".pdf;
  fi

  let epitopeLen=$((end-start+1));

  perl $toolDir/bin/fasta_to_tab.pl $DNAsDir/my_ICed_"$gene".aln.filtered.protein.fasta |\
  awk -v start="$start" -v len="$epitopeLen" '{print $1"\t"substr($2,start,len) }' |\
  awk -v antigen="$antigen" -v range="$start"_"$end" -v epitopetype="$EpitopeType" -v seq="$CommonSeq" 'FNR==NR{a[$1]=$2;next}{OFS="\t"; print antigen,range,epitopetype,seq,a[$1],$2}' $toolDir/PapeRCodes/EBV_SequencingSet_ID-Key3.txt - >> Epitopes_Table.tab

done < $InputEpitopesFile
#Annotate table with sample info
echo -e "SampleID\tViralSubtype\tSampleSource\tSampleGroup\tAntigenName\tAntigenRange\tEpitopeType\tEpitopeSequence\tSampleName\tSampleEpitopeSequence" > Epitopes_Table.txt
awk 'NR==FNR{a[$1]=$0;next}{print a[$5]"\t"$0}' $toolDir/PapeRCodes/EBV_SequencingSet_Labels.txt Epitopes_Table.tab >> Epitopes_Table.txt
rm Epitopes_Table.tab

#
# while read line;
# do
#   if [[ $line =~ "#" ]];
#     then
#     echo "Started..."
#     :
#   else
#     EpitopeType=`echo $line|cut -d " " -f1`;
#     gene=`echo $line|cut -d " " -f2`;
#     antigen=`echo $line|cut -d " " -f3`;
#     start=`echo $line|cut -d " " -f4`;
#     end=`echo $line|cut -d " " -f5`;
#     CommonSeq=`echo $line|cut -d " " -f6`;
#     echo "$antigen"_"$start"_"$end";
#
#     #for each Sample group: eBL HealthyControl
#     for group in eBL HealthyControl;
#     do
#       echo $group
#       #Subset with both types:
#       awk -v seq="$CommonSeq" -v group="$group" '{if( ($4 == group) && ($8 == seq) ) print ">"$9"\n"$10 }' \
#       Epitopes_Table.txt > "$group"_"$antigen"_"$start"_"$end".fasta
#
#       weblogo -A protein -c chemistry -D fasta -F pdf \
#       --fineprint "$group"_"$antigen"_p."$start"-"$end"_"$EpitopeType"_"$CommonSeq" \
#       --title-fontsize 6 --resolution 400 --yaxis 5 \
#       -f "$group"_"$antigen"_"$start"_"$end".fasta -o "$group"_"$antigen"_"$start"-"$end".pdf;
#
#       rm "$group"_"$antigen"_"$start"_"$end".fasta
#
#       #for each EBVsubtype: Type1 Type2
#       for type in Type1 Type2;
#       do
#         echo $group $type;
#
#         #subset the Epitopes table:
#         awk -v seq="$CommonSeq" -v group="$group" -v type="$type" '{if( ($2 == type) && ($4 == group) && ($8 == seq) ) print ">"$9"\n"$10 }' \
#         Epitopes_Table.txt > "$group"_"$type"_"$antigen"_"$start"_"$end".fasta
#
#         weblogo -A protein -c chemistry -D fasta -F pdf \
#         --fineprint "$group"_"$type"_"$antigen"_p."$start"-"$end"_"$EpitopeType"_"$CommonSeq" \
#         --title-fontsize 6 --resolution 400 --yaxis 5 \
#         -f "$group"_"$type"_"$antigen"_"$start"_"$end".fasta -o "$group"_"$type"_"$antigen"_"$start"-"$end".pdf;
#
#         rm "$group"_"$type"_"$antigen"_"$start"_"$end".fasta
#       done
#     done
#
#   fi
#
# done < $InputEpitopesFile

#How to run this script:
#bash ~/Dropbox/codes/ViralGenomeAssembly/scripts/EpitopeProfiler.sh Epitope.data_v4.txt
