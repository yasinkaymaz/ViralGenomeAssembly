#!/bin/bash

#scriptToRun=$1
#sampleListFile=$2
sampleListFile=$1

while read line
do

SAMPLE_NAME=`echo $line|cut -d " " -f1`
Kmer_start=`echo $line|cut -d " " -f3`
Kmer_end=`echo $line|cut -d " " -f4`
RunTime=`echo $line|cut -d " " -f5`
type=`echo $line|cut -d " " -f6`
echo $SAMPLE_NAME
echo -e "RunTime is $RunTime"

nt=4
#cd $SAMPLE_NAME
#rm "$SAMPLE_NAME"_err.*.txt "$SAMPLE_NAME"_out.*.txt
#-R "span[hosts=1]" \

bsub -q long \
-n $nt \
-R rusage[mem=20000] \
-R "select[tmp>1000]" \
-W $RunTime \
-e "$SAMPLE_NAME"_err.%J.txt \
-o "$SAMPLE_NAME"_out.%J.txt \
~/codes/ViralGenomeAssembly/scripts/AssemblyQualityChecker.sh "$SAMPLE_NAME"_genome.fa 1 2 "$SAMPLE_NAME" $type

#$scriptToRun $type "$SAMPLE_NAME"*.gatk_recal.bam "$SAMPLE_NAME"*_noMQ0_sorted.bam $SAMPLE_NAME
#~/project/OTHERS/dnaSeq_Tools/scripts/EBV_SequenceAnalysisPipeline.sh 1 $SAMPLE_NAME $Kmer_start $Kmer_end $nt n50 $type
#~/project/OTHERS/dnaSeq_Tools/scripts/EBV_SequenceAnalysisPipeline_SNPAnalysis.sh 1 $SAMPLE_NAME $nt $type
#cd -;

sleep 1;

done < $sampleListFile
