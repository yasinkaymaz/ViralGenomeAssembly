#!/bin/bash

FASTAfile=$1

/project/umw_jeffrey_bailey/share/bin_sync/kalign/kalign \
-distance pair \
-tree nj \
-f clu \
$FASTAfile \
${FASTAfile%.fasta}.aln \
> ${FASTAfile%.fasta}.kout

#This is how to run this script: Fasta file has to end with .fasta
#bsub -q long -n 6 -R rusage[mem=20000] -R "span[hosts=1]" -R "select[tmp>1000]" -W 240:00 -e err.%J.txt -o out.%J.txt ~/codes/EBVseq/msa.sh All.Jijoye.fasta
