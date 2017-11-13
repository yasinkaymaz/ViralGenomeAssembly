#!/bin/bash

#Input file has to end with .fasta
FASTAfile=$1
TR=$2

/project/umw_jeffrey_bailey/share/bin_sync/mafft-7.215-with-extensions/scripts/mafft \
--thread $TR \
--auto $FASTAfile > ${FASTAfile%.fasta}.aln.fasta

#This is how to run this script: Fasta file has to end with .fasta
#bsub -q long -n 24 -R rusage[mem=12000] -R "span[hosts=1]" -R "select[tmp>1000]" -W 240:00 -e err.%J.txt -o out.%J.txt ~/codes/EBVseq/msa_mafft.sh All.Jijoye.fasta TR
