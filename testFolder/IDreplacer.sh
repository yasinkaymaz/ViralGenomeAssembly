#!/bin/bash

awk 'NR==FNR {a[$1]=$2;next} {for ( i in a) gsub(i,a[i])}' EBV_SequencingSet_ID-Key.txt \
227.Genomes.aln.fasta > 227.Genomes.aln_IDchanged.fasta
#Code from https://www.unix.com/unix-for-dummies-questions-and-answers/143054-best-method-replacing-multiple-strings-multiple-files-sed-awk-most-simple-preferred.html
