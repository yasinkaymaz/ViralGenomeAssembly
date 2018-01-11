#!/bin/bash
InputFasta=$1

toolDir='/Users/yasinkaymaz/Dropbox/codes/ViralGenomeAssembly'



$toolDir/bin/fasta_to_tab.pl $InputFasta |\
awk 'FNR==NR{a[$1]=$2;next}{print ">"a[$1],"\n",$2}' $toolDir/PapeRCodes/EBV_SequencingSet_ID-Key3.txt - \
> "${InputFasta%.fasta}".IDreplaced.fasta


#Alternative
#awk 'NR==FNR {a[$1]=$2;next} {for ( i in a) gsub(i,a[i])}1' EBV_SequencingSet_ID-Key.txt \
#227.Genomes.aln.fasta > 227.Genomes.aln_IDchanged.fasta
#Code from https://www.unix.com/unix-for-dummies-questions-and-answers/143054-best-method-replacing-multiple-strings-multiple-files-sed-awk-most-simple-preferred.html
#awk 'NR==FNR {a[$1]=$4;next} {for ( i in a) gsub(i,a[i])}1' key.labeled3.txt dnds2/LMP1.aln.filtered.rc.tab|awk '{print ">"$1"\n"$2}' > geneFastas/LMP1.aln.filtered.fasta
# for gene in BYRF1.aln.filtered EBNA3A.aln.filtered EBNA3B.aln.filtered EBNA3C.aln.filtered ;
# do
#   awk 'NR==FNR {a[$1]=$4;next} {for ( i in a) gsub(i,a[i])}1' key.labeled3.txt dnds2/"$gene".tab |\
#   awk '{print ">"$1"\n"$2}' > geneFastas/"$gene".fasta
# done
