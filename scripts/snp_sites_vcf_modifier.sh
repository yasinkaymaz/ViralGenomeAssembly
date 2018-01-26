#!/bin/bash
#An alignment file is an input.
InputFastaAlignment=$1
inputvcfFile=${InputFastaAlignment%.fasta}.vcf

rm "${inputvcfFile%.vcf}"_modified.vcf

#snp-sites needs to be installed.
snp-sites -v -o $inputvcfFile $InputFastaAlignment

while read line
do
  #if line starts with #, don't do anything other than passing line to output
  if [[ $line =~ "#" ]];
    then
    echo "$line" >> "${inputvcfFile%.vcf}"_modified.vcf
    :
  elif [[ $line =~ "*" ]];
    then
    First5cols=`echo "$line"|cut -f1-5|sed 's/*/N/g'|sed 's/N,//g'|sed 's/,N//g'`

  #  alt=`echo "$line"|cut -f6`
    index=`echo "$line"|cut -f5|awk '{ gsub(/\,/,"",$1); print index($1, "*") }'`
    fixed6toRestcols=`echo "$line"|cut -f6-|sed 's/'"$index"'/./g'`

    echo -e "$index"
    echo "$line"
    echo -e "$First5cols\t$fixed6toRestcols"
    echo -e "$First5cols\t$fixed6toRestcols" >> "${inputvcfFile%.vcf}"_modified.vcf
  else
    echo "$line" >> "${inputvcfFile%.vcf}"_modified.vcf

  fi


done < $inputvcfFile
