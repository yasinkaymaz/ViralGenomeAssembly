#!/bin/bash

inputvcfFile=$1
rm "${inputvcfFile%.vcf}"_outfile.vcf

while read line
do
  #if line starts with #, don't do anything other than passing line to output
  if [[ $line =~ "#" ]];
    then
    echo "$line" >> "${inputvcfFile%.vcf}"_outfile.vcf
    :
  elif [[ $line =~ "*" ]];
    then
    First5cols=`echo "$line"|cut -f1-5|sed 's/*/N/g'`

  #  alt=`echo "$line"|cut -f6`
    index=`echo "$line"|cut -f5|awk '{ gsub(/\,/,"",$1); print index($1, "*") }'`
    fixed6toRestcols=`echo "$line"|cut -f6-|sed 's/'"$index"'/./g'`

    echo -e "$index"
    echo "$line"
    echo -e "$First5cols\t$fixed6toRestcols"
    echo -e "$First5cols\t$fixed6toRestcols" >> "${inputvcfFile%.vcf}"_outfile.vcf
  else
    echo "$line" >> "${inputvcfFile%.vcf}"_outfile.vcf

  fi


done < $inputvcfFile
