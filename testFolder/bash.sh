

grep -v ">" Zta_promoter.bed.aln.rc.fasta|cut -c39,46,47,111,206,297,430,465,471 > ~/Dropbox/Papers/EBV_project/workspace/data/Z_promoter.variants
grep ">" Zta_promoter.bed.aln.rc.fasta|paste -d "\t" - ~/Dropbox/Papers/EBV_project/workspace/data/Z_promoter.variants > ~/Dropbox/Papers/EBV_project/workspace/data/Z_promoter.variants.txt

sed 's/>//g' ~/Dropbox/Papers/EBV_project/workspace/data/Z_promoter.variants.txt| awk 'NR==FNR{a[$1]=$0;next}{print a[$1]"\t"$0}' ~/Dropbox/codes/ViralGenomeAssembly/workspace/data/pop_info.txt - > ~/Dropbox/Papers/EBV_project/workspace/data/Z_promoter.variants.annotated.txt

grep Type1 ~/Dropbox/Papers/EBV_project/workspace/data/Z_promoter.variants.annotated.txt|cut -f6,7|awk '{print ">"$1"\n"$2}' > type1.Z_promoter.variants.fasta
grep Type2 ~/Dropbox/Papers/EBV_project/workspace/data/Z_promoter.variants.annotated.txt|cut -f6,7|awk '{print ">"$1"\n"$2}' > type2.Z_promoter.variants.fasta
