#!/bin/bash

toolDir='~/codes/ViralGenomeAssembly'

#Assumes that circos is already installed and ready to be executed
InputGenomeFasta=$1

mkdir CircosPlots
cd CircosPlots
cp ~/codes/ViralGenomeAssembly/resources/circos/*conf ./
cp ~/codes/ViralGenomeAssembly/resources/circos/*txt ./
ln -s ../$InputGenomeFasta ./
#Create Contigs file from a fasta alignment file:
python ~/codes/ViralGenomeAssembly/bin/MSA2contigs_bed.py $InputGenomeFasta 1

#awk '{if(NR%2 == 0) print }' $InputGenomeFasta |sed 's/[N|-]//g'|awk '{print length}' > lens.tmp
perl ~/codes/ViralGenomeAssembly/bin/fasta_to_tab.pl $InputGenomeFasta |cut -f2|sed 's/[N|-]//g'|awk '{print length}' > lens.tmp

n=0
sr=90
head -109 circos.template.conf > circos.EBV.conf
for genomeName in `grep ">" $InputGenomeFasta|sed 's/>//g'| paste -d "\t" - lens.tmp |sort -k2,2gr|cut -f1`;
do
  let r0=89-$n;
  let r1=$r0+1;
  echo "$genomeName";
  echo -e "
<plot>
file = "${InputGenomeFasta}"."${genomeName}".contigs.bed
type = tile
layers = 20
layers_overflow = grow
layers_overflow_color = dblue
margin = 0.0u
thickness = 6
padding = 6
orientation = out
stroke_thickness = 1
stroke_color = grey
color = green
r0 = 0."$r0"r
r1 = 0."$r1"r
</plot>
"  >> circos.EBV.conf

  let n="$n"+1;
done
tail -75 circos.template.conf >> circos.EBV.conf

circos -config ./circos.EBV.conf


# <plot>
# type = histogram
# file = /home/kaymazy/EBV_Capture/Circos/EBV/1151_77_alignment_wt_log.bg
# extend_bin = no
# fill_under = yes
# fill_color = lblue
# color = blue
# thickness = 0
# r0 = 0.55r
# r1 = 0.75r
# orientation = out
# #min = 0
# #max = 1
# </plot>
#
# <plot>
# type = histogram
# file = /home/kaymazy/EBV_Capture/Circos/EBV/wt_GC_log.txt
# extend_bin = no
# fill_under = yes
# fill_color = red
# color = blue
# thickness = 0
# r0 = 0.35r
# r1 = 0.55r
# orientation = out
# </plot>
#
# <plot>
# file = EBV/1151_77_alignment_wt_snps.txt
# type = scatter
# stroke_thickness = 1
# fill_color       = grey
# stroke_color     = black
# glyph            = circle
# glyph_size       = 15
# r0 = 0.80r
# r1 = 0.85r
# #max   = 0.013
# #min   = 0
# </plot>



# # plot labels for repeats
# <plot>
# type = text
# file = repeats_label.txt
# #show_links = yes
# #link_dims = 5p,4p,8p,4p,0p
# #link_thickness = 2p
# #link_color = dgrey
# label_size = 30p
# label_font = condensed
# label_snuggle = yes
# max_snuggle_distance = 30p
# snuggle_sampling          = 2
# snuggle_tolerance         = 0.25r
# snuggle_link_overlap_test = yes
# snuggle_link_overlap_tolerance = 2p
# snuggle_refine            = yes
# padding = 0p
# rpadding = 0p
# color = black
# r0 = 0.80r
# r1 = 0.85r
# </plot>



# <plot>
# type = text
# file = EBVgenes_labelNEW.txt
# show_links = yes
# link_dims = 5p,4p,8p,4p,0p
# link_thickness = 2p
# link_color = dgrey
# label_size = 20p
# label_font = condensed
# label_snuggle = yes
# max_snuggle_distance = 30p
# snuggle_sampling          = 2
# snuggle_tolerance         = 0.25r
# snuggle_link_overlap_test = yes
# snuggle_link_overlap_tolerance = 2p
# snuggle_refine            = yes
# padding = 0p
# rpadding = 0p
# color = black
# r0 = 0.92r
# r1 = 0.98r
# </plot>
