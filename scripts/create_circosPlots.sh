#!/bin/bash

toolDir='/home/yk42w/codes/ViralGenomeAssembly'

#Assumes that circos is already installed and ready to be executed

#export PERL5LIB=/home/baileyj1/share/bin_sync/perl5lib/lib/perl5/


/home/baileyj1/share/bin_sync/circos-0.69-3/bin/circos -config ~/results/DNAseq/EBV/Capture_seq/Batch3/Circos_workspace/CircosPlots/configureCovNew.conf
/home/baileyj1/share/bin_sync/circos-0.69-3/bin/circos -config ~/results/DNAseq/EBV/Capture_seq/Batch3/Circos_workspace/CircosPlots/configureCovNew2.conf

/home/baileyj1/share/bin_sync/circos-0.69-3/bin/circos -config ~/results/DNAseq/EBV/Capture_seq/Batch3/Circos_workspace/CircosPlots/configureContigCov.conf
