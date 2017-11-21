#!/bin/bash
#BSUB -n 4
#BSUB -R rusage[mem=10000]
#BSUB -R "span[hosts=1]"
#BSUB -q short
#BSUB -W 04:00
#BSUB -o "/project/umw_jeffrey_bailey/yk42w/std_out/align_out.%J.out"
#BSUB -e "/project/umw_jeffrey_bailey/yk42w/std_err/align_err.%J.err"

source /home/yk42w/project/share/bin_sync/mugsy_x86-64-v1r2.3/mugsyenv.sh
#cd to dir of genomes with separate fasta files.

DIR=`pwd`

mugsy --directory ./ --prefix mugsy_align *.fa


