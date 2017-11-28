#!/bin/bash

toolDir='/home/yk42w/codes/ViralGenomeAssembly'
module unload python/2.7.5
module load python/2.7.9
module load python/2.7.9_packages/pandas/0.17.1
#for genome in `grep ">" sequences.aln.fasta|grep -v NC_|sed 's/>//g'`;
# for genome in Jijoye_CP_genome_fixed_assembly_Fixed Jijoye_J100_PCRsWGA_genome_fixed_assembly_Fixed LN827800.1_Jijoye;
# do
#   echo $genome;
#   python $toolDir/bin/MSA_parser_cleaner.py $genome NC_007605 sequences.aln.fasta 1;
#   python $toolDir/bin/MSA_parser_cleaner.py $genome NC_009334 sequences.aln.fasta 2;
# done
#
# for genome in Daudi_CellLine_genome_fixed_assembly_Fixed Daudi_D100_PCRsWGA_genome_fixed_assembly_Fixed;
# do
#   python $toolDir/bin/MSA_parser_cleaner.py $genome LN827545.1_Daudi sequences.aln.fasta 1;
# done
#
# for genome in Jijoye_CellLine_genome_fixed_assembly_Fixed Jijoye_CP_genome_fixed_assembly_Fixed Jijoye_J100_PCRsWGA_genome_fixed_assembly_Fixed;
# do
#   python $toolDir/bin/MSA_parser_cleaner.py $genome LN827800.1_Jijoye sequences.aln.fasta 2;
# done

python $toolDir/bin/MSA_parser_cleaner.py LN827800.1_Jijoye Jijoye_CellLine_genome_fixed_assembly_Fixed sequences.aln.fasta 2;
python $toolDir/bin/MSA_parser_cleaner.py LN827800.1_Jijoye Jijoye_J100_PCRsWGA_genome_fixed_assembly_Fixed sequences.aln.fasta 2;
python $toolDir/bin/MSA_parser_cleaner.py LN827800.1_Jijoye Jijoye_CP_genome_fixed_assembly_Fixed sequences.aln.fasta 2;
python $toolDir/bin/MSA_parser_cleaner.py LN827545.1_Daudi Daudi_D100_PCRsWGA_genome_fixed_assembly_Fixed sequences.aln.fasta 1;
python $toolDir/bin/MSA_parser_cleaner.py LN827545.1_Daudi Daudi_CellLine_genome_fixed_assembly_Fixed sequences.aln.fasta 1;
python $toolDir/bin/MSA_parser_cleaner.py KF717093.1_Raji Raji_CellLine_longRead_genome_fixed_assembly_Fixed sequences.aln.fasta 1;
python $toolDir/bin/MSA_parser_cleaner.py KF717093.1_Raji Raji_Rep1_PCRsWGA_genome_fixed_assembly_Fixed sequences.aln.fasta 1;
python $toolDir/bin/MSA_parser_cleaner.py KF717093.1_Raji Raji_Rep2_PCRsWGA_genome_fixed_assembly_Fixed sequences.aln.fasta 1;
python $toolDir/bin/MSA_parser_cleaner.py KF717093.1_Raji Raji_CellLine_shortRead_genome_fixed_assembly_Fixed sequences.aln.fasta 1;
python $toolDir/bin/MSA_parser_cleaner.py KF717093.1_Raji Raji_GenomiPhi_genome_fixed_assembly_Fixed sequences.aln.fasta 1;
python $toolDir/bin/MSA_parser_cleaner.py BL573_eBLFNAtumorDNA_genome_fixed_assembly_Fixed BL573_eBLFNAtumorDNA_wga_genome_fixed_assembly_Fixed sequences.aln.fasta 1;
python $toolDir/bin/MSA_parser_cleaner.py BL576_eBLFNAtumorDNA_genome_fixed_assembly_Fixed BL576_eBLFNAtumorDNA_wga_genome_fixed_assembly_Fixed sequences.aln.fasta 2;
python $toolDir/bin/MSA_parser_cleaner.py BL578_eBLFNAtumorDNA_genome_fixed_assembly_Fixed BL578_eBLFNAtumorDNA_wga_genome_fixed_assembly_Fixed sequences.aln.fasta 2;
python $toolDir/bin/MSA_parser_cleaner.py BL590_eBLFNAtumorDNA_genome_fixed_assembly_Fixed BL590_eBLFNAtumorDNA_wga_genome_fixed_assembly_Fixed sequences.aln.fasta 1;
python $toolDir/bin/MSA_parser_cleaner.py BL607_eBLFNAtumorDNA_genome_fixed_assembly_Fixed BL607_eBLFNAtumorDNA_wga_genome_fixed_assembly_Fixed sequences.aln.fasta 1;
python $toolDir/bin/MSA_parser_cleaner.py KF717093.1_Raji AB828191.1_LGY-Raji sequences.aln.fasta 1;
