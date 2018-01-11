#!/usr/bin/python
"""
This code is a collaborative effort with Sevtap Duman Northeastern University, systems security lab.
"""
import os
import pandas  as pd
import sys
from Bio import AlignIO
import Bio.Align
dir = os.path.dirname(__file__)
print os.getcwd()
print dir
alndata = []
#sys.argv[1] -> sample 1
#sys.argv[2] -> sample 2
#sys.argv[3] -> Alignment File
#sys.argv[4] -> type of virus
#
if len(sys.argv) < 3:
	print "This code computes the similarities between a given genome and both type 1 and type 2 reference EBV genomes separately using an multiple sequence alignment file."
	print "Please provide required arguments in proper order:"
	print "MSA_parser_cleaner.py Sample1 AlignmentFile"
	print "Feed with an alignment file in fasta format."
	print "The alignment should also include two reference EBV genomes, NC_007605 and NC_009334."
	sys.exit(1)
# Ref=''
# #Make sure that repeat files are in the same directory
# if sys.argv[4] == '1':
Ref1='NC_007605'
Repeat_inputFile = '%(directory)s/../resources/Annotation/Type1/NC_007605_repeatMask.bed' %{"directory":dir}
# else:
Ref2='NC_009334'
# 	Repeat_inputFile = '%(directory)s/../resources/Annotation/Type2/NC_009334_miropeat_default_run_repeats.bed' %{"directory":dir}
#
# #Find the genomic locations fall into miropeats repeat regions
RepeatList = []
with open(Repeat_inputFile) as repeatFile:
	for line in repeatFile:
		repstart = int(line.strip().split("\t")[1])
		repend = int(line.strip().split("\t")[2])
		for i in range(repstart,repend):
			RepeatList.append(i)
RepeatList = set(RepeatList)
print len(set(RepeatList))
#
#
outfile1 = open(sys.argv[1]+"_Mismatch_positions_with_type1.bed","w")
outfile2 = open(sys.argv[1]+"_Mismatch_positions_with_type2.bed","w")

# outfile2 = open(sys.argv[1]+sys.argv[2]+"_Mismatch_error_rates.txt","w")
#
tmpoutfile = open(sys.argv[2]+".tmp.file.aln","w")
alignment = AlignIO.read(open(sys.argv[2]), "fasta")
for record in alignment:
	#print record.id, record.seq
	tmpoutfile.write(str(record.id)+"\t"+str(record.seq)+"\n")
tmpoutfile.close()
#
#
UncoveredPos_1=[]
UncoveredPos_2=[]
with open(sys.argv[2]+".tmp.file.aln", "r") as alnfile:
	df = pd.read_csv(alnfile,sep="\t",header=None)
	dfseq = []
#     # take the names of sequences to index
	index = df[0]
	str_seq1 = df.ix[1][1]
	print "len : " +  str(len(str_seq1))
#     # create an empty list of length of sequence letters
	columns = list(range(len(str_seq1)))
#     #new data frame
	new_df = pd.DataFrame(index=index, columns= columns)
    #create new dataframe
	for i in range(len(df)):
#         # change a string of sequence to a list of sequence
		dfseq = list(df.ix[i][1])
#         # You can edit a subset of a dataframe by using loc:
#         # df.loc[<row selection>, <column selection>]
#         # place the list of sequence to the row that it belongs to in the new data frame
		new_df.loc[i:,:] = dfseq

#
#     # for each row:
	Ref_pos = 0
#	MatchCount = 0
	MissMatchCount = 0
	for i in range(len(str_seq1)):
#         # place the the sequence in a position to a set
		set_seq = set(new_df.loc[:,i])

		if new_df.loc[Ref1,i] != '-':
			Ref_pos = Ref_pos +1
		else:
			pass
# 		#put the bases at the location of the query sequence 1 and 2 into a set
		set_pair = set(new_df.loc[ [sys.argv[1],Ref1 ],i  ])

		if 'n' in set_pair:
			UncoveredPos_1.append(Ref_pos)
			print set_pair
		#if there is a mismatch error and this position is not in repeat regions, count as sequencing error.
		elif len(set_pair)>1 and Ref_pos not in RepeatList and '-' not in set_pair:
			print set_pair, Ref_pos
			MissMatchCount = MissMatchCount +1
			outfile1.write(str(Ref1)+"\t"+str(Ref_pos)+"\t"+str(Ref_pos+1)+"\n")

		######## Do the same things for Ref2, which is type 2 reference genome.
		set_pair = set(new_df.loc[ [sys.argv[1],Ref2 ],i  ])
		if 'n' in set_pair:
			UncoveredPos_2.append(Ref_pos)
			print set_pair
		#if there is a mismatch error and this position is not in repeat regions, count as sequencing error.
		elif len(set_pair)>1 and Ref_pos not in RepeatList and '-' not in set_pair:
			print set_pair, Ref_pos
			MissMatchCount = MissMatchCount +1
			outfile2.write(str(Ref1)+"\t"+str(Ref_pos)+"\t"+str(Ref_pos+1)+"\n")
#

outfile1.close()
outfile2.close()
print "Now smoothing"
outfile1 = open(sys.argv[1]+"_Mismatch_positions_with_type1_smooth.bed","w")
outfile2 = open(sys.argv[1]+"_Mismatch_positions_with_type2_smooth.bed","w")

#smoothing:
#Define a window, typically 100 bases.
win=1000
#Define a smoothing interval.
sm=200
genomeLen=Ref_pos
#genomeLen=171725

#Read the list of genomic positions that are mismatch between the aligned two genomes. Relative to Reference.
MismatchPositions =[]
with open(sys.argv[1]+"_Mismatch_positions_with_type1.bed","r") as snpfile:
    for line in snpfile:
        MismatchPositions.append(int(line.strip().split("\t")[1]))

outfile1.write("Chrom"+"\t"+"start"+"\t"+"end"+"\t"+"Sim2Type1"+"\n")
for i in range(0,genomeLen-win+sm,sm):
    windowPoss = range(i, i+win)
    #print windowPoss
    dissim=len(set(windowPoss)&set(MismatchPositions))
    CoveredLen=win-len(set(windowPoss)&set(UncoveredPos_1))+1
    PercentSim=100-(100*dissim/CoveredLen)
    #print i, dissim
#    outfile1.write("NC_007605"+"\t"+str(i)+"\t"+str(i+sm)+"\t"+str(100-dissim/10)+"\n")
    outfile1.write("NC_007605"+"\t"+str(i)+"\t"+str(i+sm)+"\t"+str(PercentSim)+"\n")

###### ----- ######
#Do the same things for Ref2, which is type 2 reference genome.

MismatchPositions =[]
with open(sys.argv[1]+"_Mismatch_positions_with_type2.bed","r") as snpfile:
    for line in snpfile:
        MismatchPositions.append(int(line.strip().split("\t")[1]))

outfile2.write("Chrom"+"\t"+"start"+"\t"+"end"+"\t"+"Sim2Type2"+"\n")
for i in range(0,genomeLen-win+sm,sm):
    windowPoss = range(i, i+win)
    #print windowPoss
    dissim=len(set(windowPoss)&set(MismatchPositions))
    CoveredLen=win-len(set(windowPoss)&set(UncoveredPos_1))+1
    PercentSim=100-(100*dissim/CoveredLen)
    #print i, dissim
    outfile2.write("NC_007605"+"\t"+str(i)+"\t"+str(i+sm)+"\t"+str(PercentSim)+"\n")


outfile1.close()
outfile2.close()
