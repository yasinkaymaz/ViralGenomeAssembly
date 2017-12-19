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

if len(sys.argv) < 5:
	print "Welcome to new era of MSA."
	print "Please provide required arguments in proper order:"
	print "MSA_parser_cleaner.py Sample1 sample2 AlignmentFile TypeOfEBV"
	print "Feed with an alignment file in fasta format."
	print "The alignment should also include two reference EBV genomes, NC_007605 and NC_009334."
	sys.exit(1)
Ref=''
#Make sure that repeat files are in the same directory
if sys.argv[4] == '1':
	Ref='NC_007605'
	Repeat_inputFile = '%(directory)s/../resources/Annotation/Type1/NC_007605_miropeat_default_run_repeats.bed' %{"directory":dir}
else:
	Ref='NC_009334'
	Repeat_inputFile = '%(directory)s/../resources/Annotation/Type2/NC_009334_miropeat_default_run_repeats.bed' %{"directory":dir}

#Find the genomic locations fall into miropeats repeat regions
RepeatList = []
with open(Repeat_inputFile) as repeatFile:
	for line in repeatFile:
		repstart = int(line.strip().split("\t")[1])
		repend = int(line.strip().split("\t")[2])
		for i in range(repstart,repend):
			RepeatList.append(i)
RepeatList = set(RepeatList)
print len(set(RepeatList))

outfile1 = open(sys.argv[1]+sys.argv[2]+"_Mismatch_error_positions.bed","w")
outfile2 = open(sys.argv[1]+sys.argv[2]+"_Mismatch_error_rates.txt","w")

tmpoutfile = open(sys.argv[3]+".tmp.file.aln","w")
alignment = AlignIO.read(open(sys.argv[3]), "fasta")
for record in alignment:
	#print record.id, record.seq
	tmpoutfile.write(str(record.id)+"\t"+str(record.seq)+"\n")
tmpoutfile.close()

#with open(sys.argv[3], "r") as alnfile:
with open(sys.argv[3]+".tmp.file.aln", "r") as alnfile:
	df = pd.read_csv(alnfile,sep="\t",header=None)
	dfseq = []
    # take the names of sequences to index
	index = df[0]
	str_seq1 = df.ix[1][1]
	print "len : " +  str(len(str_seq1))
    # create an empty list of length of sequence letters
	columns = list(range(len(str_seq1)))
    #new data frame
	new_df = pd.DataFrame(index=index, columns= columns)
    #print new_df
    #create new dataframe
	for i in range(len(df)):
        # change a string of sequence to a list of sequence
		dfseq = list(df.ix[i][1])
        # You can edit a subset of a dataframe by using loc:
        # df.loc[<row selection>, <column selection>]
        # place the list of sequence to the row that it belongs to in the new data frame
		new_df.loc[i:,:] = dfseq

    # for each row:
	Ref_pos = 0
	MatchCount = 0
	MissMatchCount = 0
	for i in range(len(str_seq1)):
        # place the the sequence in a position to a set
		set_seq = set(new_df.loc[:,i])

		if new_df.loc[Ref,i] != '-':
			Ref_pos = Ref_pos +1
		else:
			pass
		#put the bases at the location of the query sequence 1 and 2 into a set
		set_pair = set(new_df.loc[ [sys.argv[1],sys.argv[2] ],i  ])

		#skip indels and Ns in the alignment
		if 'n' in set_pair or '-' in set_pair:
			pass
		#if there is a mismatch error and this position is not in repeat regions, count as sequencing error.

		elif len(set_pair)>1 and Ref_pos not in RepeatList:
			print set_pair, Ref_pos
			MissMatchCount = MissMatchCount +1
			outfile1.write(str(Ref)+"\t"+str(Ref_pos)+"\t"+str(Ref_pos+1)+"\n")

		#If there is a perfect match and the position is not in the repeat regions, count as match.
		elif len(set_pair) == 1 and Ref_pos not in RepeatList:
			MatchCount = MatchCount +1

outfile2.write(sys.argv[1]+"\t"+sys.argv[2]+"\t"+str(MissMatchCount)+"\t"+str(MatchCount)+"\t"+str(float(MissMatchCount)/(MissMatchCount+MatchCount)))
print "Match Count:", MatchCount
outfile1.close()
outfile2.close()
