#!/usr/bin/python
"""
This code is a collaborative effort with Sevtap Duman Northeastern University, systems security lab.
"""
import os
import pandas as pd
import sys
import collections
from Bio import AlignIO
import Bio.Align
from Bio import SeqIO
dir = os.path.dirname(__file__)
print os.getcwd()
print dir

alndata = []
#sys.argv[1] -> Alignment File
#sys.argv[2] -> type of virus
"""
The purpose of this script is to patch sequences partially missing in de novo assembly genomes.
"""

if len(sys.argv) < 3:
	print "Please provide required arguments in proper order:"
	print "MSA2contigs_bed.py AlignmentFile ReferenceType[1/2] "
	sys.exit(1)

Ref=''
#Make sure that repeat files are in the same directory
if sys.argv[2] == '1':
	Ref='NC_007605'
	Repeat_inputFile = '%(directory)s/../resources/Annotation/Type1/NC_007605_miropeat_default_run_repeats.bed' %{"directory":dir}

else:
	Ref='NC_009334'
	Repeat_inputFile = '%(directory)s/../resources/Annotation/Type2/NC_009334_miropeat_default_run_repeats.bed' %{"directory":dir}


compare = lambda x, y: collections.Counter(x) == collections.Counter(y)
#gap_ignore = gi
gi=100
alignment = AlignIO.read(open(sys.argv[1]), "fasta")
for record in alignment :
	outfile= open(sys.argv[1]+'.'+record.id+'.contigs.bed', 'w')
	coords = []
	for i in range(len(record.seq)):
		#if record.seq[i] != 'N' and record.seq[i] != '-' and i != len(record.seq)-1:
		if record.seq[i] != 'N' and record.seq[i] != '-' and i != len(record.seq)-1:
			#print record.seq[i:i+10], record.seq[i:i+10].count('N')
			coords.append(i)
		elif len(coords) != 0 and (record.seq[i:i+gi].count('N') == gi or record.seq[i:i+gi].count('-') == gi):
			print Ref, min(coords), max(coords)
			outfile.write(Ref+'\t'+str(min(coords))+'\t'+str(max(coords))+'\n' )
			coords=[]
		else:
			pass
	outfile.close()
	# 	#
		# if record.seq[i] != 'N':
		# 	end=i
		# 	print record.id, str(start), str(end), record.seq[i]
		#
		# else:
		# 	start=i
		# 	print record.id, str(start), str(end), record.seq[i]
		#


#
#
#
# PatchBedfile = sys.argv[4]
#
# patchrecords = list(SeqIO.parse(sys.argv[3], "fasta"))
#
#
#
# print("Alignment length %i" % alignment.get_alignment_length())
#
# #bad alignment positions = bap
# bap = []
#
# edited = alignment[:,1:2]
# print edited
#
#
# Ref_pos = 0
# for record in alignment:
#     if record.id == Ref:
#         for i in range(len(record.seq)):
#             if record.seq[i] == '-':
#                 Ref_pos = Ref_pos
#                 bap.append(i)
#             else:
# 				Ref_pos = Ref_pos + 1
#
# print "BAP:", len(bap)
#
# #Here parse patch file. put the coordinates in a dict.
#
#
#
# def PullLoci(BedFile):
#
#     PatchDict={}
#     with open(BedFile) as bfile:
#         r=0
#         for line in bfile:
#             PositionList = []
#             BaseDict = {}
#             sample = str(line.strip().split("\t")[0])
#             start = int(line.strip().split("\t")[1])
#             end = int(line.strip().split("\t")[2])
#             b=0
#             for i in range(start,end+1):
#                 PositionList.append(i)
#                 BaseDict[i] = patchrecords[r].seq[b]
#                 b=b+1
#
#             PatchDict[sample] = BaseDict
#             r=r+1
#
# 	return PatchDict
#
# PatchLoci = PullLoci(PatchBedfile)
# #print PatchLoci
#
# for i in range(alignment.get_alignment_length()):
#     apbc = []
#     for record in alignment:
#         #print record.seq[i], record.id
#         x = record.seq[i]



	#alignment Position Base Composition = apbc

#	for record in alignment:
	#	apbc.append(record.seq[i])
    #check sampleid to see if in patch seq file:
    #need patchseq start position on ref
    #patchseq end pos is str+len(patchseq)

	#if compare(set(apbc), set(['-','n']) ) or compare(set(apbc), set(['-']) ) or compare(set(apbc), set(['n']) ) or (i in set(bap)):
	#	print record.seq[i], record.id, set(apbc)
	#	pass
	#else:
	#	edited = edited + alignment[:,i:i+1]


#edited = edited[:,1:]
#print edited

#AlignIO.write(edited, "my_ICed_"+sys.argv[1], "fasta")

#print len(bap)




#outfile1.close()
#outfile2.close()
