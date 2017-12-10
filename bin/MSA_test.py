#!/usr/bin/python
"""
This code is a collaborative effort with Sevtap Duman Northeastern University, systems security lab.
"""
import os
import pandas  as pd
import sys
from Bio import AlignIO
from Bio import SeqIO
import Bio.Align
dir = os.path.dirname(__file__)
print os.getcwd()
print dir
from collections import defaultdict, Counter
#from collections import namedtuple
import collections
import numpy
import phylopandas as pd

alndata = []

#sys.argv[1] -> Alignment File

if len(sys.argv) < 2:
	print "Welcome to new era of MSA."
	print "Please provide required arguments in proper order:"
	print "MSA_Gap-Singleton_Editor.py AlignmentFile"
	print "Feed with an alignment file in fasta format."
	sys.exit(1)


alignment = pd.read_fasta(sys.argv[1])
#print alignment
print alignment
print alignment['sequence'][1][1]
print alignment['sequence'][3][1]

NumofSeqs = len(alignment['sequence'])
alignmentSeqLen = len(alignment['sequence'][1])
#for i in range(alignmentSeqLen):
for i in range(11):
    alignedBases = []
    for x in range(NumofSeqs):
        a=1
        alignedBases.append(alignment['sequence'][x][i])
    nuccompos = collections.Counter(alignedBases)
    del nuccompos['n']
    del nuccompos['-']
    print alignedBases
    print nuccompos

#edited = alignment[:,1:2]
#
# #for i in range(alignment.get_alignment_length()):
# for i in range(10):
#
#     nuccompos = collections.Counter(alignment[:,i])
#     del nuccompos['n']
#     del nuccompos['-']
#
#     for record in alignment[:,i:i+1]:
# #	if compare(set(apbc), set(['-','n']) ) or compare(set(apbc), set(['-']) ) or compare(set(apbc), set(['n']) ):
#
#     #    print i, nuccompos.most_common(4), nuccompos.most_common(4)[1][0], alignment[:,i:i+1]
#         x=0
#         print "Original Composition:", alignment[:,i]
#         if len(nuccompos.most_common(4))> 1 and nuccompos.most_common(4)[1][1] == 1 and record.seq == nuccompos.most_common(4)[1][0]:
#     #        print("%s %s %i" % (record.seq, record.id, len(record)))
#     #        if record.seq == nuccompos.most_common(4)[1][0]:
#                 a=1
#                 editedpos.loc[x:,i:i] = 'n'
#     #            edited[x,i] = 'n'
#     #            print "singleton is:", record.seq, "and id is:", record.id, alignment[x:,i:i]
#         else:
#                 editedpos.loc[x:,i:i] = record.seq
#     #            edited[x,i] = record.seq
#                 x +=1
#         #print edited
#
#         #edited[:,i:i+1] = editedpos.loc[:,i:i]
#         #edited = edited + editedpos.loc[:,i:i]
# # edited = edited[:,1:]
# print editedpos.loc[:,1:9]
#

#alignment.to_fasta('new_alignment.fasta', sequence_col='alignment')
