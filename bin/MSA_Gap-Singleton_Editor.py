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


tmpoutfile = open(sys.argv[1]+".tmp.file.aln","w")
alignment = AlignIO.read(open(sys.argv[1]), "fasta")
for record in alignment:
	#print record.id, record.seq
	tmpoutfile.write(str(record.id)+"\t"+str(record.seq)+"\n")
tmpoutfile.close()


with open(sys.argv[1]+".tmp.file.aln", "r") as alnfile:
    df = pd.read_csv(alnfile,sep="\t",header=None)
    dfseq = []
    # take the names of sequences to index
    index = df[0]
    str_seq1 = df.ix[1][1]
    print "len : " +  str(len(str_seq1))
    # create an empty list of length of sequence letters
    columns = list(range(len(str_seq1)))

    editedpos = pd.DataFrame(index=index, columns= columns)

edited = alignment[:,1:2]
apbc=[]
for i in range(alignment.get_alignment_length()):

    nuccompos = collections.Counter(alignment[:,i])
    del nuccompos['n']
    del nuccompos['-']
    #print "Original Composition:", alignment[:,i]
    x=1
    for record in alignment[:,i:i+1]:

        if len(nuccompos.most_common(4))> 1 and nuccompos.most_common(4)[1][1] == 1 and record.seq == nuccompos.most_common(4)[1][0]:
            x=0
        else:
            pass
    if x == 1:
        edited = edited + alignment[:,i:i+1]
    else:
        pass

AlignIO.write(edited, "my_edited_"+sys.argv[1], "fasta")
