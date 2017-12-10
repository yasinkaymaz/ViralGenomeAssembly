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
dir = os.path.dirname(__file__)
print os.getcwd()
print dir

alndata = []
#sys.argv[1] -> Alignment File
#sys.argv[2] -> type of virus

if len(sys.argv) < 3:
	print "Please provide required arguments in proper order:"
	print "MSA_InsertCleaner.py AlignmentFile TypeOfEBV"
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

alignment = AlignIO.read(open(sys.argv[1]), "fasta")
print("Alignment length %i" % alignment.get_alignment_length())
for record in alignment :
#    print(record.seq + " " + record.id)
    print(record.id), record.seq.count('-'), record.seq.count('n')

#bad alignment positions = bap
bap = []

edited = alignment[:,1:2]
print edited


Ref_pos = 0
for record in alignment:
    if record.id == Ref:
        for i in range(len(record.seq)):
            if record.seq[i] == '-':
                Ref_pos = Ref_pos
                bap.append(i)
            else:
				Ref_pos = Ref_pos + 1

print "BAP:", len(bap)

for i in range(alignment.get_alignment_length()):

	#alignment Position Base Composition = apbc
	apbc = []

	for record in alignment:
		apbc.append(record.seq[i])

	if compare(set(apbc), set(['-','n']) ) or compare(set(apbc), set(['-']) ) or compare(set(apbc), set(['n']) ) or (i in set(bap)):
		print record.seq[i], record.id, set(apbc)
		pass
	else:
		edited = edited + alignment[:,i:i+1]


edited = edited[:,1:]
print edited

AlignIO.write(edited, "my_ICed_"+sys.argv[1], "fasta")

print len(bap)




#outfile1.close()
#outfile2.close()
