#!/usr/bin/python

###standard modules
import argparse #parse initial arguements from command line
import re		
import os
import os.path
import pickle
import sys
import subprocess
import textwrap
import glob
import csv
from collections import defaultdict, Counter
#from collections import namedtuple
import collections
#key science modules
import pysam #for reading sam files as well as other files
import numpy as np
import pybedtools as pybt
#import numpy.string_
#import csv 
## key bio python modules 

#print pybt.check_for_bedtools()
import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import AlignIO
from Bio import Phylo #not used yet but good for reading dendograms from clustalw

###other biologic modules
import pybedtools as pybt



#with open('%s' %GenomeFasta_name, "r") as GenomeFastaFile:
ref=0
aln=0
fileout=open("Type1_to_type2_reflocmatch.txt","w")
#with open("/home/yk42w/codes/EBVseq/NC_009334_aligned2type1.fa", "r") as GenomeFastaFile:
with open("/home/yk42w/codes/EBVseq/NC_007605_aligned2type2.fa", "r") as GenomeFastaFile:
        FastaLines = GenomeFastaFile.readlines()
	for line in FastaLines:
		for i in range(len(line)):
			if line[i] != "-":
				ref=ref+1
				aln=aln+1

#				fileout.write("NC_009334:"+str(ref)+"\t"+"NC_009334:"+str(aln)+"\n")
				fileout.write("NC_007605:"+str(ref)+"\t"+"NC_007605:"+str(aln)+"\n")
			else:
				aln=aln+1

#				fileout.write("NC_009334:"+str(ref)+"\t"+"NC_009334:"+str(aln)+"\n")
				fileout.write("NC_007605:"+str(ref)+"\t"+"NC_007605:"+str(aln)+"\n")




