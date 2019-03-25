#!/usr/bin/python
__requires__='pysam==0.8.1'
import pkg_resources

###standard modules
import argparse #parse initial arguements from command line
import re
import os, time
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
print "pysam version is: ", pysam.__version__
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
#pybedtools.load_path_config( 'bedtoolspath.txt')#{

###global variable to allow for documentation in location for subcommands
###this is populated as code is traversed initially until main is launched at end

##This work will all be in as start zero and end 1 coordinates UCSC unless noted
## e.g.	in ACTG ,	 CT would be designated as 1,3 (this is the same as python)

##built in defaults to incorporated key programs specific to install##
## these never change on an install unlike the config file that can require
## modificaitons between builds, organisms, etc.

#note to programmers: pysam access is all base [0,1) for tabix and fasta files
# Fastafile.fetch(self, reference=None, start=None, end=None, region=None)
#	(reference = None, start = None, end = None, region = None)

###############################################################################
# INSTALL SPECIFIC ENVIRONMENTAL PARAMETERS	##################################
###############################################################################

#binpath='/data/bailey/CNV_MIP_R01/pipeline_design/bin'
#path_clustalw='/usr/bin/clustalw'
#pybt.set_tempdir('/scratch')
#pybt.load_path_config( {	'bedtools' , '/share/pkg/bedtools/2.22.0/bin/',	 'tabix' , '/share/pkg/tabix/0.2.6/',	 'r', '/share/pkg/R/2.15.0/bin/R',	'tabix', '/share/pkg/tabix/0.2.6/'})

###############################################################################
# INITIALIZE GLOBAL VARIABLES #################################################
###############################################################################

subcommands={} #collects data at the point
#https://www.docker.io/learn_more/


###############################################################################
# LIST OF NEEDED DATA SETS ####################################################
###############################################################################

# Common SNPs(137) - SNPs with >= 1% minor allele frequency (MAF), mapping only once to reference assembly.
# Flagged SNPs(137) - SNPs < 1% minor allele frequency (MAF) (or unknown), mapping only once to reference assembly, flagged in dbSnp as "clinically associated" -- not necessarily a risk allele!
# Mult. SNPs(137) - SNPs mapping in more than one place on reference assembly.

	#All SNPs(137) - all SNPs from dbSNP mapping to reference assembly.
#ALL SNPS (138)
#gunzip -c snp138.txt.gz	| cut -f 2- | bgzip >	snp138.bed0.gz
#tabix -s 1 -b 2 -e 3 -0 snp138.bed0.gz
###Get all text UCSC file ###
#gunzip -c simpleRepeat.txt.gz	| cut -f 2- | bgzip >	simpleRepeat.bed0.gz
#tabix -s 1 -b 2 -e 3 -0 simpleRepeat.bed0.gz


####All SNPs is the most useful since it contains everything else and can be parsed to our specifications


###############################################################################
#	MAIN				###############################################################
###############################################################################

aparser=argparse.ArgumentParser()


def main( args ):
	"""Main allows selection of the main subcommand (aka function).
	Each subcommand launches a separate function. The pydoc subcommand
	launches pydoc on this overall program file.
	:param args: the main command line arguments passed minus subcommand
	"""
	#print globals().keys()

	if len(args)	== 0 or args[0] in ["h", "help", "-h", "--h", "--help","-help"] :
		verbosity= 'shortDesc'
		if args[0] in ["help" , "--help", "-help"]:
			verbosity = 'longDesc'
		program_name=os.path.basename(__file__)
		print "USAGE:",program_name, "[-h] subcommand [suboptions]"
		print "DESCRIPTION: A collection of tools to create and manipulate MIPs"
		print "SUBCOMMANDS:"
		#tw=TextWrap()
		for k in subcommands.keys():
			text=subcommands[k][verbosity]
			text= textwrap.dedent(text)
			if text:
				text =	"%s:	 %s " %(k, text )
				print textwrap.fill(text,77, initial_indent='', subsequent_indent='				 ')
		print "HELP:"
		print "pydoc			detailed documentation of program structure"
		print "-h/-help	 short / long	subcommand descriptions"
		print "For specific options:",program_name,"[subcommand] --help"
	elif args[0] == 'pydoc':
		os.system( "pydoc " + os.path.abspath(__file__) )
	elif args[0] in subcommands.keys():
		#execute sub command function
		globals()[args[0]](args[1:])
	else:
		print "unknown subcommand (" + args[0] + ") use -h for list of subcommands!"
		sys.exit(-1)
	sys.exit(0) #normal exit

#------------------------------------------------------------------------------
###############################################################################
####	SUBCOMMAND DESIGN_REGION ################################################
###############################################################################
#------------------------------------------------------------------------------
shortDescText="Create a concensus genome sequence using reads aligned to Ref genome"
longDescText="""Create a concensus genome sequence using reads aligned to Ref genome"""
subcommands['Reads2ConcensusGenome'] = { 'shortDesc':shortDescText, 'longDesc': longDescText }


def Reads2ConcensusGenome(args):
	"""This subcommand scans a bam file against a region or genome and merges read sequences resolving discrepancies based on multiple different parameters including distance from end of contig, number of differences from reference, depth of reads, etc.
	"""
	aparser.add_argument("-n", "--nameChr", required=True, help='indexed contig bam file against reference (added by Cliff and Yasin)')
	aparser.add_argument("-r", "--read_bam", required=True, help='indexed read bam file against reference')
	aparser.add_argument("-f", "--fasta_reference", required=True, help='indexed read bam file against reference')
	aparser.add_argument( "--regions_to_mask", required=True, help='repetitive or other regions to exclude')
	aparser.add_argument("--verbose", required=False, action="store_true",help='provide additional output for debugging')
	aparser.add_argument("-o","--outfilebase", help='output base path and filename for generated data')
	##fix need to change the program name##
	fullname= aparser.prog + " Reads2ConcensusGenome"
	aparser.prog=fullname
	args=aparser.parse_args(args=args)

	###load reference sequence ###
	from Bio import SeqIO
	refseqs = list( SeqIO.parse(args.fasta_reference, "fasta"))
	print refseqs[0].id
	print refseqs[0]
	newref=[None] * (len(refseqs[0].seq)+1)

	###load mask regions ###

#	tmp=FlexBed.readfile(args.regions_to_mask, args.nameChr)
#	maskbed=FlexBed(tmp)
	if args.read_bam:
		readbam=pysam.Samfile(args.read_bam, 'rb')
		print "Mapped # of reads:", readbam.mapped
		print "Reference of reads mapped to is ", readbam.references
		print "Reference length :", readbam.lengths

		NucleotidesfromReads = {}
		for readcol in readbam.pileup(readbam.references[0]):
			readcolnucleotides = []
			deletedReads=[]
#			print ("\ncoverage at base %s = %s" %(readcol.pos, readcol.n))

			for pileupread in readcol.pileups:
	#		if not pileupread.is_del and not pileupread.is_refskip:
            # query position is None if is_del or is_refskip is set.
				if pileupread.indel > 0:
#					print pileupread.indel, readcol.reference_pos, pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+pileupread.indel+1]
					readcolnucleotides.append(pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+pileupread.indel+1])
				elif pileupread.is_del:
					print readcol.reference_pos
					print pileupread.query_position
					deletedReads.append(readcol.reference_pos)
					readcolnucleotides.append(pileupread.alignment.query_sequence[pileupread.query_position])
				else:
					readcolnucleotides.append(pileupread.alignment.query_sequence[pileupread.query_position])
				readnucleotidecomposition=collections.Counter(readcolnucleotides)
		#		print readnucleotidecomposition
				ConsensusReadBase=''
				if len(readnucleotidecomposition.most_common(2)) > 1:
					if readnucleotidecomposition.most_common(2)[0][1] == readnucleotidecomposition.most_common(2)[1][1]:
#						print readcol.pos, readnucleotidecomposition.most_common(1)[0], readnucleotidecomposition.most_common(2)[1], refseqs[0].seq[readcol.reference_pos]
						ConsensusReadBase=refseqs[0].seq[readcol.reference_pos]
					else:
						ConsensusReadBase=readnucleotidecomposition.most_common(1)[0][0]
				else:
					ConsensusReadBase=readnucleotidecomposition.most_common(1)[0][0]
#					print readcol.pos, readnucleotidecomposition.most_common(1)[0]
			if len(deletedReads) > readcol.n/2:
#				print len(deletedReads), readcol.reference_pos
				ConsensusReadBase=''
			else:
				pass
#		print readnucleotidecomposition
#		print(readnucleotidecomposition.most_common(1))
#		print readcol.reference_pos, readnucleotidecomposition.most_common(1)[0][0]
			NucleotidesfromReads[readcol.reference_pos]=ConsensusReadBase

#	print NucleotidesfromReads
#				print ('\tbase in read %s = %s' %(pileupread.alignment.query_name,pileupread.alignment.query_sequence[pileupread.query_position]))

	start = time.time()
	consseqfile = open ( args.outfilebase + "_consensusSequence.fa", 'w')
	ConsSeq=''
	NucleotidesfromReads_keyset=set(NucleotidesfromReads.keys())

	for i in range(0,len(refseqs[0].seq)+1):
#	for i in range(0,12000):
		if i in NucleotidesfromReads_keyset:
#		if i in NucleotidesfromReads.keys():
			ConsSeq+=str(NucleotidesfromReads[i])
#			print i, NucleotidesfromReads[i]
		else:
			ConsSeq+='N'
#			print "N"
	consseqfile.write(">ConcensusSeq"+"\n")
	consseqfile.write(ConsSeq+"\n")
	consseqfile.close()
	print round(time.time() - start,1)

def calculate_endposition(readposition, readlength):
	endpos=0 # zero is the edge , positive bases are left, negative are negative
	if readposition < readlength/2 :
		endpos=readposition + 1	#left side will be positive bases from end
	else:
		endpos=-1*(readlength-readposition)	#right side will be negative bases from end
	return endpos


class ContigSegment():
	""" A Data structure for Segments
	"""
	def __init__(self):
		self.seq=None
		self.ref_begin=None
		self.ref_end=None
		self.name=None
		self.length=None
		self.rightmost_endpos=None
		self.leftmost_endpos=None
		self.left_mismatches=0
		self.right_mismatches=0
		self.left_last_good=None
		self.right_last_good=None
		self.left_bases=[]
		self.right_bases=[]
		self.left_aligned=0
		self.right_aligned=0
		self.left_trim=None	#number of bases to trim back
		self.right_trim =None #base corrections
		self.fixes=[] # [[start,end,newbases]] newbases="" for deletion
		self.data=None
		self.base_aligned_unique=None




##----------------------------------------------------------------------------
class FlexBed:
	"""flexible format bed class for in memory manipulations
	note: beds are just lists but this provides static methods
	to operate on such lists.
	the only fixed data requirements are:
	#col(1)seq/chr name
	#col(2)begin position (zero coordinate)
	#col(3)end position (one cooridnate)
	#col(4)name
	note: a bed is only defined seq:start:end:orient	(no more required)
	seq
	"""
	def __init__(self,bed):
		"""sets up a standard bed format file"""
		self.bed=bed
		self.current_positions=[0,0,0,0,0]
		self.previous_sequences=[ [] ] * len(self.current_positions)
		#0=seq,#1=start,#2=stop,#3=name, #4=score,#5=orientation
	def bed_set (self, bed):
		""" set the bed with in the bed object , provides no checking """
		self.bed_clear()
		self.bed=bed
	def bed_clear(self):
		self.bed=[]
		self.current_positions= [position] * len(self.current_positions)
		self.previous_sequences=[ [] ] * len(self.current_positions)

	def walker_reset(self,counter=None, position=0):
		""" walker reset will reset the current_positions to index position zero
		by default """
		if counter!=None:
			self.current_positions[counter]=position
			self.previous_sequences[counter]=[]
		else:
			self.current_positions= [position] * len(self.current_positions)
			self.previous_sequences=[ [] ] * len(self.current_positions)
	def walker_get_counters(self):
		return list(self.current_positions)
	def bed_sort(self):
		#TODO: need to actually add a sort function
		if len(self.bed)==0:
			return
		from operator import itemgetter
		#print self.bed[0:4]
		#check to be sure that first element
		self.bed.sort(key=itemgetter(0,1,2))
		return

	def writefile(self,filename):
		outfile=open(filename, 'w')
		#print outfile, len(self.bed)
		for i	in self.bed:
			outline=outfile.write("\t".join([str(x) for x in i])+ "\n")
		outfile.close()
		return
	def walker_step_upto(self, seq, start, counter=0):
	#	for i in xrange (self.current_positions[counter],len(self.bed)):
	#		feature=self.bed[i]
	#		if seq!=feature[0]:
				return


	def walker_get_range(self, seq, start, end, counter=0,trim=False, \
											 autoincrement=False):
		"""walks along chromosome in zero one coordinates
		this requires that comparisons are done in order
		sorted sequences and shorted coordinates"""
		newbed = []
		for i in xrange ( self.current_positions[counter], len(self.bed) ):
			feature=self.bed[i]
			#print feature
			if seq==feature[0] :	#bed feature on same seq as our range
					if feature[2] < start: #bed feature less than start of range
						if autoincrement==True:
							self.current_positions[counter]=i+1 #so increment counter
						continue #and go to next
					newbegin=max(start,feature[1])
					newend=min(end,feature[2])
					if	newend-newbegin	> 0: #bed feature within our range
						newbed.append(feature)
						continue
					else:	 # query feature << bed feature
						break # already past it so stop
			if seq < feature[0]:	 #stop and bail	current seq position is greater than seq
				break
		return newbed

		# return list (empty if nothing intersecting)
		return []

	def calc_bp_overlap2interval(self,seq,begin,end ):
		totalbases=0
		for bedsegment in self.bed:
			intersect=Interval.intersect2seq( seq,begin,end, bedsegment[0], bedsegment[1],bedsegment[2] )
			#print intersect[3]
			totalbases+=intersect[3]
		return totalbases
	@staticmethod
	def readfile(filepath, check='chr', checknum=True, seq=0, start=1, stop=2,
							 name=3,score=4,orient=5):
		"""read a tab-delimited bed file and return bed list of lists"""

		lines = open(filepath, "r").readlines()
		print lines
		bed= FlexBed.readlist(lines, check=check,checknum=checknum,seq=seq,start=start,
											stop=stop,name=name,score=score,orient=orient)
		return bed

	@staticmethod
	def readlist(textlist, check='chr',checknum=True, seq=0,
							 start=1, stop=2,	name=3,score=4,orient=5):
		bed=[]
		#print "textlist",textlist
		usedcolumns=[seq,start,stop,name,score,orient]
		for line in textlist:
			if line.startswith("#"):
				continue
			x=line.strip().split("\t")
			y=[ x[seq],x[start],x[stop],None,None,None ]
			if name!=None:
				y[3]=x[name]
			if score!=None:
				y[4]=x[score]
			if orient!=None:
				y[5]=x[orient]
			for i in xrange(0,len(x)):
				if i in usedcolumns:
					continue
				else:
					y.append(x[i])
			if checknum:
				y[1]=int(y[1])
				y[2]=int(y[2])
			if not x[0].startswith(check):
				print >> sys.stderr, x[self.chrom]+ " is not a proper seq name with : ("+check+")"
				sys.exit()
			bed.append(y)
		return bed

	#write a bed
		#bed.split("\t")
		#print bed
	#		bed.check(check='chr',checknum=True)
	#def bed.readlist(self,tab_delimited_list):
		#"""parse list of tab-delimited text into bed list of lists"""
	#tablist.strip()
		#return
	#def bed.check(tab_delmited_bed_list):

###THIS is now used in my code ###
class	Interval:
	""" ucsc genome coordinates [a,b) """
	@staticmethod
	def intersect2interval(begin1,end1, begin2,end2):
		"intersect [a,b)	[c,d) returns new begin and end of intersect"
		newbegin=max(begin1,begin2)
		newend=min(end1,end2)
		intersect=newend-newbegin
		if intersect <=0:
			return (0,'','')
		return intersect, newbegin, newend
	@staticmethod
	def intervalWellWithinInterval2(begin1,end1,begin2,end2,flank):
		"must be ensconced by a certain distance "
		#TODO: finish this
	@staticmethod
	def intervalWithinInterval2(begin1,end1,begin2,end2):
		"is [begin1,end1) within [begin2,end2)"
		if begin1 >= begin2 and end1 <= end2:
			return True
		return False

	@staticmethod
	def intersect2seq( seq1,begin1,end1, seq2,begin2,end2):
		if seq1!=seq2:
			print "MISS", seq1 , seq2
			return [None,None,None,0]
		newbegin=max(begin1,begin2)
		newend=min(end1,end2)
		intersect=newend-newbegin
		if intersect<=0:
			return [None,None,None,0]
		return (seq1,newbegin,newend,intersect)



###WHAT IS THIS CLASS DOING IN MY CODE!!!!!
class coordmapper:
	"""given a mapping file this class allows cross mapping.
	Mapping that isn't unique will cause failure."""

	def __init__(self):
		"""list of coordinates"""
		self.coord=[]
		self.reverseNames={}	#pushes list positions
		self.forwardNames={}	#pushes list positions
	def load_coord(self,coordlist):
		"""seq, begin, end, name, seq, begin, end, orient,
		First seq alway forward/plus orientation.
		Requires 0,1 coordinates (aka bed file)
		No data checks currently!
		"""
		self.coord=coordlist
		return
	def append(self, coordarray):
		"""[seq1,begin1,end1, seq2,begin2,end2,orient, name1,name2,data...]
		No data checks currently!
		"""
		self.coord.append(coordarray)
		return
	def recalculateSeqLookup(self):
		"""this recalculates dictionaries for more rapid lookups"""
	def checkData(self):
		"""NOT IMPLEMENTED:	could check consistency of the data"""
	def forwardMap(self, seq,begin,end=None):
		"""given seq1,begin1,[end1], return seq2 equilvalent"""
		return
	def reverseMap(self,seq,begin,end=None):
		"""given se2,begin,[end2], lookup seq1 equilvaent..."""
		return


###############################################################################
###############################################################################
#------------------------------------------------------------------------------
#main is left to the bottom so that global variables can be populated
if __name__ == "__main__":
	if len (sys.argv)==1:
		sys.argv.append("-h")	#if no command then it is a cry for help
	main(sys.argv[1:])
