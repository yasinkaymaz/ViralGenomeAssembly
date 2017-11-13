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


###############################################################################
####	SUBCOMMAND: MERGE BAM ALIGNED CONTIGS AGAINST REFERENCE #################
###############################################################################
shortDescText="merge ref aligned contigs that overlap when placed on reference"
longDescText="""merge contigs that have been aligned to a reference--
finds overlaps and deals with poorly trimmed ends (containing extraneous seq)
as well as with poor depth and resolves sequence differences due to these
artifcats
"""
subcommands['ref_merge_contigs'] = { 'shortDesc':shortDescText, 'longDesc': longDescText }


def ref_merge_contigs(args):
	"""This subcommand scans a bam file against a region or genome and merges sequences (contigs) resolving discrepancies based on multiple different parameters including distance from end of contig, number of differences from reference, depth of reads, etc.
	"""

	aparser.add_argument("-c", "--contig_bam", required=True, help='indexed contig bam file against reference')
	aparser.add_argument("-n", "--nameChr", required=True, help='indexed contig bam file against reference (added by Cliff and Yasin)')
	aparser.add_argument("-r", "--read_bam", required=False, help='indexed read bam file against reference')
	aparser.add_argument("-f", "--fasta_reference", required=True, help='indexed read bam file against reference')
	aparser.add_argument("-u", "--unmapped_trim_length" ,required=True, type=int, help='number of bases to remove when end not aligned')
	aparser.add_argument("-t", "--mapped_trim_max_length", required=True, type=int, help='number of base pairs from end of contig to potentially trim')

	aparser.add_argument("-g", "--good_end_max_mismatch", required=True, type=int, help='number of mismatches to tolerate before calling bad end')
	aparser.add_argument("--min_unique_bp", required=True,type=int ,help = ' number of unique bases')
	aparser.add_argument( "--regions_to_mask", required=False, help='repetitive or other regions to exclude')
	aparser.add_argument("--wipe", required=False, action="store_true",help='wipe any intermediate tmp files (fully_reruns process)')
	aparser.add_argument("--verbose", required=False, action="store_true",help='provide additional output for debugging')
	aparser.add_argument("--check",required=False, action="store_true",help='initial check of genome databases for integrity')
	aparser.add_argument("-o","--outfilebase", help='output base path and filename for generated data')
	##fix need to change the program name##
	fullname= aparser.prog + " ref_merge_contigs"
	aparser.prog=fullname
	args=aparser.parse_args(args=args)
	contigbam=pysam.Samfile(args.contig_bam, 'rb')

	###load reference sequence ###
	from Bio import SeqIO
	refseqs = list( SeqIO.parse(args.fasta_reference, "fasta"))
	print refseqs[0].id
	print refseqs[0]
	newref=[None] * (len(refseqs[0].seq)+1)

	###load mask regions ###
	start = time.time()
	tmp=FlexBed.readfile(args.regions_to_mask, args.nameChr)
	maskbed=FlexBed(tmp)
	if args.read_bam:
		readbam=pysam.Samfile(args.read_bam, 'rb')
		print "Mapped # of reads:", readbam.mapped
		print "Reference of reads mapped to is ", readbam.references
		print "Reference length :", readbam.lengths
	#	Read_pile = readbam.pileup(readbam.references[0])
		NucleotidesfromReads = {}
		for readcol in readbam.pileup(readbam.references[0]):
			readcolnucleotides = []
			deletedReads=[]
#			print ("\ncoverage at base %s = %s" %(readcol.pos, readcol.n))
			readdepth=0

			starttime = time.time()
			for pileupread in readcol.pileups:
	#		if not pileupread.is_del and not pileupread.is_refskip:
            # query position is None if is_del or is_refskip is set.

				if pileupread.indel > 0:
#					print pileupread.indel, readcol.reference_pos, pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+pileupread.indel+1]
					readcolnucleotides.append(pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+pileupread.indel+1])
				elif pileupread.is_del:
			#		print readcol.reference_pos
					deletedReads.append(readcol.reference_pos)
					readcolnucleotides.append(pileupread.alignment.query_sequence[pileupread.query_position])
				else:
					readcolnucleotides.append(pileupread.alignment.query_sequence[pileupread.query_position])
				readnucleotidecomposition=collections.Counter(readcolnucleotides)

#				print readnucleotidecomposition,
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
				readdepth=readdepth+1
				#Use a threshold for depth of coverage for sampling nucleotide composition.
				if readdepth == 1000:
					break
			if len(deletedReads) > readcol.n/2:
#				print len(deletedReads), readcol.reference_pos
				ConsensusReadBase=''
			else:
				pass
			print ("coverage at base %s = %s" %(readcol.pos, readcol.n)), ConsensusReadBase, round(time.time() - starttime, 1)
#		print readnucleotidecomposition
#		print(readnucleotidecomposition.most_common(1))
#		print readcol.reference_pos, readnucleotidecomposition.most_common(1)[0][0]
			NucleotidesfromReads[readcol.reference_pos]=ConsensusReadBase

#	print NucleotidesfromReads

	NucleotidesfromReads_keyset=set(NucleotidesfromReads.keys())
	print round(time.time() - start,1)
	### output basic statistics for the bam file
 # print	"NO COORDINATE #:", contigbam.nocoordinate
	print "Contigs MAPPED #:", contigbam.mapped
 # print	"REFERENCE #: contigbam.nreferences
	print "REFERENCES:", contigbam.references
	print "REF LENGTHS:", contigbam.lengths


	### NEW REFERENCE GENOME ###
 # refbases= [None]* (len refseqs[0].Seq +10000)

	### GENERATE LIST OF CONTIGS ###
	print "READING CONTIGS..."
	contigs=[]
	contigdict={}
	for s in contigbam:
		a=ContigSegment()
#		print "a:",a
		a.data=s
		a.name = a.data.query_name
		a.length = a.data.query_length
		contigdict[a.name]=a
		contigs.append(a)
		if a.data.is_unmapped:
			###add code here to dump unmapped contigs
			continue
		else:
			a.start = a.data.reference_start
			a.end = a.data.reference_end


	### SORT CONTIGS DESCENDING ORDER ###
	print "SORTING CONTIGS ON THE BASIS OF LENGTH LARGEST TO SMALLEST..."
	contigs = sorted(contigs, key=lambda cntg:cntg.data.query_length, reverse=True)
	print "Sorting done."
	#############################################################################
	### RUN THROUGH THE ENTIRE PILEUP ###########################################
	pile = contigbam.pileup(contigbam.references[0])
	for col in pile:
		ref_base=''
		contig_bases=[]
		for read in col.pileups:
			readbase = read.alignment.seq[read.query_position]
			ref_base = refseqs[0].seq[col.reference_pos]
#			print read.query_position, readbase, ref_base
#			print read.query_position
			### the entire inserted sequence ###
			readseq = read.alignment.query_sequence[read.query_position: read.query_position+read.indel+1]
			qname= read.alignment.query_name
			### DETERMINE BASES FROM EITHER END (NEGATIVE -- RIGHT)(POSITIVE --LEFT)
			endpos=calculate_endposition(read.query_position, read.alignment.query_length)
			if contigdict[qname].leftmost_endpos==None:
				contigdict[qname].leftmost_endpos=endpos
			contigdict[qname].rightmost_endpos=endpos
			#print col.pos, ref_base, readbase, readseq, read.indel
			if	abs (endpos) <= int(args.mapped_trim_max_length):
				d= [endpos,read.indel,readbase, ref_base, read.query_position ,readseq ]
				#print qname, endpos, read.query_position, read.alignment.query_length, "	 ",
				if endpos <0:
					contigdict[qname].right_bases.append(d)
					contigdict[qname].right_aligned+=1
					if readbase!=ref_base:
						contigdict[qname].right_mismatches+=1
				else:
					contigdict[qname].left_bases.append(d)
					contigdict[qname].left_aligned+=1
					if readbase!=ref_base:
						contigdict[qname].left_mismatches+=1

			 # print qname, contigdict[qname].left_mismatches,":" , contigdict[qname].left_aligned, "	",contigdict[qname].right_mismatches,						":",	contigdict[qname].right_aligned
	#############################################################################
	print "CONTIG PROCESSING FOR DUP/REPEAT CONTENT AND FOR END ALIGNMENT QUALITY..."
	for c in contigs:
		### find repetitive content
		#print c.leftmost_endpos, c.rightmost_endpos, c.length
		if c.data.is_unmapped==False:
			qname = c.name
			repeatoverlap=maskbed.calc_bp_overlap2interval(contigbam.references[0],c.start, c.end)
			c.uniquebp = c.end-c.start - repeatoverlap
			#print c.start, c.end, "LEN", c.length, "UNIQUE",c.uniquebp
		#check right alignment
		for x in c.right_bases:
			if x[1]==0 and x[2]==x[3]:
				c.right_last_good=x[0]
			else:
				break
		#check left alignment
		c.left_bases=sorted(c.left_bases , key= lambda x:x[0],reverse=True)
		for x in c.left_bases:
			if x[1]==0 and x[2]==x[3]:
				c.left_last_good=x[0]
			else:
				break
		#print c.left_bases
		#print c.left_last_good
	#############################################################################
	pile = contigbam.pileup(contigbam.references[0]) #TODO: learn out to reset
	print "PICK THE BEST BASE AND RECORD CHANGES..."
	start = time.time()
	for col in pile:
		ref_base=''
		contig_bases=[]
		for read in col.pileups:
			readbase = read.alignment.query_sequence[read.query_position]
			ref_base = refseqs[0].seq[col.reference_pos]
			readseq = read.alignment.query_sequence[read.query_position: read.query_position+read.indel+1]
		#	print col.pos, "R",ref_base, readbase, readseq

			### DETERMINE BASES FROM END ###
			endpos=calculate_endposition(read.query_position, read.alignment.query_length)
			valid=True
			repeat=False
			if abs (endpos) <= int(args.mapped_trim_max_length):
				last_good_base = contigdict[read.alignment.query_name].left_last_good
				if endpos < 0:
					last_good_base = contigdict[read.alignment.query_name].right_last_good
				if last_good_base == None or abs(endpos) < abs(last_good_base):
					valid=False
				if contigdict[read.alignment.query_name].uniquebp < args.min_unique_bp:
					repeat=True
			readdata = [readbase, readseq, endpos, read.indel, read.query_position, read.alignment.query_length, read.alignment.query_name, valid, repeat]
			contig_bases.append(readdata)

		### out of contig loop #################
		allbasecounts = [ c[1] for c in contig_bases ]
		allbasecounts = Counter(allbasecounts)
		#VALID CONTIGS ARE ANCHORED IN NONREPETIVE SEQUENCE
		contig_bases_valid=[c for c in contig_bases if c[7]==True and c[8]==False]
		validbasecounts=[c[1] for c in contig_bases_valid ]
		validbasecounts=Counter(validbasecounts)
		validbases=validbasecounts.keys()
#		print "here is valid bases", validbases, col.reference_pos
		newbase="NOT DETERMINED"

		if len (validbases) == 1:
			### we have single sequence in all our contigs that are valid ###
#			if ref_base == validbases[0] and col.reference_pos in NucleotidesfromReads_keyset: #contigs same as reference AND READ PILE UP IS PRESENT
#				newbase=NucleotidesfromReads[col.reference_pos]
#			elif ref_base != validbases[0] and col.reference_pos in NucleotidesfromReads_keyset:
			#	newbase= validbases[0].lower() #lowercase (for reference)
			if col.reference_pos in NucleotidesfromReads_keyset:
				newbase=NucleotidesfromReads[col.reference_pos]#*************************************************
			else: #contigs represent a substitution
				newbase=validbases[0] #uppercase
		elif len (validbases) > 1:
			if len(validbases[0]) > 1 or len(validbases[1]) > 1:
				print "This should be inserted: ", validbases[0], validbases[1]

			if ref_base in validbases:
				#MULTIPLE WELL PLACED VALID CONTIGS#
				maxlen=max( [c[5] for c in contig_bases_valid ] )
				lendict={}
				for c in contig_bases_valid:
					lendict[c[5]] = c[1]
				#WARNINGS START
				if lendict[maxlen] != ref_base and args.verbose:
					print col.pos, ref_base
					print "WARNING ref_base", ref_base, " is not the largest contig", c[5], lendict
					print "	", contig_bases
					for c in contig_bases_valid:
						if len(c[1])>10:
							print col.pos, ref_base, c[1],"WARNINGSKIPPING LARGE INSERTION"
				#WARNINGS END
				if col.reference_pos in NucleotidesfromReads_keyset:
					newbase=NucleotidesfromReads[col.reference_pos]#*************************************************
					#HERE WHERE YOU NEED TO PICK MAJOR VARIANT!!!
	#!#				print "Here you should pick Major var:", NucleotidesfromReads[col.reference_pos], col.reference_pos
				else:
					newbase=ref_base.lower()#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			else:
				if args.verbose:
					print "WARNING: NO REF BASE DESPITE MULTIPLE BASES--PICK LARGER CONTIG"
				contig_bases=sorted(contig_bases, key=lambda cntg:cntg[5], reverse=True)
				if col.reference_pos in NucleotidesfromReads_keyset:
					newbase=NucleotidesfromReads[col.reference_pos]#*************************************************
				else:
					newbase=contig_bases[0][1]
				#HERE WHERE YOU CAN PICK MAJOR VARIANT!!!
		elif len (validbases) ==0:
			check=[]
			for c in contig_bases:
				span=abs(c[2])+c[3]
				if c[2] > 0:	#left with
					#print "LEFT TRIM ",span, args.unmapped_trim_length
					if span > int(args.unmapped_trim_length):
						contigdict[c[6]].left_trim=args.unmapped_trim_length
						starttrim=contigdict[c[6]].left_trim-c[2]+1
						newbase=c[1][starttrim:]
						###add a removal step
						#print "fixiing trim", starttrim, contigdict[c[6]].left_trim, newbase
					else:
						 newbase=''
				else:	 #C[2] <0	 left
					if len(c[1]) > 5:
						print contig_bases
						print "	",allbasecounts
						print "RIGHT ",span, args.unmapped_trim_length
						print
					else: #toss it
						newbase=''
				check.append(newbase)
			check=sorted(check,key=lambda x:len(x),reverse=True)
			if len(check) > 1:
				print "here is check: ", check
			if check[0]=='':
				check[0]=None
			newbase=None
		### OUT OF THE LOOP
		newref[col.pos]=newbase
		if newbase == "NOT DETERMINED" :
			print "NEWBASE",newbase
			print ref_base, col.pos
			print contig_bases
			print allbasecounts
			print contig_bases_valid
			sys.exit(0)
	print round(time.time() - start,1)
	#######################################################################3#####
	### DUMP THE OUTPUT ########################################################
###	HERE WE CAN START LOOKING AT GENOMIC REGIONS WITHOUT CONTIGS TO SEE IF THERE IS READ SUPPORT FROM READ ALIGNMENT FILE,
###	from i=0 to i=endOfGenome: IF i is not col.reference_pos; then look at alignment file to see if there is any read aligned to i. in other words, check if col.reference_pos in NucleotidesfromReads_keyset:
###	if if col.reference_pos in NucleotidesfromReads_keyset: TRUE then take the base as newbase
#	print newref
	lastpos=None
	contig=0
	gfile = open ( args.outfilebase + "_genome.fa", 'w')
	gfile.write( '>'+args.outfilebase + "\n")
	cfile = open ( args.outfilebase + "_contig_new.fa", 'w')
	contigseq=''
	for seq in newref:
		if seq==None:
			gfile.write('N')
		else:
			gfile.write(seq)
			contigseq+=seq
		if seq==None and lastpos!=None:
			name=str(contig).zfill(3)
			cfile.write( ">"+name+"\n")
			cfile.write(contigseq+"\n")
			contigseq=''
			contig+=1
		lastpos=seq
	gfile.close()
	cfile.close()

	###### DUMP THE OLD CONTIGS ALBEIT TRIMMED ##################################
	cfile = open ( args.outfilebase + "_contig_old.fa", 'w')
	for c in contigs:
		header=">"+c.name
		seq=''
		begin=None
		end=None
		if c.data.is_unmapped:
			begin=args.unmapped_trim_length
			end=c.length-args.unmapped_trim_length
			header+="_tt"

		else:
			#left
####			print c.left_last_good, c.right_last_good, c.name
			begin=args.unmapped_trim_length
			left="t"
			if c.left_last_good !=None and begin <= args.mapped_trim_max_length:
				begin=c.left_last_good-1
				left="T"
			#right
			end=c.length-args.unmapped_trim_length
			right="t"
			if c.right_last_good != None and c.right_last_good <= args.mapped_trim_max_length:
				end= c.length + c.right_last_good + 1
				right="T"
			header+="_"+left+right


		header+=str(begin) + "-"+str(end) + " OLDLEN" + str(c.length)
####		print len (c.data.query_sequence), c.length
		seq= c.data.query_sequence[begin:end]
		cfile.write(header+ "\n"+ seq + "\n")
	cfile.close()



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
