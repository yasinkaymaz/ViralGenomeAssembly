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

subcommands={} #collects data at the point
#https://www.docker.io/learn_more/

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
	#aparser.add_argument("-n", "--nameChr", required=True, help='indexed contig bam file against reference (added by Cliff and Yasin)')
	aparser.add_argument("-r", "--read_bam", required=True, help='indexed read bam file against reference')
	aparser.add_argument("-f", "--fasta_reference", required=True, help='indexed read bam file against reference')
#	aparser.add_argument( "--regions_to_mask", required=True, help='repetitive or other regions to exclude')
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
	#tmp=FlexBed.readfile(args.regions_to_mask, args.nameChr)
	#maskbed=FlexBed(tmp)
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
			#here check to see if coverage is above a threshold:
			if readcol.n < 5:
					NucleotidesfromReads[readcol.reference_pos]='N'
			else:
					#else do this:
					readdepth=0
					starttime = time.time()

					for pileupread in readcol.pileups:
			#		if not pileupread.is_del and not pileupread.is_refskip:
		            # query position is None if is_del or is_refskip is set.
						readdepth=readdepth+1
						if readdepth == 1000:
							break
						if pileupread.indel > 0:
		#					print pileupread.indel, readcol.reference_pos, pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+pileupread.indel+1]
							readcolnucleotides.append(pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+pileupread.indel+1])
						elif pileupread.is_del:
					#		print readcol.reference_pos
					#		print pileupread.query_position
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
		#			print ("coverage at base %s = %s" %(readcol.pos, readcol.n)), ConsensusReadBase, round(time.time() - starttime, 1)

					NucleotidesfromReads[readcol.reference_pos]=ConsensusReadBase

#	print NucleotidesfromReads
#				print ('\tbase in read %s = %s' %(pileupread.alignment.query_name,pileupread.alignment.query_sequence[pileupread.query_position]))
	#print round(time.time() - start,1)

	start = time.time()
	consseqfile = open ( args.outfilebase + "_genome.fa", 'w')
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
	consseqfile.write('>'+args.outfilebase +"_Fixed"+ "\n")
	consseqfile.write(ConsSeq+"\n")
	consseqfile.close()
	print round(time.time() - start,1)


shortDescText="Create a concensus genome sequence using reads aligned to Ref genome"
longDescText="""Create a concensus genome sequence using reads aligned to Ref genome"""
subcommands['PushSecondMajorGenome'] = { 'shortDesc':shortDescText, 'longDesc': longDescText }

def PushSecondMajorGenome(args):
	"""This subcommand scans a bam file against a region or genome and merges read sequences resolving discrepancies based on multiple different parameters including distance from end of contig, number of differences from reference, depth of reads, etc.
	"""
	#aparser.add_argument("-n", "--nameChr", required=True, help='indexed contig bam file against reference (added by Cliff and Yasin)')
	aparser.add_argument("-r", "--read_bam", required=True, help='indexed read bam file against reference')
	aparser.add_argument("-f", "--fasta_reference", required=True, help='indexed read bam file against reference')
#	aparser.add_argument( "--regions_to_mask", required=True, help='repetitive or other regions to exclude')
	aparser.add_argument("--verbose", required=False, action="store_true",help='provide additional output for debugging')
	aparser.add_argument("-o","--outfilebase", help='output base path and filename for generated data')
	##fix need to change the program name##
	fullname= aparser.prog + " PushSecondMajorGenome"
	aparser.prog=fullname
	args=aparser.parse_args(args=args)

	###load reference sequence ###
	from Bio import SeqIO
	refseqs = list( SeqIO.parse(args.fasta_reference, "fasta"))
	print refseqs[0].id
	print refseqs[0]
	newref=[None] * (len(refseqs[0].seq)+1)

	###load mask regions ###
	#tmp=FlexBed.readfile(args.regions_to_mask, args.nameChr)
	#maskbed=FlexBed(tmp)
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
			#here check to see if coverage is above a threshold:
			if readcol.n < 5:
					NucleotidesfromReads[readcol.reference_pos]='N'
			else:
					#else do this:
					readdepth=0
					starttime = time.time()

					for pileupread in readcol.pileups:
			#		if not pileupread.is_del and not pileupread.is_refskip:
		            # query position is None if is_del or is_refskip is set.
						readdepth=readdepth+1
						if readdepth == 1000:
							break
						if pileupread.indel > 0:
		#					print pileupread.indel, readcol.reference_pos, pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+pileupread.indel+1]
							readcolnucleotides.append(pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+pileupread.indel+1])
						elif pileupread.is_del:
					#		print readcol.reference_pos
					#		print pileupread.query_position
							deletedReads.append(readcol.reference_pos)
							readcolnucleotides.append(pileupread.alignment.query_sequence[pileupread.query_position])
						else:
							readcolnucleotides.append(pileupread.alignment.query_sequence[pileupread.query_position])
					readnucleotidecomposition=collections.Counter(readcolnucleotides)
		#			print readnucleotidecomposition.most_common(2), readcol.n

					ConsensusReadBase=''
					if len(readnucleotidecomposition.most_common(2)) > 1:
						minorPercent= 100*readnucleotidecomposition.most_common(2)[1][1]/readcol.n
						if minorPercent > 25.0 and readnucleotidecomposition.most_common(2)[1][1] > 10:

							print readnucleotidecomposition.most_common(2), minorPercent, readcol.n, readnucleotidecomposition.most_common(2)[1][0]
							ConsensusReadBase=readnucleotidecomposition.most_common(2)[1][0]

						elif readnucleotidecomposition.most_common(2)[0][1] == readnucleotidecomposition.most_common(2)[1][1]:

				#			print readnucleotidecomposition.most_common(2), readcol.n, readcol.pos, readnucleotidecomposition.most_common(1)[0], readnucleotidecomposition.most_common(2)[1], refseqs[0].seq[readcol.reference_pos]
							ConsensusReadBase=refseqs[0].seq[readcol.reference_pos]
						else:
							ConsensusReadBase=readnucleotidecomposition.most_common(1)[0][0]
					else:
						ConsensusReadBase=readnucleotidecomposition.most_common(1)[0][0]
		#				print readcol.pos, readnucleotidecomposition.most_common(1)[0]

					if len(deletedReads) > readcol.n/2:
		#			print len(deletedReads), readcol.reference_pos
						ConsensusReadBase=''
					else:
						pass
		#	print readnucleotidecomposition
		#	print(readnucleotidecomposition.most_common(1))
		#	print readcol.reference_pos, readnucleotidecomposition.most_common(1)[0][0]
		#			print ("coverage at base %s = %s" %(readcol.pos, readcol.n)), ConsensusReadBase, round(time.time() - starttime, 1)

					NucleotidesfromReads[readcol.reference_pos]=ConsensusReadBase

#	print NucleotidesfromReads
#				print ('\tbase in read %s = %s' %(pileupread.alignment.query_name,pileupread.alignment.query_sequence[pileupread.query_position]))
	#print round(time.time() - start,1)

	start = time.time()
	consseqfile = open ( args.outfilebase + "_genome.fa", 'w')
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
	consseqfile.write('>'+args.outfilebase +"_Fixed"+ "\n")
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


###############################################################################
###############################################################################
#------------------------------------------------------------------------------
#main is left to the bottom so that global variables can be populated
if __name__ == "__main__":
	if len (sys.argv)==1:
		sys.argv.append("-h")	#if no command then it is a cry for help
	main(sys.argv[1:])
