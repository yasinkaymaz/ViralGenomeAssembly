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
## e.g.  in ACTG ,   CT would be designated as 1,3 (this is the same as python)

##built in defaults to incorporated key programs specific to install##
## these never change on an install unlike the config file that can require
## modificaitons between builds, organisms, etc.  

#note to programmers: pysam access is all base [0,1) for tabix and fasta files
# Fastafile.fetch(self, reference=None, start=None, end=None, region=None)
#  (reference = None, start = None, end = None, region = None)

###############################################################################
# INSTALL SPECIFIC ENVIRONMENTAL PARAMETERS  ##################################
###############################################################################

binpath='/data/bailey/CNV_MIP_R01/pipeline_design/bin'
path_clustalw='/usr/bin/clustalw'
#pybt.set_tempdir('/scratch')
#pybt.load_path_config( {  'bedtools' , '/share/pkg/bedtools/2.22.0/bin/',   'tabix' , '/share/pkg/tabix/0.2.6/',   'r', '/share/pkg/R/2.15.0/bin/R',  'tabix', '/share/pkg/tabix/0.2.6/'})

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
#gunzip -c snp138.txt.gz  | cut -f 2- | bgzip >  snp138.bed0.gz
#tabix -s 1 -b 2 -e 3 -0 snp138.bed0.gz 
###Get all text UCSC file ###
#gunzip -c simpleRepeat.txt.gz  | cut -f 2- | bgzip >  simpleRepeat.bed0.gz
#tabix -s 1 -b 2 -e 3 -0 simpleRepeat.bed0.gz


####All SNPs is the most useful since it contains everything else and can be parsed to our specifications


###############################################################################
#  MAIN        ###############################################################
###############################################################################

aparser=argparse.ArgumentParser()


def main( args ):
  """Main allows selection of the main subcommand (aka function).
  Each subcommand launches a separate function. The pydoc subcommand 
  launches pydoc on this overall program file. 
  :param args: the main command line arguments passed minus subcommand
  """
  #print globals().keys()
   
  if len(args)  == 0 or args[0] in ["h", "help", "-h", "--h", "--help","-help"] :
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
        text =  "%s:   %s " %(k, text )
        print textwrap.fill(text,77, initial_indent='', subsequent_indent='         ')
    print "HELP:"
    print "pydoc      detailed documentation of program structure"
    print "-h/-help   short / long  subcommand descriptions"
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
####  SUBCOMMAND DESIGN_REGION ################################################
###############################################################################
#------------------------------------------------------------------------------


###############################################################################
####  SUBCOMMAND: MERGE BAM ALIGNED CONTIGS AGAINST REFERENCE #################
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

  aparser.add_argument("-c",  "--contig_bam", required=True, help='indexed contig bam file against reference')
  aparser.add_argument("-n",  "--nameChr", required=True, help='indexed contig bam file against reference (added by Cliff and Yasin)')
  aparser.add_argument("-r",  "--read_bam", required=False, help='indexed read bam file against reference')
  aparser.add_argument("-f",  "--fasta_reference", required=True, help='indexed read bam file against reference')
  aparser.add_argument("-u", "--unmapped_trim_length" ,required=True, type=int, help='number of bases to remove when end not aligned')
  aparser.add_argument("-t", "--mapped_trim_max_length", required=True,  type=int ,  help='number of base pairs from end of contig to potentially trim')

  aparser.add_argument("-g", "--good_end_max_mismatch", required=True, type=int ,    help='number of mismatches to tolerate before calling bad end') 
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
  newref=[None]  * (len(refseqs[0].seq)+1)
 
  ###load mask regions ###

  tmp=FlexBed.readfile(args.regions_to_mask, args.nameChr)
  maskbed=FlexBed(tmp)
  if args.read_bam:
    readbam=pysam.Samfile(args.read_bam, 'rb')
  ### output basic statistics for the bam file
 # print  "NO COORDINATE #:", contigbam.nocoordinate
  print  "MAPPED #:", contigbam.mapped
 # print  "REFERENCE #: contigbam.nreferences
  print  "REFERENCES:", contigbam.references
  print  "REF LENGTHS:", contigbam.lengths


  ### NEW REFERENCE GENOME ###  
 # refbases= [None]* (len refseqs[0].Seq +10000)  
  
  ### GENERATE LIST OF CONTIGS ###
  print  "READING CONTIGS..."
  contigs=[]
  contigdict={}
  for s in contigbam:
    a=ContigSegment()
    print "a:",a
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
  print  "SORTING CONTIGS ON THE BASIS OF LENGTH LARGEST TO SMALLEST..."
  contigs =  sorted(contigs, key=lambda cntg:cntg.data.query_length, reverse=True)
  
  #############################################################################
  ### RUN THROUGH THE ENTIRE PILEUP ###########################################
  pile = contigbam.pileup(contigbam.references[0])
  for col in pile:
    ref_base=''
    contig_bases=[]
    for read in col.pileups:
      readbase = read.alignment.seq[read.query_position]
      ref_base =  refseqs[0].seq[col.reference_pos]
      ### the entire inserted sequence ###
      readseq =   read.alignment.query_sequence[read.query_position: read.query_position+read.indel+1]  
      qname= read.alignment.query_name
      ### DETERMINE BASES FROM EITHER END (NEGATIVE -- RIGHT)(POSITIVE --LEFT)
      endpos=calculate_endposition(read.query_position, read.alignment.query_length)
      if contigdict[qname].leftmost_endpos==None:
        contigdict[qname].leftmost_endpos=endpos
      contigdict[qname].rightmost_endpos=endpos
      #print col.pos, ref_base, readbase, readseq, read.indel
      if  abs (endpos) <= int(args.mapped_trim_max_length):
        d= [endpos,read.indel,readbase, ref_base, read.query_position ,readseq ]
        #print qname, endpos, read.query_position, read.alignment.query_length, "   ",
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

       # print qname, contigdict[qname].left_mismatches,":" , contigdict[qname].left_aligned, "  ",contigdict[qname].right_mismatches,            ":",  contigdict[qname].right_aligned
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
    for  x in c.right_bases:
      if x[1]==0 and x[2]==x[3]:
        c.right_last_good=x[0]
      else:
        break
    #check left alignment
    c.left_bases=sorted(c.left_bases , key= lambda x:x[0],reverse=True)
    for  x in c.left_bases:
      if x[1]==0 and x[2]==x[3]:
        c.left_last_good=x[0]
      else:
        break
    #print c.left_bases
    #print c.left_last_good
  #############################################################################
  pile = contigbam.pileup(contigbam.references[0]) #TODO: learn out to reset
  print "PICK THE BEST BASE AND RECORD CHANGES..."
  for col in pile:
    ref_base=''
    contig_bases=[]
    for read in col.pileups:
      readbase = read.alignment.query_sequence[read.query_position]
      ref_base =  refseqs[0].seq[col.reference_pos]
      readseq =   read.alignment.query_sequence[read.query_position: read.query_position+read.indel+1]  
      #print col.pos, "R",ref_base, readbase, readseq

      ### DETERMINE BASES FROM END ###
      endpos=calculate_endposition(read.query_position, read.alignment.query_length)
      valid=True
      repeat=False
      if  abs (endpos) <= int(args.mapped_trim_max_length):
        last_good_base = contigdict[read.alignment.query_name].left_last_good
        if endpos < 0:
          last_good_base = contigdict[read.alignment.query_name].right_last_good
        if last_good_base == None or  abs(endpos) <  abs(last_good_base):
          valid=False
        if contigdict[read.alignment.query_name].uniquebp < args.min_unique_bp:
          repeat=True
      readdata = [readbase, readseq, endpos, read.indel  , read.query_position, read.alignment.query_length, read.alignment.query_name, valid, repeat]
      contig_bases.append(readdata)
      
    ### out of contig loop #################
    allbasecounts = [ c[1] for c in contig_bases ]
    allbasecounts = Counter(allbasecounts)
    #VALID CONTIGS ARE ANCHORED IN NONREPETIVE SEQUENCE
    contig_bases_valid=[c for c in contig_bases if c[7]==True and c[8]==False]
    validbasecounts=[c[1] for c in contig_bases_valid ]
    validbasecounts=Counter(validbasecounts)
    validbases=validbasecounts.keys()
    newbase="NOT DETERMINED" 
     
    if len (validbases) == 1:
      ### we have single sequence in all our contigs that are valid ###
      if ref_base == validbases[0]:       #contigs same as reference
        newbase= validbases[0].lower()           #lowercase (for reference)
      else:                               #contigs represent a substitution
        newbase=validbases[0]                    #uppercase
    elif len (validbases) > 1:
      if  ref_base in validbases:
        #MULTIPLE WELL PLACED VALID CONTIGS#
        maxlen=max( [c[5] for c in contig_bases_valid ] )
        lendict={}
        for c in contig_bases_valid:
          lendict[c[5]] = c[1]
        #WARNINGS START
        if lendict[maxlen] != ref_base and args.verbose:
          print col.pos, ref_base
          print "WARNING ref_base", ref_base, " is not the largest contig", c[5], lendict
          print "  ",contig_bases 
          for c in contig_bases_valid:
            if len(c[1])>10: 
              print col.pos, ref_base, c[1],"WARNINGSKIPPING LARGE INSERTION"
        #WARNINGS END
        newbase=ref_base.lower()
      else:
        if args.verbose:
          print "WARNING: NO REF BASE DESPITE MULTIPLE BASES--PICK LARGER CONTIG"
        contig_bases=sorted(contig_bases, key=lambda cntg:cntg[5], reverse=True)
        newbase=contig_bases[0][1]
        
    elif len (validbases) ==0:
      check=[]
      for c in contig_bases: 
        span=abs(c[2])+c[3]
        if c[2] > 0:  #left with
          #print "LEFT TRIM ",span, args.unmapped_trim_length
          if span > int(args.unmapped_trim_length):
            contigdict[c[6]].left_trim=args.unmapped_trim_length
            starttrim=contigdict[c[6]].left_trim-c[2]+1
            newbase=c[1][starttrim:]
            ###add a removal step
            #print "fixiing trim", starttrim, contigdict[c[6]].left_trim, newbase
          else:
             newbase=''
        else:   #C[2] <0   left
          if len(c[1]) > 5:
            print contig_bases
            print "  ",allbasecounts
            print "RIGHT ",span, args.unmapped_trim_length
            print 
          else: #toss it
            newbase=''
        check.append(newbase)
      check=sorted(check,key=lambda x:len(x),reverse=True)
      if len(check) > 1:
        print check
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
  
  #######################################################################3#####
  ### DUMP THE OUTPUT ########################################################
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
      print c.left_last_good, c.right_last_good, c.name
      begin=args.unmapped_trim_length
      left="t"
      if  c.left_last_good !=None and  begin <= args.mapped_trim_max_length:
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
    print len (c.data.query_sequence), c.length
    seq= c.data.query_sequence[begin:end]
    cfile.write(header+ "\n"+ seq + "\n")
  cfile.close()
    
    

  
    


def calculate_endposition(readposition, readlength):
  endpos=0  # zero is the edge , positive bases are left, negative are negative
  if   readposition < readlength/2 :
    endpos=readposition + 1  #left side will be positive bases from end
  else:
    endpos=-1*(readlength-readposition)  #right side will be negative bases from end 
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
    self.left_trim=None  #number of bases to trim back
    self.right_trim =None #base corrections
    self.fixes=[] # [[start,end,newbases]] newbases="" for deletion 
    self.data=None
    self.base_aligned_unique=None
    



###############################################################################
# SUBCOMMAND: CHECK KEY MODULES ###############################################
###############################################################################

shortDescText="run functional checks on key required external modules"
longDescText=""" simple checks on the following modules: pysam
It also contains commented documentation in terms of implementing it. 
"""
subcommands['check_modules'] = { 'shortDesc':shortDescText, 'longDesc': longDescText }

def check_modules(args):
  print "*Checking pybedtools..."
  a = pybt.example_bedtool('a.bed')
  b = pybt.example_bedtool('b.bed')
  print a.intersect(b)

###############################################################################
# SUBCOMMAND: DESIGN MIPS IN A GIVEN REGION ###################################
###############################################################################

shortDescText="creates paralogous MIPs from given paralogous regions and other targets"
longDescText=""" takes input of paralogous region coordinates, aligns and designs
pMIPs and incorporates other targets based on the given input while excluding troublesome 
regions within such as tandem repeats, yount transposable lements"""
subcommands['design_region'] = { 'shortDesc':shortDescText, 'longDesc': longDescText }

def design_region(args):
  """This subcommand function runs overarching code to design mips 
  in a given region (including SNPs, paralogs, indes, deletions, as well
  as specific include regions and exclude regions.  It is a large 
  overaching routine drawing on many subroutines and modules to
  provide the needed flexibility to rapidly design accurate MIPs.
  """
  #parse the subcommand
  #TODO: potentially could move aparser to a separate routine
  aparser.add_argument("-r",  "--region_info", required=True, help='info of region targets and paralogs')
  aparser.add_argument("-c",  "--config_file", required=False, help='a config file that can load fixed genome/tabix options')
  aparser.add_argument("--genome_fasta_bgzip", required=False,   help='indexed fasta bgzipped file of genome (all of genome including allelic')
  aparser.add_argument("--genome_highdepth_bam_dir", required=False, help='directory containing high depth genomes for removing SNPs')
  
  #aparser.add_argument("--genome_bowtie2db_nonallelic", required=False, help='bowtie (bt) formateted genome database lacking allelic seqs')
  aparser.add_argument("--genome_rm_tabix", required=False, help='repeatmasker tab-delimited files indexed by tabix')
  aparser.add_argument("--genome_trf_tabix", required=False, help='tandem repeat finder (trf) tab-delimited files indexed by tabix')
  aparser.add_argument("--genome_snps_tabix", required=False, help='ucsc snp file (all of dbSNP) indexed by tabix')
  aparser.add_argument("--action_include", required=False, help='default action if not given in region info file (span:snp:cds)' )
  aparser.add_argument("--novelseq_fasta_bgzip", required=False, help='additonal novel sequence (e.g. extra paralogs) need to design paralogous MIPs')
  #extra options#
  aparser.add_argument("-t", "--intermediate_tmp_dir", required=False,     help='dir to store intermediate files')
  aparser.add_argument("--wipe", required=False, action="store_true",help='wipe any intermediate tmp files (fully_reruns process)')
  aparser.add_argument("--check",required=False, action="store_true",help='initial check of genome databases for integrity')
  aparser.add_argument("-o","--outfilebase", help='output base path and filename for generated data')
  args=aparser.parse_args(args=args)

  ### add args from config file if present ###
  loaded_opts=None
  if (args.config_file != None):
    print "Loading configuration file..."
    loaded_opts=config_file_load(args.config_file)
  ### define args if in config file into args ###
  for a in args.__dict__:
    #print "VAL", args.__dict__[a]
    if args.__dict__[a] == None:
      if loaded_opts!=None and a in loaded_opts:
        args.__dict__[a] = loaded_opts[a]
      elif a.startswith('genome') or a.startswith('tabix'):
        #exit because genome/tabix must be set
        print >>sys.stderr, "ERROR: Argument:",a, \
            " -- not defined on command line or in config file!"
        sys.exit(1)
    if a.startswith('genome') or a.startswith('tabix'):
      val=args.__dict__[a]
      #print "testing.. .", a,val,
      if not os.path.exists(val):
        print >> sys.stderr, a,val, "ERROR: file/dir does not exist!"
        sys.exit(1)
      
  ### load args
  highdepth_genomes=[]
  if (args.genome_highdepth_bam_dir):
    highdepth_genomes=glob.glob(args.genome_highdepth_bam_dir+"/*.bam")
  
  print "Initializing key data sources through pysam..."
  pysam_fgenome=pysam.Fastafile(args.genome_fasta_bgzip) #f is for fasta and fast
  tabix_rmasker=pysam.Tabixfile(args.genome_rm_tabix)
  tabix_snps=pysam.Tabixfile(args.genome_snps_tabix)
  tabix_trf=pysam.Tabixfile(args.genome_trf_tabix)
  
  
  ############################################################################  
  ### additonal argument preprocessing, checks and setup ####
  ### and to see if data is loading 
  if args.check==True :
    
    #genome_fasta_bg_zip 
    print "Quick Checks of Data Sources (by no means comprehensive)..."
    print "\n*CHECK genome_fasta_bgzip", args.genome_fasta_bgzip , "..."
    ##pysam pysam_fgenome.fetch operates in 0,1 coordiantes
    tmpseq=pysam_fgenome.fetch(reference='chr1',start=999999,end=1000000)
    if tmpseq == "T": #base 0,1)
      print "   OK hg19  base 1,000,000 ", tmpseq, " == T" #in hg19 this should be a T
    #seq=pysam_fgenome.fetch(pysam_fgenome,'chr1',999999,1000000,'+') #base 1,1]  

    #snps 
    print "\n*CHECK genome_snps_tabix", args.genome_snps_tabix, "..."
    #print " SNPS rs334 for sickle cell chr11:5248232-5248232 in hg19..."
    for bedline in tabix_snps.fetch('chr11', 5248231,5248232):
      if "rs334" in bedline:
        print textwrap.fill("OK hg19"+bedline,75, initial_indent='  ', subsequent_indent='  ')

    #trf
    print "\n*CHECK genome_trf_tabix", args.genome_trf_tabix, "..."
    #print " repeat chr11:5247112-5247146 is ATTTATATGCAGAAAT "
    for bedline in tabix_trf.fetch('chr11', 5247144,5247145):
      if "ATTTATATGCAGAAAT" in bedline:
        print textwrap.fill("OK hg19"+bedline,75, initial_indent='  ', subsequent_indent='  ')
    print "\n*CHECK genome_rm_tabix", args.genome_rm_tabix, "..."
    for bedline in tabix_rmasker.fetch('chr11', 5247410,5247415):
      if "LINE" in bedline:
        print textwrap.fill("OK hg19 "+bedline,75, initial_indent='  ', subsequent_indent='  ')

    #high coverage
    print "\n*CHECK genome_high_depth_bam_dir\n   ",args.genome_highdepth_bam_dir
    print textwrap.fill(str(highdepth_genomes),78, initial_indent='  ', subsequent_indent='  ')  
    ###bail out and quit
    sys.exit(0)
    
  #############################################################################
  outbase=args.outfilebase
  if outbase==None:
    outbase=args.region_info
    outbase=re.sub(r'\.rinfo',"",outbase,1) + "_out"
  print "BASENAME for output is ",outbase
  
  logstream=open( (outbase+".log"),'w')
  print >> logstream , "help what is going on\n"
  #logstream.write("cat\nhelp"+'text')
  ###load info file ###
  
  print "LOADING REGION INFORMATION FILE..."
  rinfo=info_file_load(args.region_info)
  print "SINGLEVALUE PARAMETERS", rinfo['#VAL'].keys()
  print "MULTIVALUE DATA TYPES",   rinfo['COL'].keys()
  print ""
  flatdata=[]
  print "\nANALYZING THE REGION:",rinfo['NAME']
 
  ########## SEGMENT PROCESSING ###############################################
 
  segments= [x for x in rinfo["REGION"]]
  segments.sort()
  print   "INITIAL ANALYSIS OF SEGMENTS:", segments
  for sname in segments:
    seg=Segment(sname)
    flatdata.append(seg) 
    print "*EXTRACTING SEG:", seg.name
    seg.copies=[ x for x in rinfo["REGION"][seg.name] ]
    seg.copies.sort()
    print "**LOADING COPY DATA AND SEQ... ", seg.copies
    seg.clustal_infile = "tmp."+ rinfo['#VAL']['NAME'] +".seg."+ seg.name +".fa"
    fasta_handle=open(seg.clustal_infile,'w')
    seg.copydata=rinfo["REGION"][seg.name]
    for copy in seg.copies:
      #note s and p are text
      cdata=seg.copydata[copy] 
      print "**"+copy, cdata#pointer to abbreviate
      rec=indexfa_grab_seq(pysam_fgenome,cdata['seqn'],cdata['seqb'], \
          cdata['seqe'],cdata['orient'])
      rec.name= copy +"."+ cdata['copyname'] #make proper sequence record
      SeqIO.write(rec, fasta_handle , "fasta")  
      cdata['dna']=rec.seq
    fasta_handle.close()
    #cmd="fasta_length3.pl " + seg.clustal_infile
    #print cmd
    #out=run_command (cmd)
    #print out 
    
    ###### POPULATE REFERENCE SEGMENT DATA ####################################
    refpnt=rinfo["REGION"][seg.name][seg.copies[0]]
    seg.refdna=refpnt['dna']
    seg.refseqname=refpnt['seqn']
    seg.refseqbegin=refpnt['seqb']
    seg.refseqend=refpnt['seqe']
    seg.refseqcopyname=refpnt['copyname']
    
  #############################################################################
  #############################################################################
  print "\nANALYZING MIP DIRECTIONS FROM REGION INFO FILE..."
  for featuretype in ['FORCE','TARGET','EXCLUDE']:
    if not featuretype in rinfo.keys() or len(rinfo[featuretype])==0:
      print "*",featuretype, " : does not exist"
      continue
    print "*",featuretype, ": assigning to proper copy..."
    print "**(checking for orphan targets outside of a segment)"
  #do not be so bold as to presume that these are TARGETS are sorted
    for i, feat in list (enumerate (rinfo[featuretype])):
      hit=0
      hits=[]
      for seg in flatdata:
         for c in seg.copies: 
           cd= seg.copydata[c]
           intersect=Interval.intersect2seq(cd['seqn'],cd['seqb'],cd['seqe'], \
               feat['seqn'],feat['seqb'],feat['seqe'] )
           if intersect :
             #print intersect
             feat['segment']=seg.name
             feat['copy']=c
             hits=[seg.name,copy]
             hit+=1
             if not featuretype in cd:
               cd[featuretype]=[]
             cd[featuretype].append(i)
      if hit == 0:
        print "*WARNING!: no segment contains feature#",i, "===>", feat
      if hit >= 2: 
        print "*WARNING!: possible boundary issue with feature in multiple regions  feature#" \
                ,i, "===>", feat, " overlaps these copies ", hits
  
  #############################################################################
  #############################################################################
  print "\nGENERATE AND PROCESS CLUSTALW..."
  for seg in flatdata:
    print  "*"+ seg.name, "FLAT REF (", seg.refseqname, seg.refseqbegin, seg.refseqend, ")"
    if len(seg.copies) ==1:
      print "*SEGMENT ",seg.name," has only one copy -- no paralogs"
      seg.msalign=None
      seg.pdiffs=FlexBed([])
      seg.pdashes=FlexBed([])
    else:
      seg.has_paralogs=True
      print "*GENERATING CLUSTALW with >1 copy -- paralogs present ..." 
      seg.clustal_outfile = seg.clustal_infile + "sta" #.fasta  
      command=path_clustalw+" -INFILE="+seg.clustal_infile \
          +" -TYPE=DNA -OUTPUT=FASTA " + " -OUTORDER=INPUT "
      if args.wipe:
        os.remove(seg.clustal_outfile)
      if os.path.isfile(seg.clustal_outfile):
        print "**CLUSTALW skipped as ", seg.clustal_outfile, \
            "exists: use -wipe to rerun."
      else:
        print command
        out = run_command(command)
        if out[0]:
          print sys.stderr >> "ERROR CLUSTALW:",out[0]
          sys.exit(1)
      #########################################################################
      print "*LOADING CLUSTALW MULTIPLE ALIGNMENT..."
      seg.msalign=MultipleAlign()
      seg.msalign.loadAlignment(seg.clustal_outfile,"fasta")
      #if seg.name=="S2":
        #system.exit(1)
        #seg.msalign.generateInternalCoords(debug=True)
      #else:
      seg.msalign.generateInternalCoords()
      region_coords=[ [seg.copydata[c]['seqn'], seg.copydata[c]['seqb'], \
              seg.copydata[c]['seqe'], seg.copydata[c]['orient'], \
              seg.copydata[c]['copyname']   ]   for c in seg.copies ]
      #print region_coords
      #print "region", (region_coords[0][2]- region_coords[0][1])
      seg.msalign.generateExternalCoords(region_coords)
      #########################################################################
      print "*CALCULATING  DIFFERENCES AND PRIMER EXCLUDES BASED ON MULTIPLE ALIGNMENT..."
      seg.msalign.getSeqCoord(0,seg.msalign.getAlignLength()-1, coord='external')
      #########################################################################
      pdiffs= seg.msalign.getDifferences(seqout=0, types="base", external=True)
      print "**COPY/PARALOG DIFFS FOR MIP TARGETING:", len(pdiffs)
      #print pdiffs[0]
      seg.pdiffs=FlexBed(pdiffs)
     ##########################################################################
      pdashes=seg.msalign.getDifferences(seqout=0, types="dash", external=True)
      print "**EXCLUDE DASH POSITIONS FROM PRIMERS:", len(pdashes)
      #print pdashes[0]
      seg.pdashes=FlexBed(pdashes)

  #############################################################################
  #############################################################################
 #  for seg in flatdata:  
 #   print "LOADING SNPs for the region into FlexBed"
    #This should bring in all SNPs of interest#
    #may also want to pull out SNPs from other paralog copies#
    #TODO: pull out other snps and cross translate them onto ref paralog
  #  lines=[]
   # for r in tabix_snps.fetch(seg.refseqname,seg.refseqbegin,seg.refseqend):
    #  snp=r.split("\t")
     # chr_count=2
      #if snp[22] != '':
       # chr_count=sum(map(float,snp[22].rstrip(",").split(",") ))
      #if chr_count<4:
       #   continue      
      #lines.append(r)
      #print r
    #print len(lines)
    #lines=FlexBed.readlist(lines)
    #seg.allsnps=FlexBed(lines)
  
  #############################################################################  
  #############################################################################  
  print "LOADING SNPs from reference and other paralogous regions."
  for seg in flatdata:
    allsnps=[]
    includesnps=[]
    to_include=['coding-synon','missense','cds-indel', 'nonsense']
    print "*", seg.name, "processing SNPs"
   # print seg.has_paralogs
   # print seg.copies
   # print seg.copydata
    #if seg.has_paralogs:
    exclude_names=defaultdict(int)
    for cval in xrange (0, len(seg.copies)):
      #print "SEG", seg.name,"COPYVAL",cval
      cname=seg.copies[cval]
      cdata = seg.copydata[cname]
      #print cdata
      for r in tabix_snps.fetch(cdata['seqn'],cdata['seqb'],cdata['seqe']):
        snp=r.split("\t")
        chr_count=2
        if snp[22] != '':  ###counts
          chr_count=sum(map(float,snp[22].rstrip(",").split(",") ))
        if chr_count<10:
          continue
        #print snp[14],snp[22], snp[1],snp[2]
        #if not (cdata['seqb']<= int(snp[1])<=cdata['seqe']):
          #print "outside"
        if seg.has_paralogs:  
          oldname=cname+":"+snp[0]+','+ snp[1]+ ','+ snp[2]
          (c,b,e)= seg.msalign.ext2extInterval(cval, snp[0],int(snp[1]),int(snp[2]), 0 )
          snp[0]=c
          snp[1]=b
          snp[2]=e
          snp[4]=oldname
        allsnps.append(snp)
        #print "ALLSNPS", len(allsnps)
        if snp[14] in to_include:
          includesnps.append(snp)
        else:
          exclude_names[snp[14]]+=1
          
    print "**", "Include: ",len(includesnps), " All Snps: ", len(allsnps)
    print "**", "Excluded: ", exclude_names.keys()
    seg.includesnps=FlexBed(includesnps)
    seg.allsnps=FlexBed(allsnps)
    seg.allsnps.bed_sort()
    seg.includesnps.bed_sort()    
  #############################################################################
  #############################################################################
  print "LOADING PREVIOUS HIGH DEPTH PDIFF GENOME DATA...\n*" , 
  for seg in flatdata:
    seg.highdepth_tmp= "tmp."+rinfo['#VAL']['NAME']+".seg."+seg.name+".highdepth"
    if os.path.isfile(seg.highdepth_tmp):    
      print  seg.name, 
      seg.pdiffs=pickle.load(open(seg.highdepth_tmp, 'rb') )
  print " "
  #############################################################################
  #############################################################################
  if not os.path.isfile(seg.highdepth_tmp):  
    print "PROCESSING POTENTIAL PARALOGOUS DIFFERENCES FOR REPEATS..."
    for seg in flatdata:
      for p in seg.pdiffs.bed:
        features={'rm': None}
        p.append(features)
        #######################################################################
        for r in tabix_rmasker.fetch(p[0],p[1],p[2]):  
          features['rm']=r.split("\t")
          #keep only the repeats that are less than 20% divergence#
  #############################################################################
  if not os.path.isfile(seg.highdepth_tmp):  
    print "PROCESSING POTENTIAL PARALOGOUS DIFFERENCES FOR SNPS..."
    for seg in flatdata:
      for p in seg.pdiffs.bed:
        features=p[6]
        features['snp']=None
        count=0
        #TODO: PROBLEM SNPS ONLY FROM REFERENCE COPY 
        for s in tabix_snps.fetch(p[0],p[1],p[2]):  
          snp=s.split("\t")
          chr_count=2
          if snp[22] != '':
            chr_count=sum(map(float,snp[22].rstrip(",").split(",") ))
          #features['snp']=
          ##need to pick one snp per position##
          ##toss snps only seen in a single genome or two. 
          if chr_count<10:
            continue
          features['snp']=snp
          if int(snp[2])-int(snp[1]) > 1:
            count=count+1
        if count>1:
          sys.exit(0)
  #############################################################################
  if not os.path.isfile(seg.highdepth_tmp):  
    print "CAPTURING POTENTIAL PARALOGOUS DIFFERENCES FOR SUPPORT IN HIGH DEPTH GENOMES..."
    for seg in flatdata:
      for p in seg.pdiffs.bed:
        features=p[6]
        features['highdepth']=[]
      bamcount=0
      for bamfile in highdepth_genomes:
        print bamfile
        bamcount+=1
        bamin=pysam.Samfile(bamfile, 'rb') #print bamin #print bamin.header
        for p in seg.pdiffs.bed:
          features=p[6]
          if features['rm'] != None:
            continue
          for pileupcolumn in bamin.pileup(reference=p[0],start=p[1],end=p[2]):
            #print pileupcolumn
            if pileupcolumn.pos < p[1]:
              continue
            #print pileupcolumn.pos, pileupcolumn.n
            bases={'A':0,'C':0,'T':0,'G':0, 'N':0}
            for pileupread in pileupcolumn.pileups:
              bases[ pileupread.alignment.seq[pileupread.qpos]]+=1
              #break #needed to avoid continuing through pileup
            features['highdepth'].append(bases)
            break
          #print "features", len(features['highdepth']), features['highdepth']
          #print "  READS",counting
            #print "  ",x.tid, bamin.getrname(x.tid)
            #print "     ", x.qname, x.pos, x.aend, "(",x.alen,")
        bamin.close()
      #save a pickle version of all the pdifferences
      print "pickling depth (", seg.highdepth_tmp, ")"
      pickle.dump(seg.pdiffs, open( seg.highdepth_tmp, "wb") )

  #----------------------------------------------------------------------------
  print "DUMPING DATA TO BED FILES (see bed* pdiffs,pindels,include,allsnps) ..."
  print "* ",  
  for seg in flatdata:
    print seg.name,
    basename= "bedraw." + seg.name + "."  
    if seg.has_paralogs:
      seg.pdiffs.writefile((basename+"pdiffs"))
      seg.pdashes.writefile((basename+"pindels"))   
    seg.includesnps.writefile((basename+"include"))
    seg.allsnps.writefile((basename+"allsnps"))
  print
  ############################################################################
  print "CALCULATING PROBABLE PARALOGS..."
  for seg in flatdata:
    for p in seg.pdiffs.bed:
      features=p[6]
      if len (features['highdepth'])==0:
        features['highdepth']=None
        continue
      ## calculate the number of genomes with all differences (usually both)
      paralog_count=0
      pd=p[3].split(":")
      #baseorder=(sorted())
      pbases={}
      for char in pd[2]:
        if char in pbases.keys():
          pbases[char]+=1
        else:
          pbases[char]=1
      pbases_order= sorted (pbases, key=pbases.__getitem__, reverse=True)
      #print pbases_order
      for bases in features['highdepth']:
        #print bases
        all_present=True
        for char in pbases_order:
          if bases[char]==0:
            all_present=False
        if all_present==True:
          paralog_count+=1
      #print paralog_count, " (number of genomes with all paralogs)"
      features['highdepth']=paralog_count
      #print p
      #this is very crude and should be replaced by an algorithm that 
      #trys to predict whether given paralogs carry a gene based on the
      #add qc in terms of whether this fits with the number of known copies#
    
    
  #############################################################################
  print "CALCULATING CUTOFF FOR LIKELY PARALOGY.."
  for seg in flatdata:
    print seg.name,
    #snps should represent rare bases for the most part
    #paralog psvs will represent common for the most part
    paralogpresent = [  x[6]['highdepth'] for x in seg.pdiffs.bed  ]
    #print paralogpresent
    print len(paralogpresent),
    onlyvals= [  x for x in paralogpresent if x is not None ]
    #print onlyvals
    print len(onlyvals),
    if len(onlyvals)>0:    
      seg.pdiffCutoffHighDepth=int( np.median(onlyvals))
    else: 
      seg.pdiffCutoffHighDepth=None
    print seg.pdiffCutoffHighDepth     
    
  
  #############################################################################
  ###### add the includes #####################################################
    #includes= [x for x in rinfo["INCLUDE"]]
    #print includes
  ###MAP COL DATA TARGET, FORCE EXCLUDE ONTO THE THIS DATA
  ### create  seg.col[TARGET] ###
  ### create  seg.col[FORCE] ###
  ### create  seg.col[EXCLUDE] ###
  ### create  seg.col[EXCLUDE] ###
     
  #############################################################################   
  print "GENERATING POTENTIAL MIPS"

  mipinsertsize=140
  ext_arm_range=[16,17,18,19,20]
  lig_arm_range=[24,23,22,21,20]
  ext_arm_max=max(ext_arm_range)
  lig_arm_max=max(lig_arm_range)
  arm_max=max(ext_arm_max,lig_arm_max)
  arm_array_size=xrange(len(ext_arm_range))
  
  for seg in flatdata:
    print seg.name
    seg.mips=[]
    #try to capture
    #print seg.pdiffs
    # print dir(seg.pdiffs)
    #print seg.pdashes
    ## targeted paralogous differences and non targeted regions ##
 
    seg.pdiffs.walker_reset()
    print "PDIFFS",len(seg.pdiffs.bed)
  
    seg.includesnps.walker_reset()
    seg.includesnps.bed_sort()    
    print "INCLUDES", len(seg.includesnps.bed)
    print "   counters:",seg.includesnps.walker_get_counters()
    
    seg.pdashes.walker_reset()
    print "PDASHES",len(seg.pdashes.bed)
    
    seg.allsnps.walker_reset()
    print "ALLSNPS",len(seg.allsnps.bed)
    
    print "REFSEQ",seg.refseqname, seg.refseqbegin, seg.refseqend, len(seg.refdna)
    for i in xrange(len (seg.refdna)):
      if i < arm_max:
        continue
      body_begin= i + seg.refseqbegin
      body_end=body_begin + mipinsertsize
      if body_end + arm_max >= seg.refseqend:
        continue
      ### check the body ####################

      print i, body_begin, '-', body_end, body_end-body_begin
     
     ### move counters ###
      psvs=seg.pdiffs.walker_get_range(seg.refseqname,body_begin,body_end, counter=0)
      #psvs=remove_psvs(psvs,withRepeats=False,withSnps=False,minHighDepth=None)
      includes=seg.includesnps.walker_get_range(seg.refseqname,body_begin,body_end,counter=0)      
      print "PSVS",len(psvs), "INCLUDES", len(includes)
      if len(psvs)==0 and len(includes)==0:
        continue
      pgaps=seg.pdashes.walker_get_range(seg.refseqname,body_begin,body_end,counter=0)
      print  "GAPS", len(pgaps)
      number_gaps=1
      max_size_gaps=2
      gapbase=1  #TODO: calculate gap bases
      if len(pgaps)>0 and gapbases < max_size_gaps:
        continue
      #TODO: CHECK FOR OVERLAP WITH EXCLUDE REGIONS
      #check for anything falling within a region#
      #check for anything falling within a 
      
      print psvs
      for strand in ["+","-"]:
         print strand
      ### check for a potential target arms
      ### calculate score
      ### check if quality is good enough 
      
      ###walk to appropriate position
      
      #keep=False

      ######################################
            
      
      #for v in arm_array_sizes:
      
      
      #sys.exit()
      #walk all the walkers to the 
      #counter 0, is for the body
      #counter 1, is for the left arm
      #counter 2 is for the right arm
      #check if pdiffs or includes are within
      #pdiffs_fbed.walker_get_range(ref_chr,body_chr_begin, body_chr_end,counter=1)
        # set keep=True
        #check if pdashes or excludes are within
        # set keep=False
        
        #design algorithm (2 mips + and - strand)        
        #choose the best one on each strand based on
        #scoring system
        

 #       pdbed=pdiffs_fbed.walker_get_range(ref_chr,body_chr_begin, body_chr_end,counter=1)
 #       tsbed=targetsnps_fbed.walker_get_range( ref_chr, body_chr_begin, body_chr_end, counter=1)
  #      if pdbed or tsbed:
  #        print pdbed
   #       print tsbed
  #      else:
  #        continue
        
        #--is there anything worthwhile?
        #--is there a PSV?
        #--is there a SNP?
        #if (not (pnd or snp)) then:
          #  continue
        # if a potential mip then toss
        ### evaluate the different potential ends ###
        ### score the different arms and keep one ###
        ### check for SNPs or excludes within arms ###
        ### yes/no  (larry corrie)
        #print rseq[body_begin:body_end]
        #print len(rseq[body_begin:body_end] )
        
        
        ### check the body ###
        # check for things to exclude#
        # skip if in the body ###
        
        #sys.exit(0)
        #calculate body--see if it is encompassing anything useful#
        #is it of any interest?#
          #next if nothing of interest then toss#
        
        #if body of interest then check if any reason to toss it?
            #e.g. more than a single indel
        
        ###make sure the sequence is not horrible###
        
        #pick best F and R based on melting scores#
          #check SNPs or indels within arms if so bail out#
        
        ###search regions####
        
        #check body
        #allow single SNP in one end in either ligation or extension#
        #check the melting temp
        #check the ends#
        #names all the same check both strands.
        
      
      
      
      
      
      
      #combine excludes
      #
      
      #badprimer=
      #dout=multiple_alignment_traverse(marray, regions=region_coords,
          #  seq_coord=0)
      #have array and dictionary
      
    pysam_fgenome.close()
    
  
  
  #this analysis will be at the duplicon level 
  #alignments may be segmented but should be in linear order
  #with regard to the first sequence which is the master
  #the master should usually be the one with the most insertions
  #but all the data will be analyzed in all orientations. 
  ###load the segments from the multiple alignments
  ###run initial scan masking out bad regions###
  ###ordered segments 0,1,2,3,4,5
  ### data structures###


  

##############################################################################
def design_mips(mip):
  #given the beginnings of a mip, finish selecting the best arms, 
  return



def config_file_load (infile_path):
  """ load a standard config file 
  the configuration file allows comments with #
  the key value pairs are separated by an = sign
  """
  config={}
  #print type(config)
  pattern=re.compile(r"\s+")
  inreader=open(infile_path,'r')
  for line in inreader:
    #print "L",line
    line=line.strip()
    line=pattern.sub('',line,0)
    if line =='' or line.startswith("#"):
      continue
    #print "l",line
  
    line=re.split(r"=",line)
    if len(line)!=2:
      print >> sys.stderr, line
      print >> sys.stderr, "Improper variable and value pair!"
      sys.exit(1)
    name=line[0].lower()
    config[name]=line[1]
    #print "config:", name, line[1]
  return config
        
###############################################################################
def info_file_load ( infile_path):
  """Helper function to load need region information into a dictionary
  This file contains:
  *regions to target with MIPs
  *paralogs to align
  *SNPs,INDELs and other MIPs to include
  """
  line_count=0
  legalsingleparameter=['NAME','DESC']
  legalmultiparameter=['REGION','FORCE','TARGET','EXCLUDE']
  rinfo = recursivedefaultdict()
  regex=re.compile(r"\s+")
  inreader=open(infile_path,'r')
  lines=inreader.readlines()
  for l in lines:
    line_count+=1 #count lines for error reporting
    #print line_count
    l=l.strip( ) #remove any spaces
    if re.match(r"^#",l) or l == '':
      continue
    c=re.split(r"\s+",l,0)
    if c[0] in legalsingleparameter: #check if col zero is a legal parameter (single value)
      rinfo['#VAL'][c[0]]=c[1]
      #print "legal %s with value (%s)\n " % (c[0],c[1])
    else:
      #parse first column to look for legal multivalue (array data)
      level=re.split(r":",c[0])
      #print level
      if level[0] == 'COL':
        if level[1] in legalmultiparameter:
          rinfo['COL'][ level[1] ]=c[1:]
          #print c[1:]
          if level[1]!='REGION':
            rinfo[level[1]]=[]
        else:
          print >> sys.stderr, "Line:%s Error, %s , not legal col argument\n"  % (line_count, l[1] )
          sys.exit(1)
      elif level[0] in legalmultiparameter:
        if not level[0] in rinfo['COL']:
          print >> sys.stderr, "Data type (", level[0], ") has no COL header information!"
          sys.exit(1)
        d= dict (zip (rinfo['COL'][level[0]], c[1:]) )
        #print d
        info_clean_line(d,line_count)
        #print d
        if level[0] == 'REGION':
          rinfo['REGION'][level[1]][level[2]]=d
        else:
          #print level[0],  c[1:]
          l=rinfo[level[0]]
          l.append(d)
      else: 
        print >> sys.stderr, "Line:%s Illegal option (%s) \n" % (line_count,level[0])
        print >> sys.stderr,  "Legal info single value input:", legalsingleparameter
        print >> sys.stderr, " Legal info multiple value input:", legalmultiparameter
        sys.exit(1)
  return rinfo

def info_clean_line(d,line_count):
  if d['orient'] in ['F','f','P','+']:
    d['orient']="+"
  elif d['orient'] in ['R','r','N','-']: 
    d['orient']="-"
  else:
    print >> sys.stderr , "Line:%s Error: the orient(%s) is not legitimate" \
        % (line_count, d['orient'])
    sys.exit(1)
  d['seqb']=int(d['seqb'])
  d['seqe']=int(d['seqe'])
  d['seqb']=int(d['seqb'])
  if 'length' in d:
    d['length']=int(d['length'])
  return d

##############################################################################  
def run_command (command,shell=True):
  proc = subprocess.Popen(command,stdin=None,stdout=subprocess.PIPE, \
            stderr=subprocess.PIPE,shell=shell)
  (odata,edata) = proc.communicate()
  #print odata
  if len(edata)>0:
      edata=0
  return [edata,odata]


#def indexfa_grab_seq_region (faidx,chrom, start, end, orient):
#  """this is in base 1,1 coordinates"""
#  region=chrom+":"+str(start)+"-"+str(end)
#  seqtext=faidx.fetch(region= region)
#  seq=Seq(seqtext,generic_dna)
#  if orient == '-':
#    seq=seq.reverse_complement()
#  rec=SeqRecord(seq,region, '', '')
#  return rec

def indexfa_grab_seq (faidx,chrom,start,end,orient):
  seqtext=faidx.fetch(reference=chrom,start=start,end=end)
  seq=Seq(seqtext,generic_dna)
  if orient == '-':
    seq=seq.reverse_complement()
  region=chrom+":"+str(start)+"-"+str(end)
  rec=SeqRecord(seq,region, '', '')
  return rec  
  
### subroutines 
def dict_count_keywithmaxval(d):
  v=list(d.values())
  k=list(d.keys())
  return k[v.index(max(v))]
def dict_count_keywithminval(d):
  v=list(d.values())
  k=list(d.keys())
  return k[v.index(min(v))]
  

###############################################################################
###############################################################################
def remove_psvs(psvs,withRepeats=False,repeatDiv=None,withSnps=False,minHighDepth=None):
  """This script is dependent on psv data structure
   withRepeats -- set to True to remove all Repeats 
   repeatDiv -- set minDiv to keep
   withSNPs -- set to True to toss any called SNPs underlying 
   minHighDepth -- number of genomes from highdepth that required to call
  """
  print "REMOVE",psvs
  for p in list(psvs):
    print p[6]['highdepth']
    if withRepeats and p[6]['rm']!=None:
      psvs.remove(p)
    elif repeatDiv!=None and int(p[6]['rm'][1]) < repeatDiv:
      psvs.remove(p)
    elif withSnps and p[6]['snp']!=None:
      psvs.remove(p)
    elif minHighDepth != None and (p[6]['highdepth']==None or int(p[6]['highdepth'])< minHighDepth ):
      psvs.remove(p)
  return psvs
    

###############################################################################
###############################################################################

#------------------------------------------------------------------------------
class recursivedefaultdict(defaultdict):
  """ Allows for rapid creation without initializaiton of
    multidepth  dictionary datastructures """
    
  def __init__(self):
    """ the only function to steal everything from default dictionary"""
    self.default_factory = type(self)  
  

### class pairwise

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
  note: a bed is only defined seq:start:end:orient  (no more required)
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
    for i  in self.bed:
      outline=outfile.write("\t".join([str(x) for x in i])+ "\n")
    outfile.close()
    return
  def walker_step_upto(self, seq, start, counter=0):
  #  for i in xrange (self.current_positions[counter],len(self.bed)):
  #    feature=self.bed[i]      
  #    if seq!=feature[0]: 
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
      if seq==feature[0] :  #bed feature on same seq as our range
          if feature[2] < start: #bed feature less than start of range
            if autoincrement==True:
              self.current_positions[counter]=i+1 #so increment counter
            continue #and go to next
          newbegin=max(start,feature[1]) 
          newend=min(end,feature[2])
          if  newend-newbegin  > 0: #bed feature within our range
            newbed.append(feature)
            continue
          else:   # query feature << bed feature
            break # already past it so stop
      if seq < feature[0]:   #stop and bail  current seq position is greater than seq
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
               start=1, stop=2,  name=3,score=4,orient=5):
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
  #    bed.check(check='chr',checknum=True) 
  #def bed.readlist(self,tab_delimited_list):
    #"""parse list of tab-delimited text into bed list of lists"""
  #tablist.strip()
    #return 
  #def bed.check(tab_delmited_bed_list):
###THIS is now used in my code ### 
class  Interval:
  """ ucsc genome coordinates [a,b) """
  @staticmethod
  def intersect2interval(begin1,end1, begin2,end2):
    "intersect [a,b)  [c,d) returns new begin and end of intersect"
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
    self.reverseNames={}  #pushes list positions
    self.forwardNames={}  #pushes list positions
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
    """NOT IMPLEMENTED:  could check consistency of the data"""
  def forwardMap(self, seq,begin,end=None):
    """given seq1,begin1,[end1], return seq2 equilvalent"""
    return
  def reverseMap(self,seq,begin,end=None):
    """given se2,begin,[end2], lookup seq1 equilvaent..."""
    return


#-----------------------------------------------------------------------------
class MultipleAlign():
  """A more functional multiple alignment class wrapping Bio.MultipleAlignment
  """  
  def __init__(self):
    self.position=0
    self.malign=None #stores the Bio.AlignIO format 
    self.seqstr=None
    self.internalCoords=None  
    self.externalCoords=None
    self.internalLookup=[] #array of hashes seq position to alignment position. 
    self.externalLookup=[]
    self.externalSeqs=None
    self.externalNames=None
  def loadAlignment(self,filepath,alignformat):
    "loads a biopython alignment"
    self.malign=AlignIO.read(filepath, alignformat)
    #tmplen=self.malign.get_alignment_length()
    #for s in self.malign:
    #  print len(str(s.seq))
    #  self.seqstr.append(str(s.seq))
    #print tmplen
    return
  def generateInternalCoords(self,debug=False):
    "generates coordinates for the alignment"
    self.internalCoords=[]
    self.internalLookup=[]
    if debug:
      print "DEBUGGING"
    for s in self.malign:
      print s.name
      count=0
      apos=0
      positions=[]
      dictlookup={}
      for i in s.seq:
      
        if i in ['A','C','T','G', 'N']:
          positions.append(count)
          dictlookup[count]=apos      
          count+=1
          
        elif i in ['-']:
          positions.append(count-0.5)
        else:
          print >>sys.stderr, i, "is not ACTGN-"
          sys.exit(1)
        
       # if debug and count >33333:
       #     print i, count,
        #print apos, positions[-1], i
        apos+=1
      #print "POSLEN",len (positions)
     # print "DICT", len (dictlookup)
      #print positions[-30:]
      self.internalCoords.append( positions)
      self.internalLookup.append( dictlookup)
    return
  def generateExternalCoords(self, externalRegions):
    """"list [chr,b,e,orient,name], [ ...], .... to put alignment in context.
    Order of alignment positions."""

    self.externalCoords=[]
    self.externalLookup=[]
    self.externalSeqs=[]
    self.externalNames=[]
    if len(externalRegions) != len(self.malign): 
      print >> sys.stderr, "coordinate regions not equal to size of alignment"
      sys.exit(1)
    for i in range(len(externalRegions)):
      seq=externalRegions[i][0]
      begin=externalRegions[i][1]
      end=externalRegions[i][2]
      orient=externalRegions[i][3]
      name=externalRegions[i][4]
      #print seq,begin,end, orient, name
      positions=[]
      dictlookup={}
      if orient=='+':
        positions=[ (pos+begin) for pos in self.internalCoords[i] ]
        #print len(positions)
        for k,v in self.internalLookup[i].iteritems():
          extcoord=k+begin
          dictlookup[k+begin]=v
      elif orient=='-':
        positions=[(end-pos-1) for pos in self.internalCoords[i] ]
        for k,v in self.internalLookup[i].iteritems():
          dictlookup[end-k-1]=v
      else:
        print >> sys.stderr, "orientations ", orient, "is not + or -!"
        sys.exit(1)
      self.externalNames.append(name)
      self.externalSeqs.append(seq)
      self.externalCoords.append(positions)
      self.externalLookup.append(dictlookup)
     
      
  def getSeqCoord(self,seqnum,position,coord='internal'):
    """Input alignment seqnum and position
    Returns: sequence coordiante
    coord='internal' is default 'external'
    """
    if coord=='external':
      #print self.externalCoords[seqnum]
      return self.externalCoords[seqnum][position]
    elif coord=='internal':
      return self.internalCoords[seqnum][position]
    else:
      print >>sys.stderr, "coord=",coord," must be 'internal' or 'external'!"
      sys.exit(1)
  def getSeqBase(self, seqnum, position, coord='internal'):
    "input: seqname and position, return base characterstring"
    return self.malign [seqnum,position]
  def getAlignColumn(self,position, coord='internal'):
    "input: position coord=internal/external output: string of bases"
    return self.malign[:,position]
  #def getExternalCoord(self,seqname,position):
  def getDifferences(self, seqout=-1, external=False, types="base:dash"):
    """returns a bed like list of positions where the multiple alignment
    differs between the bases
    seq begin end alignement_info(a:alignposition#:bases  """
    seqnum=self.getSeqCount()
    data=[]
    nodash=True
    nobase=True
    if 'base' in types:
      nobase =False
    if 'dash' in types:
      nodash=False
    if nodash==True and nobase==True:
      print >>sys.stderr, "types for get Differences excludes both base and dash!"
      sys.exit(1)
    for i in range(   self.getAlignLength()   ):
      column=self.getAlignColumn(i)
      if seqnum!=column.count(column[0]): #checks if first is same
        if column.count("-")>0: #there is a dash
          if nodash:
            continue
        else: #there is no dash (presumably base count)
          if nobase:
            continue
        seqname='align'
        #print column
        seqbegin=i
      
        score=1
        
        name="a:"+str(i)+":"+column
        if seqout>=0:
          seqname=seqout
          seqbegin=self.internalCoords[seqout][i]
        if external==True:
          seqname=self.externalSeqs[seqout]
          seqbegin=self.externalCoords[seqout][i]
        suns=[]
        for s in range (len(column)):
          if column.count(column[s])==1:
            suns.append(str(s))
        name+=":"+",".join(suns)
        data.append([seqname,seqbegin,seqbegin+1,name,  score, '+' ])
    return data
  def ext2extInterval(self,incopynum,inchr, inbegin,inend,outcopynum):   
    (c,b)= self.ext2extPosition(incopynum, inchr,inbegin, outcopynum )
    (c,e)= self.ext2extPosition(incopynum, inchr,int(inend)-1, outcopynum )
    if b > e:
      (b,e)=(e,b+1)
    else:
      e=e+1
    return c,b,e
  def ext2extPosition(self, incopynum, inchr, inposition,outcopynum):
    #TODO: Check for text --need int
    #TODO: Need check on the chromosome
    #print "SELF",self.externalSeqs[incopynum]
    #print incopynum, inchr, inposition, outcopynum
    #print self.externalLookup[incopynum]
    internalposition=self.externalLookup[incopynum][int(inposition)]
    outposition=self.externalCoords[outcopynum][internalposition]
    #return  self.externalSeqs[outcopynum], outposition
    return self.externalSeqs[outcopynum],outposition
  def printer(self):
    x=str(self.malign)
    print x
    return
  def getAlignLength(self):
    return len(self.malign[0])
  def getSeqCount(self):
    return len(self.malign)

#------------------------------------------------------------------------------
class Segment:
  def __init__(self,name):
    self.__dict__['initialized']=0
    self.name=name  #descriptive name for segment 
    self.has_paralogs=False #if more than one copy then True
    self.copies=[] # this is the ... names of the copyies???
    self.copydata={}  ### deep dictionary from rinfo with names
    self.clustal_infile=None
    self.clustal_outfile=None
    self.clustal_outfile=None
    self.highdepth_tmp=0 ##
    self.refdna='' ##actual dna seq of the reference paralog
    self.refseqname=''
    self.refseqbegin=-1
    self.refseqend=-1
    self.refseqcopyname=''
    self.msalign=None #multiple sequence alignment   
    self.pdiffs=None
    self.pdiffCutoffHighDepth=0
    self.pdashes=None
    self.mips=None
    self.allsnps=None
    self.othersnps=None
    self.includesnps=None
    self.initialized=1
  def __setattr__(self,k,v):
    if  self.initialized==1:
      if k in self.__dict__.keys():
        self.__dict__[k]= v
      else:
        raise AttributeError, "Attempt to set a non-init attribute:"+k
    else: 
      self.__dict__[k]= v
#------------------------------------------------------------------------------



class Mip():
  __slots__ = [ 'name', 'strand', 'extarmSeq' ]
  def __init__(self):

    self.name=name
    self.strand
    self.extarmSeq
    self.extarmBegin
    self.extarmEnd
    self.extarmCopy
    self.ligarmSeq
    self.ligarmBegin
    self.ligarmEnd
    self.ligarmCopy
    self.bodySeq
    self.bodyBegin
    self.bodyEnd
    self.featureBegin
    self.featureEnd
    self.featureScore
 


###############################################################################
###############################################################################
#------------------------------------------------------------------------------
#main is left to the bottom so that global variables can be populated
if __name__ == "__main__":
  if len (sys.argv)==1:
    sys.argv.append("-h")  #if no command then it is a cry for help
  main(sys.argv[1:])





