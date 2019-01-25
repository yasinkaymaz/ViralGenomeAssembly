#!/usr/bin/python
import os
dir = os.path.dirname(__file__)
import argparse
import glob
import sys
import textwrap

import subroutines as sr

argparser = argparse.ArgumentParser(description='Process some integers.')
subcommands={}

def main(args):
    if len(args)  == 0 or args[0] in ["h", "help", "-h", "--h", "--help","-help"] :
        verbosity= 'shortDesc'
        if args[0] in ["help" , "--help", "-help"]:
            verbosity = 'longDesc'
        program_name=os.path.basename(__file__)
        print "\n","USAGE:",program_name, "[-h] subcommand [suboptions]","\n", "-----"
        print "DESCRIPTION: A collection of tools to analyze genetic diversity and distances between whole genome sequences."
        print "-----------"
        print "SUBCOMMANDS:"
        print "-----------"
        for k in subcommands.keys():
            text=subcommands[k][verbosity]
            text= textwrap.dedent(text)
            if text:
                text =  "%s:   %s " %(k, text )
                print textwrap.fill(text, 75, initial_indent='', subsequent_indent='                    ')
        print "\n","HELP:"
        print "----"
        print "-h/-help   short / long  subcommand descriptions","\n"
        print "For specific options:", program_name, "[subcommand] --help", "\n"
    elif args[0] == 'pydoc':
        os.system( "pydoc " + os.path.abspath(__file__) )
    elif args[0] in subcommands.keys():
        globals()[args[0]](args[1:])
    else:
        print "unknown subcommand (" + args[0] + ") use -h for list of subcommands!"
        sys.exit(-1)
    sys.exit(0)

    return


shortDescText="Calculate average genetic distance of a given msa alignment."
longDescText="""This function calculates genetic distance using Kimura 2-parameter method."""
subcommands['MeanPairwiseKimuraDist'] = { 'shortDesc':shortDescText, 'longDesc': longDescText }

def MeanPairwiseKimuraDist(args):
    argparser.add_argument("-af", "--alignmentFile", required=True, help='A multiple sequence alignment file in fasta format. Please provide with full directory.')
    argparser.add_argument("-rn", "--regionName", required=True, help='A multiple sequence alignment file in fasta format. Please provide with full directory.')
    args=argparser.parse_args(args=args)
    alignmentFile = glob.glob(args.alignmentFile)[0]
    regionName = args.regionName

    dist = sr.DNAmsaDist(alignmentFile, regionName)
    print dist.MeanPairwiseKimuraDist()
    return



shortDescText="Calculate average dN/dS ratio of a given msa alignment."
longDescText="""This function calculates dN/dS ratio."""
subcommands['MeandNdS'] = { 'shortDesc':shortDescText, 'longDesc': longDescText }

def MeandNdS(args):
    argparser.add_argument("-af", "--alignmentFile", required=True, help='A multiple sequence alignment file in fasta format. Please provide with full directory.')
    argparser.add_argument("-rn", "--regionName", required=True, help='A multiple sequence alignment file in fasta format. Please provide with full directory.')
    args=argparser.parse_args(args=args)
    alignmentFile = glob.glob(args.alignmentFile)[0]
    regionName = args.regionName

    dist = sr.DNAmsaDist(alignmentFile, regionName)

    print dist.MeandNdS()
    return




if __name__ == "__main__":
    if len (sys.argv)==1:
        sys.argv.append("--help")  #if no command then it is a cry for help
    main(sys.argv[1:])
