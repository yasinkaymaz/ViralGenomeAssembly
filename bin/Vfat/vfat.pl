#!/share/pkg/perl/5.18.1/bin/ perl -w

# Copyright © 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

use strict;
use Getopt::Long;
my %option = (
	minlongcont		=> 100,
	mincontlen 	=> 350,		# Minimum length of contigs that will be analyzed
	readfa			=> '',		# Input file containing the reads in fasta format
	readq			=> '',		# Input file containing the read quals
	readfa2				=> '',		# Input file containing the reads in fasta format (2nd mate if paired)
	readq2				=> '',		# Input file containing the read quals  (2nd mate if paired)
	readfq				=> '',		# Input file containing the reads in fastq format
	readfq2				=> '',		# Input file containing the reads in fastq format (2nd mate if paired)
	maxseggap		=> 30,		# contigMerged.pl opt for max gap length before splitting segments
	maxsegins		=> 60,
	minseglen		=> 50,		# contigMerged.pl opt for min length of a segment
	maxsegdel		=> 0.25,
	longcont			=> 1500,
	contigcov			=> 2,
	allcont				=> '',
	sequencer		=> 'illumina',
	genelist		=> '',
	pepfolder		=> '',
	ref				=> '',
	contigs			=> '',
	o				=> '',
	maxannogaplen	=> 50, #annotate.pl option for max gap length before splitting gene
	minhomosize	=> 3,
	readwindow	=> 5,
	nofix		=> '',
	fixhomo		=> '',
	forcefix	=> '',
	virus		=> '',
	amps		=> '',
	h			=> '',
	minpctid	=> 0.4,
	minpctlen	=> 0.2,
	fakequals	=> 30,
	details			=> '',
	noannot		=> '',
	sampname	=> '',
);

GetOptions(
	"minlongcont=i"		=> \$option{minlongcont},
	#"mincontiglen=i"	=> \$option{mincontiglen},
	"mincontlen=i"       => \$option{mincontlen},
	"readfa=s"			=> \$option{readfa},
	"readq=s"				=> \$option{readq},
	"readfa2=s"			=> \$option{readfa2},
	"readq2=s"			=> \$option{readq2},
	"readfq=s"			=> \$option{readfq},
	"readfq2=s"			=> \$option{readfq2},
	"maxseggap=i"		=> \$option{maxseggap},
	"maxsegins=i"		=> \$option{maxsegins},
	"minseglen=i"		=> \$option{minseglen},
	"maxsegdel=f"		=> \$option{maxsegdel},
	"longcont=i"			=> \$option{longcont},
	"contigcov=i"		=> \$option{contigcov},
	"allcont"				=> \$option{allcont},
	"sequencer=s"		=> \$option{sequencer},
	"genelist=s"			=> \$option{genelist},
	"pepfolder=s"		=> \$option{pepfolder},
	"ref=s"					=> \$option{ref},
	"contigs=s"			=> \$option{contigs},
	"o=s"					=> \$option{o},
  	"maxannogaplen=i"	=> \$option{maxannogaplen}, 
 	"minhomosize=i"		=> \$option{minhomosize},
 	"readwindow=i"		=> \$option{readwindow},
 	"nofix"				=> \$option{nofix},
  	"fixhomo"			=> \$option{fixhomo},
  	"forcefix=f"		=> \$option{forcefix},
	"virus=s"			=> \$option{virus},
	"amps=s"			=> \$option{amps},
 	"h"					=> \$option{h},
	"minpctid=f"		=> \$option{minpctid},
	"minpctlen=f"	=> \$option{minpctlen},
  	"fakequals=i"		=> \$option{fakequals},
	"details"			=> \$option{details},
	"noannot"		=> \$option{noannot},
	"sampname=s"		=> \$option{sampname},
);

### DEPENDENCIES AT BROAD:
### setenv WISECONFIGDIR /seq/annotation/bio_tools/GeneWise/wise2.2.0/wisecfg/
### use R-2.9

my $scriptpath = "/project/umw_jeffrey_bailey/share/bin_sync/Vfat/";
my $refDataPath = "/project/umw_jeffrey_bailey/share/EBV/Annotation/Type1/Vfat/EBV/";
my $mosaikParam = " -hgop 20 -gop 40 -gep 10 -bw 29 -st illumina";
if($option{sequencer} eq '454')
{
	$mosaikParam = " -hgop 4 -gop 15 -gep 6.66 -bw 0 -st 454";
}
if($option{fakequals})
{
	$mosaikParam .= " -fakequals ".$option{fakequals};
}

if($option{h})
{
	print "vfat.pl - Post-Assembly toolkit for orienting, merging, frameshift correcting, annotate and calculate statistics on assembly data\n\n";

	print "Scripts requirement : orientContig.pl, contigMerger.pl, fixFrameshifts.pl, annotate.pl, runMosaik2.pl, QA_stats.pl\n\n";

	print "Paths requirement (either have these in your environment variables or specify through configPaths.pl):\nMosaik v2.1.33 or higher\n";
	print "R-2.9 or higher\n";
	print "Genewise (wise2.2)\n";
	print "The location of all the scripts, muscle, and the reference data (if using the -virus option) can be modified if required in the scripts headers in the appropriate variables\n\n";

	print "Command lines ([] means optional):\n";
	print "Without -virus :\nperl vfat.pl -contigs <contigs.fa> [-readfa <reads.fa> -readq <reads.qual>] -ref <ref.fa> -genelist <ref_genelist.txt> -amps <ref_amplicons.txt> -pepfolder <peptidesFolder> -o <outputBasename>\n";
	print "With -virus :\nperl vfat.pl -contigs <contigs.fa> -virus <virus> [-readfa <reads.fa> -readq <reads.qual>] -o <outputBasename>\n\n";

	print "Command line switches (see doc for specific file formats):\n\nInput Files:\n";
	print "-contigs\tContigs file in multi-fasta format\n";
	print "-ref\t\tReference file in fasta format\n";
	print "-readfa\t\tReads file in fasta format\n";
	print "-readq\t\tReads quality file in fasta format (spaced integers)\n";
	print "-readfa2\t\tReads file in fasta format (2nd mate if paired)\n";
	print "-readq2\t\tReads quality file in fasta format (spaced integers) (2nd mate if paired)\n";
	print "-readfq\t\tReads file in fastq format\n";
	print "-readfq2\t\tReads file in fastq format (2nd mate if paired)\n";
	print "-genelist\tReference gene positions list (see doc for format)\n";
	print "-amps\t\tReference amplicons positions list (see doc for format)\n";
	print "-pepfolder\tFolder containing fastas of gene peptides (see doc for format)\n";
	print "-o\t\tOutput base name, all output files will begin with this\n";
	print "-virus\t\tAllows you to specify a virus name instead of the reference, genelist, amps and pepfolder options (See doc for details)\n";
	print "\t\tAvailable viruses for you are:";
	if(-d $refDataPath)
	{
		opendir(REFDATA, $refDataPath);
		while(my $file = readdir(REFDATA))
		{
			if($file eq '.' || $file eq '..'){next;}
			if(-d $refDataPath.$file)
			{
				print " $file";
			}
		}
		print "\n\n";
	}else{
		print " None, need to setup your reference data path (see doc)\n\n";
	}
	print "Options (default value in bracket, if any):\n";
	print "-sequencer(illumina)\tGeneral parameter. Sets which sequencer was used to determine which version of Mosaik to run. 	Currently supports illumina, 454\n";
	print "-fakequals()\trunMosaik2.pl parameter. Fake all quality scores to a given value in qlx files. This won't affect the results of any script in this package and will speed it up, but qlx files will not be valid for other softwares using them like V-Phaser\n";
	print "-mincontlen(350)\torientContig.pl filter. Minimum length of contigs that will be analyzed\n";
	print "-minlongcont(100)\torientContig.pl filter. Minimum length for the longest continuous stretch of the contig aligning to reference\n";
	print "-maxseggap(30)\t\tcontigMerger.pl parameter. Maximum gap length before splitting segments\n";
	print "-maxsegins(60)\t\tcontigMerger.pl parameter. Maximum insertion length before splitting segments\n";
	print "-minseglen(50)\t\tcontigMerger.pl parameter. Minimum length of a segment\n";
	print "-maxsegdel(0.25)\t\tcontigMerger.pl parameter. Max % of deletion allowed in a segment\n";
	print "-longcont(1500)\t\tcontigMerger.pl parameter. Min length to be a 'long contig', which will always be considered\n";
	print "-contigcov(2)\t\tcontigMerger.pl parameter. Contigs will be considered as long as you have a coverage below contigcov\n";
	print "-allcont\t\tcontigMerger.pl parameter. Analyze all contigs, regardless of other filters\n";
	print "-minhomosize(3)\t\tfixFrameshifts.pl parameter. Minimum length to consider a string of bases an homopolymer. Recommend 2 for 454, 3 for Illumina\n";
	print "-readwindow(5)\t\tfixFrameshifts.pl parameter. Length of sequence on each side of a frameshift to look at in reads to find a non-frameshifted correction\n";
	print "-nofix\t\t\tfixFrameshifts.pl parameter. Using -nofix will still look for frameshifts and write the log, but it will not automatically correct them in the assembly\n";
	print "-fixhomo\t\tfixFrameshifts.pl parameter. Forces correction of frameshifts in homopolymer regions. Will work without reads supplied, but requires at least 1 read support to fix if reads are supplied\n";
	print "-forcefix\t\tfixFrameshifts.pl parameter. Will fix frameshifts as long as one of the correct length windows is present above # (supplied, i.e. 0.25) fraction of the reads. Off by default\n";
	print "-maxannogaplen(50)\tannotate.pl parameter. Maximum gap in the alignment before splitting gene\n";
	print "-minpctlen(0.2)\t\tannotate.pl parameter. Minimum % length of a gene aligning to consider it partial instead of missing\n";
	print "-minpctid(0.4)\t\tannotate.pl parameter. Minimum % identity of a gene aligning to consider it partial instead of missing\n";
	print "-sampname()\t\tannotate.pl parameter. Specifies the sample name to write in NCBI submission files. If none is specified, uses the output name (-o)\n";
	print "\nSome option wording might be confusing due to length constraints. See full documentation for more details\n";
	
	exit();
}



if($option{virus})
{
	unless($option{ref})
	{
		$option{ref} = $refDataPath.$option{virus}."/".$option{virus}."_Reference.fasta";
	}
	unless($option{genelist})
	{
		$option{genelist} = $refDataPath.$option{virus}."/".$option{virus}."_Reference_genelist.txt";
	}
	unless($option{amps})
	{
		$option{amps} = $refDataPath.$option{virus}."/".$option{virus}."_Reference_amplicons.txt";
	}
	unless($option{pepfolder})
	{
		$option{pepfolder} = $refDataPath.$option{virus}."/".$option{virus}."_Peptides";
	}
}

my $inputcontIn = $option{contigs};
my $refIn = $option{ref};
my $output = $option{o};

open(COMMANDLOG, ">".$output."_commandlog.txt");

### ASSIGN ALL OPTIONS TO COMMAND LINES
my $optCont2Ass = '';
my $orientcontOpt = '';
my $contigmergerOpt = '';
my $fsOpt = '';
my $annoOpt = '';
my $QAopt = '';


if($option{minlongcont} != 100)
{
	$optCont2Ass .= " -minlongcont ".$option{minlongcont};
	$orientcontOpt .= " -minlongcont ".$option{minlongcont};
}
if($option{mincontlen} != 350)
{
	$orientcontOpt .= " -mincontlen ".$option{mincontlen};
	$optCont2Ass .= " -mincontlen ".$option{mincontlen};
	$contigmergerOpt .= " -mincontlen ".$option{mincontlen};
}
if($option{maxseggap} != 30)
{
	$contigmergerOpt .= " -maxseggap ".$option{maxseggap};
	$optCont2Ass .= " -maxseggap ".$option{maxseggap};
}
if($option{maxsegins} != 60)
{
	$contigmergerOpt .= " -maxsegins ".$option{maxsegins};
	$optCont2Ass .= " -maxsegins ".$option{maxsegins};
}
if($option{minseglen} != 50)
{
	$contigmergerOpt .= " -minseglen ".$option{minseglen};
	$optCont2Ass .= " -minseglen ".$option{minseglen};
}
if($option{maxsegdel} != 0.25)
{
	$contigmergerOpt .= " -maxsegdel ".$option{maxsegdel};
	$optCont2Ass .= " -maxsegdel ".$option{maxsegdel};
}
if($option{longcont} != 1500)
{
	$contigmergerOpt .= " -longcont ".$option{longcont};
	$optCont2Ass .= " -longcont ".$option{longcont};
}
if($option{contigcov} != 2)
{
	$contigmergerOpt .= " -contigcov ".$option{contigcov};
	$optCont2Ass .= " -contigcov ".$option{contigcov};
}
if($option{allcont})
{
	$contigmergerOpt .= " -allcont ".$option{allcont};
	$optCont2Ass .= " -allcont ".$option{allcont};
}

if($option{readfa})
{
	$contigmergerOpt .= " -readfa ".$option{readfa};
	$optCont2Ass .= " -readfa ".$option{readfa};
}
if($option{readq})
{
	$contigmergerOpt .= " -readq ".$option{readq};
	$optCont2Ass .= " -readq ".$option{readq};
}
if($option{readfa2})
{
	$contigmergerOpt .= " -readfa2 ".$option{readfa2};
	$optCont2Ass .= " -readfa2 ".$option{readfa2};
}
if($option{readq2})
{
	$contigmergerOpt .= " -readq2 ".$option{readq2};
	$optCont2Ass .= " -readq2 ".$option{readq2};
}
if($option{readfq})
{
	$contigmergerOpt .= " -readfq ".$option{readfq};
	$optCont2Ass .= " -readfq ".$option{readfq};
}
if($option{readfq2})
{
	$contigmergerOpt .= " -readfq2 ".$option{readfq2};
	$optCont2Ass .= " -readfq2 ".$option{readfq2};
}
if($option{sequencer} ne 'illumina')
{
	$contigmergerOpt .= " -sequencer ".$option{sequencer};
	$optCont2Ass .= " -sequencer ".$option{sequencer};
}
if($option{fakequals} != 0)
{
	$contigmergerOpt .= " -fakequals ".$option{fakequals};
	$optCont2Ass .= " -fakequals ".$option{fakequals};
}
if($option{readwindow} != 5)
{
	$fsOpt .= " -readwindow ".$option{readwindow};
	$optCont2Ass .= " -readwindow ".$option{readwindow};
}
if($option{minhomosize} != 3)
{
	$fsOpt .= " -minhomosize ".$option{minhomosize};
	$optCont2Ass .= " -minhomosize ".$option{minhomosize};
}
if($option{nofix})
{
	$fsOpt .= " -nofix";
	$optCont2Ass .= " -nofix";
}
if($option{fixhomo})
{
	$fsOpt .= " -fixhomo";
	$optCont2Ass .= " -fixhomo";
}
if($option{forcefix})
{
	$fsOpt .= " -forcefix ".$option{forcefix};
	$optCont2Ass .= " -forcefix ".$option{forcefix};
}
if($option{maxannogaplen} != 50)
{
	$annoOpt .= " -maxannogaplen ".$option{maxannogaplen};
	$optCont2Ass .= " -maxannogaplen ".$option{maxannogaplen};
}
if($option{minpctid} != 0.4)
{
	$annoOpt .= " -minpctid ".$option{minpctid};
	$optCont2Ass .= " -minpctid ".$option{minpctid};
}
if($option{minpctlen} != 0.2)
{
	$annoOpt .= " -minpctlen ".$option{minpctlen};
	$optCont2Ass .= " -minpctlen ".$option{minpctlen};
}
if($option{sampname})
{
	$annoOpt .= " -sampname ".$option{sampname};
	$optCont2Ass .= " -sampname ".$option{sampname};
}
if($option{virus})
{
	$QAopt .= " -virusName ".$option{virus};
	$optCont2Ass .= " -virus ".$option{virus};
}

my $hasreads = 0;
if($option{readfa} || $option{readfq})
{
	$hasreads = 1;
}


print COMMANDLOG "perl $scriptpath"."vfat.pl -contigs $inputcontIn -ref $refIn -genelist ".$option{genelist}." -amps ".$option{amps}." -pepfolder ".$option{pepfolder}." -o $output".$optCont2Ass."\n";

### Run orientContig.pl (Order and Orient)
print "Running orientContig.pl...\n";
system("perl $scriptpath"."orientContig.pl $inputcontIn $refIn $output".$orientcontOpt);
print COMMANDLOG "perl $scriptpath"."orientContig.pl $inputcontIn $refIn $output".$orientcontOpt."\n";

### Run contigMerger.pl (merging of contigs to form 1 assembly, which can be X padded)
print "Running contigMerger.pl...\n";
print COMMANDLOG "perl $scriptpath"."contigMerger.pl $output"."_orientedContigs $refIn $output"."_merger".$contigmergerOpt."\n";
system("perl $scriptpath"."contigMerger.pl $output"."_orientedContigs $refIn $output"."_merger".$contigmergerOpt);

### Run fixFrameshifts.pl (fix Frameshifts in assembly based on reference comparison)
print "Running fixFrameshifts.pl...\n";
if($hasreads)
{
	system("perl $scriptpath"."fixFrameshifts.pl -fa $output"."_merger_assembly.fa -ref $refIn -o $output -genelist ".$option{genelist}." -qlx $output"."_merger_readAlignment.qlx".$fsOpt);
	print COMMANDLOG "perl $scriptpath"."fixFrameshifts.pl -fa $output"."_merger_assembly.fa -ref $refIn -o $output -genelist ".$option{genelist}." -qlx $output"."_merger_readAlignment.qlx".$fsOpt."\n";
	### Make final read alignment
	if($option{readfa})
	{
		if($option{readfa2})
		{
			system("perl $scriptpath"."runMosaik2.pl -fa ".$option{readfa}." -qual ".$option{readq}." -fa2 ".$option{readfa2}." -qual2 ".$option{readq2}." -ref $refIn -o $output"."_alignVsRef -qlx".$mosaikParam);
			system("perl $scriptpath"."runMosaik2.pl -fa ".$option{readfa}." -qual ".$option{readq}." -fa2 ".$option{readfa2}." -qual2 ".$option{readq2}." -ref $output"."_fixed_assembly.fa -o $output"."_alignVsAssembly -qlx".$mosaikParam);
			print COMMANDLOG "perl $scriptpath"."runMosaik2.pl -fa ".$option{readfa}." -qual ".$option{readq}." -fa2 ".$option{readfa2}." -qual2 ".$option{readq2}." -ref $refIn -o $output"."_alignVsRef -qlx".$mosaikParam."\n";
			print COMMANDLOG "perl $scriptpath"."runMosaik2.pl -fa ".$option{readfa}." -qual ".$option{readq}." -fa2 ".$option{readfa2}." -qual2 ".$option{readq2}." -ref $output"."_fixed_assembly.fa -o $output"."_alignVsAssembly -qlx".$mosaikParam."\n";
		}else{
			system("perl $scriptpath"."runMosaik2.pl -fa ".$option{readfa}." -qual ".$option{readq}." -ref $refIn -o $output"."_alignVsRef -qlx".$mosaikParam);
			system("perl $scriptpath"."runMosaik2.pl -fa ".$option{readfa}." -qual ".$option{readq}." -ref $output"."_fixed_assembly.fa -o $output"."_alignVsAssembly -qlx".$mosaikParam);
			print COMMANDLOG "perl $scriptpath"."runMosaik2.pl -fa ".$option{readfa}." -qual ".$option{readq}." -ref $refIn -o $output"."_alignVsRef -qlx".$mosaikParam."\n";
			print COMMANDLOG "perl $scriptpath"."runMosaik2.pl -fa ".$option{readfa}." -qual ".$option{readq}." -ref $output"."_fixed_assembly.fa -o $output"."_alignVsAssembly -qlx".$mosaikParam."\n";
		}
	}elsif($option{readfq})
	{
		if($option{readfq2})
		{
			system("perl $scriptpath"."runMosaik2.pl -fq ".$option{readfq}." -fq2 ".$option{readfq2}." -ref $refIn -o $output"."_alignVsRef -qlx".$mosaikParam);
			system("perl $scriptpath"."runMosaik2.pl -fq ".$option{readfq}." -fq2 ".$option{readfq2}." -ref $output"."_fixed_assembly.fa -o $output"."_alignVsAssembly -qlx".$mosaikParam);
			print COMMANDLOG "perl $scriptpath"."runMosaik2.pl -fq ".$option{readfq}." -fq2 ".$option{readfq2}." -ref $refIn -o $output"."_alignVsRef -qlx".$mosaikParam."\n";
			print COMMANDLOG "perl $scriptpath"."runMosaik2.pl -fq ".$option{readfq}." -fq2 ".$option{readfq2}." -ref $output"."_fixed_assembly.fa -o $output"."_alignVsAssembly -qlx".$mosaikParam."\n";
		}else{
			system("perl $scriptpath"."runMosaik2.pl -fq ".$option{readfq}." -ref $refIn -o $output"."_alignVsRef -qlx".$mosaikParam);
			system("perl $scriptpath"."runMosaik2.pl -fq ".$option{readfq}." -ref $output"."_fixed_assembly.fa -o $output"."_alignVsAssembly -qlx".$mosaikParam);
			print COMMANDLOG "perl $scriptpath"."runMosaik2.pl -fq ".$option{readfq}." -ref $refIn -o $output"."_alignVsRef -qlx".$mosaikParam."\n";
			print COMMANDLOG "perl $scriptpath"."runMosaik2.pl -fq ".$option{readfq}." -ref $output"."_fixed_assembly.fa -o $output"."_alignVsAssembly -qlx".$mosaikParam."\n";
		}
	}
}else{
	system("perl $scriptpath"."fixFrameshifts.pl -fa $output"."_merger_assembly.fa -ref $refIn -o $output -genelist ".$option{genelist}.$fsOpt);
	print COMMANDLOG "perl $scriptpath"."fixFrameshifts.pl -fa $output"."_merger_assembly.fa -ref $refIn -o $output -genelist ".$option{genelist}.$fsOpt."\n";
}

if($option{noannot}){exit();}

### Run annotate.pl (annotation of genes based on reference using genewise)
print "Running annotate.pl...\n";
system("perl $scriptpath"."annotate.pl -fa $output"."_fixed_assembly.fa -ref $refIn -o $output"."_annotation -pepfolder ".$option{pepfolder}." -genelist ".$option{genelist}.$annoOpt);
print COMMANDLOG "perl $scriptpath"."annotate.pl -fa $output"."_fixed_assembly.fa -ref $refIn -o $output"."_annotation -pepfolder ".$option{pepfolder}." -genelist ".$option{genelist}.$annoOpt."\n";

### Run QA_stats.pl (quality assessment statistics)
if($hasreads)
{
	if($option{readfa})
	{
		if($option{readfa2})
		{
			system("perl $scriptpath"."QA_stats.pl -ref $refIn -genelist ".$option{genelist}." -amps ".$option{amps}." -assem $output"."_fixed_assembly.fa -readfa ".$option{readfa}." -readfa2 ".$option{readfa2}." -qlxref $output"."_alignVsRef.qlx -qlxassem $output"."_alignVsAssembly.qlx -mergelist $output"."_merger_mergingList.txt -annot $output"."_annotation_summary.txt -mergeR $output"."_merger_contigsMap.R -o $output"."_QA".$QAopt);
			print COMMANDLOG "perl $scriptpath"."QA_stats.pl -ref $refIn -genelist ".$option{genelist}." -amps ".$option{amps}." -assem $output"."_fixed_assembly.fa -readfa ".$option{readfa}." -readfa2 ".$option{readfa2}." -qlxref $output"."_alignVsRef.qlx -qlxassem $output"."_alignVsAssembly.qlx -mergelist $output"."_merger_mergingList.txt -annot $output"."_annotation_summary.txt -mergeR $output"."_merger_contigsMap.R -o $output"."_QA".$QAopt."\n";
		}else{
			system("perl $scriptpath"."QA_stats.pl -ref $refIn -genelist ".$option{genelist}." -amps ".$option{amps}." -assem $output"."_fixed_assembly.fa -readfa ".$option{readfa}." -qlxref $output"."_alignVsRef.qlx -qlxassem $output"."_alignVsAssembly.qlx -mergelist $output"."_merger_mergingList.txt -annot $output"."_annotation_summary.txt -mergeR $output"."_merger_contigsMap.R -o $output"."_QA".$QAopt);
			print COMMANDLOG "perl $scriptpath"."QA_stats.pl -ref $refIn -genelist ".$option{genelist}." -amps ".$option{amps}." -assem $output"."_fixed_assembly.fa -readfa ".$option{readfa}." -qlxref $output"."_alignVsRef.qlx -qlxassem $output"."_alignVsAssembly.qlx -mergelist $output"."_merger_mergingList.txt -annot $output"."_annotation_summary.txt -mergeR $output"."_merger_contigsMap.R -o $output"."_QA".$QAopt."\n";
		}
	}elsif($option{readfq})
	{
		if($option{readfq2})
		{
			system("perl $scriptpath"."QA_stats.pl -ref $refIn -genelist ".$option{genelist}." -amps ".$option{amps}." -assem $output"."_fixed_assembly.fa -readfq ".$option{readfq}." -readfq2 ".$option{readfq2}." -qlxref $output"."_alignVsRef.qlx -qlxassem $output"."_alignVsAssembly.qlx -mergelist $output"."_merger_mergingList.txt -annot $output"."_annotation_summary.txt -mergeR $output"."_merger_contigsMap.R -o $output"."_QA".$QAopt);
			print COMMANDLOG "perl $scriptpath"."QA_stats.pl -ref $refIn -genelist ".$option{genelist}." -amps ".$option{amps}." -assem $output"."_fixed_assembly.fa -readfq ".$option{readfq}." -readfq2 ".$option{readfq2}." -qlxref $output"."_alignVsRef.qlx -qlxassem $output"."_alignVsAssembly.qlx -mergelist $output"."_merger_mergingList.txt -annot $output"."_annotation_summary.txt -mergeR $output"."_merger_contigsMap.R -o $output"."_QA".$QAopt."\n";
		}else{
			system("perl $scriptpath"."QA_stats.pl -ref $refIn -genelist ".$option{genelist}." -amps ".$option{amps}." -assem $output"."_fixed_assembly.fa -readfq ".$option{readfq}." -qlxref $output"."_alignVsRef.qlx -qlxassem $output"."_alignVsAssembly.qlx -mergelist $output"."_merger_mergingList.txt -annot $output"."_annotation_summary.txt -mergeR $output"."_merger_contigsMap.R -o $output"."_QA".$QAopt);
			print COMMANDLOG "perl $scriptpath"."QA_stats.pl -ref $refIn -genelist ".$option{genelist}." -amps ".$option{amps}." -assem $output"."_fixed_assembly.fa -readfq ".$option{readfq}." -qlxref $output"."_alignVsRef.qlx -qlxassem $output"."_alignVsAssembly.qlx -mergelist $output"."_merger_mergingList.txt -annot $output"."_annotation_summary.txt -mergeR $output"."_merger_contigsMap.R -o $output"."_QA".$QAopt."\n";
		}
	}
}else{
	system("perl $scriptpath"."QA_stats.pl -ref $refIn -genelist ".$option{genelist}." -amps ".$option{amps}." -assem $output"."_fixed_assembly.fa -mergelist $output"."_merger_mergingList.txt -annot $output"."_annotation_summary.txt -mergeR $output"."_merger_contigsMap.R -o $output"."_QA".$QAopt);
	print COMMANDLOG "perl $scriptpath"."QA_stats.pl -ref $refIn -genelist ".$option{genelist}." -amps ".$option{amps}." -assem $output"."_fixed_assembly.fa -mergelist $output"."_merger_mergingList.txt -annot $output"."_annotation_summary.txt -mergeR $output"."_merger_contigsMap.R -o $output"."_QA".$QAopt."\n";
}


### Clean up files
unless(-d "$output"."_additionalFiles")
{
	system("mkdir $output"."_additionalFiles");
}

system("mv $output"."_alignPair_fixed.afa $output"."_assemVsRef.afa");
system("mv $output"."_orientedContigs_* $output"."_additionalFiles");
system("mv $output"."_orientedContigs.fa $output"."_additionalFiles");
system("mv $output"."_merger_vsRef.mfa $output"."_additionalFiles");
system("mv $output"."_merger_vsRef.afa $output"."_additionalFiles");
system("mv $output"."_merger_assembly.fa $output"."_additionalFiles");
system("mv $output"."_merger_contigsMap.R $output"."_additionalFiles");
system("mv $output"."_alignPair*.mfa $output"."_additionalFiles");
system("mv $output"."_alignPair.afa $output"."_additionalFiles");
system("mv $output"."_annotation_aligned* $output"."_additionalFiles");
system("mv $output"."_annotation_details* $output"."_additionalFiles");
#if(-d "$output"."_additionalFiles/$output"."_annotation")
#{
#	system("rm -r $output"."_additionalFiles/$output"."_annotation");
#}
system("mv $output"."_annotation $output"."_additionalFiles");

if($option{readfa} || $option{readfq})
{
	system("rm $output"."_*aligned.qlx");
	system("rm $output"."*_unmappedIDs.txt");
	system("rm $output"."*multiple.bam");
	system("mv $output"."_merger_readAlignment.qlx $output"."_additionalFiles");
	system("mv $output"."_QA_coverageVsAssembly.R $output"."_additionalFiles");
	system("mv $output"."_QA_coverageVsRef.R $output"."_additionalFiles");
	system("mv $output"."_alignVsAssembly.qlx $output"."_additionalFiles");
	system("mv $output"."_alignVsAssembly.sam $output"."_additionalFiles");
	system("mv $output"."_alignVsRef.qlx $output"."_additionalFiles");
	system("mv $output"."_alignVsRef.sam $output"."_additionalFiles");
}

#unless($option{details})
#{
#	system("rm -r $output"."_additionalFiles");
#}

