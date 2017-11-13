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
	genelist		=> '',
	ref				=> '',
	assem			=> '',
	qlxref			=> '',
	qlxassem		=> '',
	readfa			=> '',
	readfa2			=> '',
	readfq			=> '',
	readfq2			=> '',
	mergelist		=> '',
	amps			=> '',
	o				=> '',
	annot			=> '',
	mergeR			=> '',
	virusName		=> '',
	h				=> '',
);

GetOptions(
	"genelist=s"		=> \$option{genelist},
	"ref=s"				=> \$option{ref},
	"assem=s"			=> \$option{assem},
	"qlxref=s"			=> \$option{qlxref},
	"qlxassem=s"		=> \$option{qlxassem},
	"readfa=s"			=> \$option{readfa},
	"readfa2=s"			=> \$option{readfa2},
	"readfq=s"			=> \$option{readfq},
	"readfq2=s"			=> \$option{readfq2},
	"mergelist=s"		=> \$option{mergelist},
	"amps=s"			=> \$option{amps},
	"o=s"				=> \$option{o},
	"annot=s"			=> \$option{annot},
	"mergeR=s"			=> \$option{mergeR},
	"virusName=s"		=> \$option{virusName},
	"h"					=> \$option{h},
);

if($option{h})
{
	print "Usage : perl QA_stats.pl -ref reference.fasta -genelist reference_genelist.txt -amps reference_amplicons.txt -assem assembly.fa -mergelist mergingList.txt -annot annotation_summary.txt -mergeR merger_contigsMap.R -o output\n";
	print "QA_stats.pl requires R to run\n";
	print "Options:\n";
	print "If reads are supplied:\n";
	print "-readfa\tReads file in fasta format\n";
	print "-readfa2\tReads file in fasta format (2nd mate if paired)\n";
	print "-readfq\tReads file in fastq format\n";
	print "-readfq2\tReads file in fastq format (2nd mate if paired)\n";
	print "-qlxref\tRead alignment in qlx format against the reference\n";
	print "-qlxassem\tRead alignment in qlx format against the assembly\n";
	print "Other options\n";
	print "-virusName\tName of the virus (for outputs)\n";
	exit();
}

my $output = $option{o};
my $mincall = 1;

my $Rpath = "/share/pkg/R/3.0.1/bin/";

unless($option{virusName})
{
	my @refinput = split(/\//, $option{ref});
	my $reffile = $refinput[scalar(@refinput)-1];
	my @refpart = split(/\_/, $reffile);
	$option{virusName} = $refpart[0];
}

open(DETAILED, ">$output"."_StatsDetailed.txt");

printFilesHeader();

my $hasreads = 0;
if($option{readfa} || $option{readfq}){$hasreads = 1;}


my $totalReads = 0;
my $avgCovTarget = 0;
my $stdCovTarget = 0;
my $avgCovAss = 0;
my $stdCovAss = 0;
my @allGeneCov;
my @allAmpCov;
my $nonDomBase = 0;
my $nonDomDel = 0;
my $nonDomIns = 0;

#### Get Nb of Reads
if($option{readfa})
{
	open(READSFA, $option{readfa});
	while(my $line = <READSFA>)
	{
		if($line =~ />/){$totalReads++;}
	}
	close READSFA;
	if($option{readfa2})
	{
		open(READSFA, $option{readfa2});
		while(my $line = <READSFA>)
		{
			if($line =~ />/){$totalReads++;}
		}
		close READSFA;
	}
}elsif($option{readfq})
{
	open(READSFQ, $option{readfq});
	while(my $line = <READSFQ>)
	{
		chomp $line;
		if($line eq '+'){$totalReads++;}
	}
	close READSFQ;
	if($option{readfq2})
	{
		open(READSFQ, $option{readfq2});
		while(my $line = <READSFQ>)
		{
			chomp $line;
			if($line eq '+'){$totalReads++;}
		}
		close READSFQ;
	}
}

### Read in Reference
my $refseq;
open(REFFILE, $option{ref});
while(my $line = <REFFILE>)
{
	chomp $line;
	if($line =~ />.+/)
	{
		next;
	}else{
		$refseq .= $line;
	}
}
close REFFILE;

### Read in Assembly
my $assseq;
open(ASSFILE, $option{assem});
while(my $line = <ASSFILE>)
{
	chomp $line;
	if($line =~ />.+/)
	{
		next;
	}else{
		$assseq .= $line;
	}
}
close ASSFILE;

#### Read Genelist
open(GENELIST, $option{genelist});
my @genes;
my %target;
$target{start} = 0;
my $geneflag = 0;
while(my $line = <GENELIST>)
{
	chomp $line;
	if($line =~ /(.+?)\t(.+?)\t(.+)/)
	{
		if($target{start} == 0){$target{start} = $2;}
		$target{stop} = $3;
		$genes[$geneflag]{name} = $1;
		$genes[$geneflag]{start} = $2;
		$genes[$geneflag]{stop} = $3;
		$geneflag++;
	}
}
close GENELIST;

### Get amplicons positions
open(AMPS, $option{amps});
my @amps;
my $ampsflag = 0;
my %ampcov;
while(my $line = <AMPS>)
{
	chomp $line;
	if($line =~ /(.+?)\t(.+?)\t(.+)/)
	{
		$amps[$ampsflag]{name} = $1;
		$amps[$ampsflag]{start} = $2;
		$amps[$ampsflag]{stop} = $3;
		for(my $pos = $amps[$ampsflag]{start}; $pos <= $amps[$ampsflag]{stop}; $pos++)
		{
			$ampcov{$pos}++;
		}
		$ampsflag++;
	}
}
my @uniqueamps;
$ampsflag = 0;
foreach my $amp(@amps)
{
	my $start = 0;
	my $stop = 0;
	for(my $pos = $$amp{start}; $pos <= $$amp{stop}; $pos++)
	{
		unless($start)
		{
			if($ampcov{$pos} == 1){$start = $pos;}
		}
		if($ampcov{$pos} == 1){$stop = $pos;}elsif($start > 0){last;}
	}
	$uniqueamps[$ampsflag]{name} = $$amp{name};
	$uniqueamps[$ampsflag]{start} = $start;
	$uniqueamps[$ampsflag]{stop} = $stop;
#	print "Unique Amp : ".$$amp{name}."\t".$start."\t".$stop."\n";
	$ampsflag++;
}


close AMPS;

### Read Mergelist output
open(MERGELIST, $option{mergelist});
my %poscovered;
my $nbcontigs = 0;
my $nbNstr = 0;
while(my $line = <MERGELIST>)
{
	chomp $line;
	if($line =~ /(.+?)\t(.+?)\t(.+)/)
	{
		my $start = $1;
		my $stop = $2;
		my $contdetail = $3;
		unless($contdetail =~ /^[0-9]+N$/)
		{
			for(my $pos = $start; $pos <= $stop; $pos++)
			{
				$poscovered{$pos} = 1;
			}
		}else{
			$nbNstr++;
		}
	}elsif($line =~ /NB CONTIG USED:	(.+)/)
	{
		$nbcontigs = $1;
	}
}
close MERGELIST;

my $posmatch = 0;
my $nbpos = $target{stop} - $target{start} + 1;
for(my $pos = $target{start}; $pos <= $target{stop}; $pos++)
{
	if($poscovered{$pos}){$posmatch++;}
}
my $pctrefcov = sprintf("%.3f",($posmatch/$nbpos*100));

### Read annotation output
open(ANNOT, $option{annot});
my %annotStats;
$annotStats{fullstartstop} = 0;
$annotStats{fullfs} = 0;
$annotStats{fullNstr} = 0;
$annotStats{nbfull} = 0;
$annotStats{nbpart} = 0;
$annotStats{nbmiss} = 0;

my $firstline = 1;
while(my $line = <ANNOT>)
{
	chomp $line;
	if($firstline){$firstline = 0; next;}
	unless($line =~ /\t/){next;}
	my @linedata = split(/\t/, $line);
	if($linedata[2] eq 'Full'){
		$annotStats{nbfull}++;
		if($linedata[6] == 0){
			$annotStats{fullstartstop}++;
		}
		if($linedata[7] == 0){
			$annotStats{fullfs}++;
		}
		if($linedata[8] == 0){
			$annotStats{fullNstr}++;
		}
	}elsif($linedata[2] eq 'Partial'){
		$annotStats{nbpart}++;
	}elsif($linedata[2] eq 'Missing'){
		$annotStats{nbmiss}++;
	}
}
close ANNOT;

### Print general info
print DETAILED "Reference data\n\n";
print DETAILED "Virus : ".$option{virusName}."\n";
print DETAILED "Length reference : ".length($refseq)."\n";
print DETAILED "Nb Amplicons : ".scalar(@amps)."\n";
print DETAILED "Nb Genes : ".($annotStats{nbfull} + $annotStats{nbpart} + $annotStats{nbmiss})."\n";
print DETAILED "Target region : ".$target{start}."-".$target{stop}."\n";	
print DETAILED "\n---------------------------\n\n";

### Print assembly info
print DETAILED "Assembly QC General Stats\n\n";
print DETAILED "Length Assembly : ".length($assseq)."\n";
print DETAILED "% Reference Covered : $pctrefcov"."%\n";
print DETAILED "Nb Contigs Used in Assembly : $nbcontigs\n";
print DETAILED "Nb N strings inserted to merge Assembly: $nbNstr\n";
print DETAILED "Nb Genes Fully Covered : ".$annotStats{nbfull}."\n";
print DETAILED "Nb Genes Partially Covered : ".$annotStats{nbpart}."\n";
print DETAILED "Nb Genes Missing : ".$annotStats{nbmiss}."\n";
print DETAILED "Nb Full Genes with Frameshift : ".$annotStats{fullfs}."\n";
print DETAILED "Nb Full Genes missing Start/Stop codon : ".$annotStats{fullstartstop}."\n";
print DETAILED "Nb Full Genes with N-string : ".$annotStats{fullNstr}."\n";
print DETAILED "\n---------------------------\n\n";


unless($hasreads){
	### Print some output ###
	print DETAILED "*** No Reads were used for Coverage Statistics ***\n";
	writeExcelOutput();
	exit();
}

print DETAILED "Assembly QC Coverage Stats\n\n";
print DETAILED "Nb Reads in input : $totalReads\n";


my %ntfreq;
my %nbReadAln;

my %minorvar;
my @coverages;
my %avgCoverage;

my %pctCovered;
$pctCovered{'1X'} = 0;
$pctCovered{'10X'} = 0;
$pctCovered{'50X'} = 0;
$pctCovered{'200X'} = 0;

readQlx($option{qlxref}, 'ref');
#print "% Reads Aligned Ref : ".($nbReadAln{'ref'}/$totalReads*100)."\n";
my $pctReadsToRef = sprintf("%.3f",($nbReadAln{'ref'}/$totalReads*100));
print DETAILED "% Reads Aligned to Reference : $pctReadsToRef"."%\n";
readQlx($option{qlxassem}, 'assem');
my $pctReadsToAss = sprintf("%.3f",($nbReadAln{'assem'}/$totalReads*100));
print DETAILED "% Reads Aligned to Assembly : $pctReadsToAss"."%\n";

printNtOutput('ref');
printNtOutput('assem');

plotVsRef();
plotVsAss();

writeExcelOutput();

sub writeExcelOutput
{
	open(XLS, ">$output"."_StatsSummary.xls");
	print XLS "Project\tVirus\t%Reference Covered by Assembly\t#Contigs Used for Merging\t#N-strings inserted for Merging\t#Genes Full\t#Genes Partial\t#Genes Missing\t#Full Genes w/ Frameshift\t#Full Genes Missing Start/Stop codon\t#Full Genes w/ N-string";
	
	if($hasreads)
	{
		print XLS "\t#Reads\t%Reads Aligned To Ref\t%Reads Aligned To Assembly\tAvg Coverage (vs Ref, target region only)\tStdDev on Coverage";
		foreach my $amp(@amps)
		{
			print XLS "\tAvg Coverage ".$$amp{name};
		}
		foreach my $gene(@genes)
		{
			print XLS "\tAvg Coverage ".$$gene{name};
		}
		print XLS "\tAvg Coverage vs Assembly\tStdDev on Coverage vs Assembly\t%Loci with non-dominant base\t%Loci with non-dominant deletion\t%Loci with non-dominant insertion";
	}
	print XLS "\n";

	print XLS $output."\t".$option{virusName}."\t".$pctrefcov."\t".$nbcontigs."\t".$nbNstr."\t".$annotStats{nbfull}."\t".$annotStats{nbpart}."\t".$annotStats{nbmiss}."\t".$annotStats{fullfs}."\t".$annotStats{fullstartstop}."\t".$annotStats{fullNstr};
	if($hasreads)
	{
		print XLS "\t".$totalReads."\t".$pctReadsToRef."\t".$pctReadsToAss."\t".$avgCovTarget."\t".$stdCovTarget;
		my $flag = 0;
		foreach my $amp(@amps)
		{
			print XLS "\t".$allAmpCov[$flag];
			$flag++;
		}
		$flag = 0;
		foreach my $gene(@genes)
		{
			print XLS "\t".$allGeneCov[$flag];
			$flag++;
		}
		print XLS "\t$avgCovAss\t$stdCovAss\t$nonDomBase\t$nonDomDel\t$nonDomIns";
	}
	print XLS "\n";	
}

sub printFilesHeader
{
	print DETAILED "Files used:\n";
	print DETAILED "Reference File : ".$option{ref}."\n";
	print DETAILED "Assembly File : ".$option{assem}."\n";
	print DETAILED "Read Alignment vs Reference : ".$option{qlxref}."\n";
	print DETAILED "Read Alignment vs Assembly : ".$option{qlxassem}."\n";
	print DETAILED "Reference Genelist : ".$option{genelist}."\n";
	print DETAILED "Reference Amplicons : ".$option{amps}."\n";
	print DETAILED "Contig merging list (output of contigMerger.pl) : ".$option{mergelist}."\n";
	print DETAILED "Contig merging R file  (output of contigMerger.pl): ".$option{mergeR}."\n";
	print DETAILED "Annotation summary  (output of annotate.pl): ".$option{annot}."\n";
	if($option{readfa})
	{
		if($option{readfa2})
		{
			print DETAILED "Reads file : ".$option{readfa}." and ".$option{readfa2}."\n";
		}else{
			print DETAILED "Reads file : ".$option{readfa}."\n";
		}
	}elsif($option{readfq})
	{
		if($option{readfq2})
		{
			print DETAILED "Reads file : ".$option{readfq}." and ".$option{readfq2}."\n";
		}else{
			print DETAILED "Reads file : ".$option{readfq}."\n";
		}
	}else{
		print DETAILED "\n*** No read file supplied ***\n";
	}
	print DETAILED "\n---------------------------\n\n";
}

sub readQlx
{
	my $qlxfile = shift;
	my $reftype = shift;

	my $refstart;
	my $refend;
	my $readid;
	my $readstring;
	my $refstring;
	my $qualstring;
	my $state = 0;

	open(QLXFILE, $qlxfile);
	
	while (my $line = <QLXFILE>)
	{
		chomp $line;
		
		if($line =~ />Read/)
		{
			$nbReadAln{$reftype}++;
			my @query = split(/\s/, $line);
			$readid = $query[1];
			$refstart = $query[7];
			$refend = $query[8];
			$readstring = "";
			$refstring = "";
			$qualstring = "";
			$state = 1;
		}elsif($line =~ /[ATGC]/ && $state == 1)
		{
			$refstring .= $line;
			$state = 2;
		}elsif($line =~ /[ATGC]/ && $state == 2)
		{
			$readstring .= $line;
			$state = 3;
			evalReadNt($readstring, $refstring, $refstart, $reftype);
		}
	}
}

sub evalReadNt
{
	my $readstr = shift;
	my $refstr = shift;
	my $refstart = shift;
	my $reftype = shift;

	#READ NTS
	my $ntpos = $refstart;
	my $curinsert = '';
	my $goodinsert = 1;
	my $curdel = '';
	my $delstart = 0;
	my $insertstart = 0;
	
	for(my $currefpos = 0; $currefpos < length($refstr); $currefpos++)
	{
		my $curnt = substr($readstr, $currefpos, 1);
		my $curref = substr($refstr, $currefpos, 1);
		
		if($curref eq "-")
		{
			$curinsert .= $curnt;
		}elsif($curnt eq "-")
		{
			unless($curdel){$delstart = $ntpos;}
			$curdel .= '-';
		}else{
			if($curinsert)
			{
				$ntfreq{$reftype}{$ntpos - 1}{I}{$curinsert}++;
				$ntfreq{$reftype}{$ntpos - 1}{I}{count}++;
			}
			$curinsert = '';

			if($curdel)
			{
				my $del = "D".length($curdel);
				for(my $delpos = $delstart; $delpos < $ntpos; $delpos++)
				{
					if($delpos == $delstart)
					{
						$ntfreq{$reftype}{$delpos}{'-'}{count}++;
					}else{
						$ntfreq{$reftype}{$delpos}{'-'}{count}++;
					}
					$ntfreq{$reftype}{$delpos}{coverage}{count}++;
				}
				$curdel = '';
				$delstart = 0;
			}
			$ntfreq{$reftype}{$ntpos}{$curnt}{count}++;
			$ntfreq{$reftype}{$ntpos}{coverage}{count}++;
		}
		if($curref ne "-")
		{
			$ntpos++;
		}
	}
}

sub plotVsRef
{
	my $contigCol = 'blue';
	my $ampCol = 'black';
	my $coverageCol = 'red';
	my $geneCol = 'green';
	
	### Read in Contigs from Merger
	open(MERGER, $option{mergeR});
	my %contiglines;
	while(my $line = <MERGER>)
	{
		chomp $line;
		if($line =~ /segments\((.+?)\,\-(.+?)\,(.+?)\,.*/)
		{
			$contiglines{$2}{$1} = $3;
		}
	}
	
	### Print R File	
	open(ROUTPUT, ">$output"."_coverageVsRef.R");
	print ROUTPUT "x=c(";

	#print X axis values
	for(my $i = 1; $i < length($refseq); $i++)
	{
		print ROUTPUT $i.",";
	}
	print ROUTPUT length($refseq).")\n";
	

	# Fill coverage matrix
	my $maxycov = 0;
	print ROUTPUT "y=c(".$ntfreq{'ref'}{1}{coverage}{count};
	for(my $curpos = 2; $curpos <= length($refseq); $curpos++)
	{
		print ROUTPUT " ,".$ntfreq{'ref'}{$curpos}{coverage}{count};
		if($ntfreq{'ref'}{$curpos}{coverage}{count} > $maxycov){$maxycov = $ntfreq{'ref'}{$curpos}{coverage}{count};}
	}
	print ROUTPUT ")\n";
	
	my $contigSpace = int(0.02*$maxycov);
	if($contigSpace < 5){$contigSpace = 5;}
	my $ampSpace = int(0.01 * $maxycov);
	if($ampSpace < 3){$ampSpace = 3};
	my $nbAmpGeneLines = scalar(@amps) + scalar(@genes);
	
	my $nbCont = scalar(keys %contiglines);
	
	print ROUTPUT "plot(x,y,col=\"$coverageCol\", xlab=\"Position\",ylab=\"Coverage\",type=\"l\",lwd=2,ann=T,cex.lab=0.8, ylim=c(".($nbCont*(-$contigSpace)).",max(y)+".($nbAmpGeneLines*$ampSpace)."), main=\"$output Coverage vs Reference\")\n";
	
	### Print Contig Segments
	foreach my $contig (sort {$a <=> $b} keys %contiglines)
	{
		foreach my $startpos (sort {$a <=> $b} keys %{$contiglines{$contig}})
		{
			print ROUTPUT "segments(".$startpos.", -".($contig * $contigSpace).", ".$contiglines{$contig}{$startpos}.", -".($contig * $contigSpace).", col=\"$contigCol\", lwd=2)\n";
		}
	}

	### Print Amplicon Segments
	my $ampno = 0;
	foreach my $curamp(@amps)
	{
		print ROUTPUT "segments(".$$curamp{start}.",max(y)+".((scalar(@genes) + scalar(@amps) - $ampno) * $ampSpace).",".$$curamp{stop}.",max(y)+".((scalar(@genes) + scalar(@amps) - $ampno) * $ampSpace).",col=\"$ampCol\",lwd=2)\n";
		$ampno++;
	}

	### Print Gene Segments
	foreach my $curgene(@genes)
	{
		print ROUTPUT "segments(".$$curgene{start}.",max(y)+".((scalar(@genes) + scalar(@amps) - $ampno) * $ampSpace).",".$$curgene{stop}.",max(y)+".((scalar(@genes) + scalar(@amps) - $ampno) * $ampSpace).",col=\"$geneCol\",lwd=2)\n";
		$ampno++;
	}

	print ROUTPUT "pdf(\"$output"."_coverageVsRef.pdf\")\n";

	system($Rpath."R < $output"."_coverageVsRef.R --no-save --slave -q");
	system("mv Rplots.pdf $output"."_coverageVsRef.pdf");

	
	
}

sub plotVsAss
{
	my $coverageCol = 'red';

	### Print R File
	open(ROUTPUT, ">$output"."_coverageVsAssembly.R");
	print ROUTPUT "x=c(";

	#print X axis values
	for(my $i = 1; $i < length($assseq); $i++)
	{
		print ROUTPUT $i.",";
	}
	print ROUTPUT length($assseq).")\n";
	

	# Fill coverage matrix
	my $maxycov = 0;
	print ROUTPUT "y=c(".$ntfreq{'assem'}{1}{coverage}{count};
	for(my $curpos = 2; $curpos <= length($assseq); $curpos++)
	{
		print ROUTPUT " ,".$ntfreq{'assem'}{$curpos}{coverage}{count};
		if($ntfreq{'assem'}{$curpos}{coverage}{count} > $maxycov){$maxycov = $ntfreq{'assem'}{$curpos}{coverage}{count};}
	}
	print ROUTPUT ")\n";
	
	print ROUTPUT "plot(x,y,col=\"$coverageCol\", xlab=\"Position\",ylab=\"Coverage\",type=\"l\",lwd=2,ann=T,cex.lab=0.8, ylim=c(0,max(y)), main=\"$output Coverage vs Assembly\")\n";
	
	print ROUTPUT "pdf(\"$output"."_coverageVsAssembly.pdf\")\n";
		
	system($Rpath."R < $output"."_coverageVsAssembly.R --no-save --slave -q");
	system("mv Rplots.pdf $output"."_coverageVsAssembly.pdf");
}


sub printNtOutput
{
	my $reftype = shift;
	undef @coverages;
	open(NTOUTPUT, ">$output"."_ntfreq_".$reftype.".txt") || die ("can't open Nt output file $output");

	print NTOUTPUT ">Nucleotide Frequency\nPos\tConsensusNt\tCoverage\tFreqA\tFreqT\tFreqG\tFreqC\tFreqDel\tFreqInsertion\tInsertions(Count)\n";
	my $mainseq = $refseq;
	if($reftype eq 'assem'){$mainseq = $assseq;}
	
	for(my $ntpos = 1; $ntpos <= length($mainseq); $ntpos++)
	{
		print NTOUTPUT $ntpos."\t".substr($mainseq, $ntpos-1, 1)."\t";
		if($ntfreq{$reftype}{$ntpos}{coverage}{count})
		{
			if($reftype eq 'ref')
			{
				if($ntpos >= $target{start} && $ntpos <= $target{stop}){push(@coverages, $ntfreq{$reftype}{$ntpos}{coverage}{count});}
			}elsif($reftype eq 'assem')
			{
				push(@coverages, $ntfreq{$reftype}{$ntpos}{coverage}{count});
			}

			print NTOUTPUT $ntfreq{$reftype}{$ntpos}{coverage}{count}."\t";
			
			
			if($reftype eq 'ref')
			{
				if($ntpos >= $target{start} && $ntpos <= $target{stop})
				{
					if($ntfreq{$reftype}{$ntpos}{coverage}{count} >= 200)
					{
						$pctCovered{'200X'}++;
						$pctCovered{'50X'}++;
						$pctCovered{'10X'}++;
						$pctCovered{'1X'}++;
					}elsif($ntfreq{$reftype}{$ntpos}{coverage}{count} >= 50)
					{
						$pctCovered{'50X'}++;
						$pctCovered{'10X'}++;
						$pctCovered{'1X'}++;
					}elsif($ntfreq{$reftype}{$ntpos}{coverage}{count} >= 10)
					{
						$pctCovered{'10X'}++;
						$pctCovered{'1X'}++;
					}elsif($ntfreq{$reftype}{$ntpos}{coverage}{count} >= 1)
					{
						$pctCovered{'1X'}++;
					}
				}
			}
			
			my $maxres = '';
			my $maxcount = 0;
			if($ntfreq{$reftype}{$ntpos}{A}{count})
			{
				print NTOUTPUT $ntfreq{$reftype}{$ntpos}{A}{count}/$ntfreq{$reftype}{$ntpos}{coverage}{count}."\t";
				if($ntfreq{$reftype}{$ntpos}{A}{count} > $maxcount){
					$maxcount = $ntfreq{$reftype}{$ntpos}{A}{count};
					$maxres = 'A';
				}
			}else{print NTOUTPUT "0\t";}
			if($ntfreq{$reftype}{$ntpos}{T}{count})
			{
				print NTOUTPUT $ntfreq{$reftype}{$ntpos}{T}{count}/$ntfreq{$reftype}{$ntpos}{coverage}{count}."\t";
				if($ntfreq{$reftype}{$ntpos}{T}{count} > $maxcount){
					$maxcount = $ntfreq{$reftype}{$ntpos}{T}{count};
					$maxres = 'T';
				}
			}else{print NTOUTPUT "0\t";}
			if($ntfreq{$reftype}{$ntpos}{G}{count})
			{
				print NTOUTPUT $ntfreq{$reftype}{$ntpos}{G}{count}/$ntfreq{$reftype}{$ntpos}{coverage}{count}."\t";
				if($ntfreq{$reftype}{$ntpos}{G}{count} > $maxcount){
					$maxcount = $ntfreq{$reftype}{$ntpos}{G}{count};
					$maxres = 'G';
				}
			}else{print NTOUTPUT "0\t";}
			if($ntfreq{$reftype}{$ntpos}{C}{count})
			{
				print NTOUTPUT $ntfreq{$reftype}{$ntpos}{C}{count}/$ntfreq{$reftype}{$ntpos}{coverage}{count}."\t";
				if($ntfreq{$reftype}{$ntpos}{C}{count} > $maxcount){
					$maxcount = $ntfreq{$reftype}{$ntpos}{C}{count};
					$maxres = 'C';
				}
			}else{print NTOUTPUT "0\t";}
			if($ntfreq{$reftype}{$ntpos}{'-'}{count})
			{
				print NTOUTPUT $ntfreq{$reftype}{$ntpos}{'-'}{count}/$ntfreq{$reftype}{$ntpos}{coverage}{count}."\t";
				if($ntfreq{$reftype}{$ntpos}{'-'}{count} > $maxcount){
					$maxcount = $ntfreq{$reftype}{$ntpos}{A}{count};
					$maxres = 'D';
				}
			}
			if($ntfreq{$reftype}{$ntpos}{I}{count})
			{
				print NTOUTPUT $ntfreq{$reftype}{$ntpos}{I}{count}/$ntfreq{$reftype}{$ntpos}{coverage}{count}."\t";
				if($reftype eq 'assem')
				{
					if($ntfreq{$reftype}{$ntpos}{I}{count} > ($ntfreq{$reftype}{$ntpos}{coverage}{count}/2) && ($ntpos >= $target{start} && $ntpos <= $target{stop}))
					{
						$minorvar{$ntpos}{I} = 1;
						$minorvar{insert}{count}++;
					}
				}
				
				foreach my $insert(sort keys %{$ntfreq{$reftype}{$ntpos}{I}})
				{
					if($insert eq 'count'){next;}
					print NTOUTPUT $insert."(".$ntfreq{$reftype}{$ntpos}{I}{$insert}.") ";
				}
				print NTOUTPUT "\n";
			}else{print NTOUTPUT "0\t\n";}
			
			if($reftype eq 'assem')
			{
				if($maxres ne substr($assseq, $ntpos-1, 1) && ($ntpos >= $target{start} && $ntpos <= $target{stop})){
#					print $ntpos."\t".$maxres."\t".substr($assseq, $ntpos-1, 1)."\n";
					if($maxres eq 'D')
					{
						$minorvar{del}{count}++;
					}else{
						$minorvar{mismatch}{count}++;
					}
				}
			}
		}else{
			print NTOUTPUT "0\t0\t0\t0\t0\t0\t0\t0\t\n";
			$ntfreq{$reftype}{$ntpos}{coverage}{count} = 0;
		}
	}
	
	if($reftype eq 'ref')
	{
		print DETAILED "\n* Coverage Data vs Reference *\n";
		my $lentarget = $target{stop} - $target{start} + 1;
		$avgCovTarget = sprintf("%.3f",&average(\@coverages));
		$stdCovTarget = sprintf("%.3f",&stdev(\@coverages));
		my $pct1X = sprintf("%.3f",(($pctCovered{'1X'}/$lentarget)*100));
		my $pct10X = sprintf("%.3f",(($pctCovered{'10X'}/$lentarget)*100));
		my $pct50X = sprintf("%.3f",(($pctCovered{'50X'}/$lentarget)*100));
		my $pct200X = sprintf("%.3f",(($pctCovered{'200X'}/$lentarget)*100));
		
		print DETAILED "Average Coverage Target: $avgCovTarget\n";
		print DETAILED "Standard Deviation Coverage: $stdCovTarget\n";
		print DETAILED "%Reference Covered at 1X : $pct1X"."%\n";
		print DETAILED "%Reference Covered at 10X : $pct10X"."%\n";
		print DETAILED "%Reference Covered at 50X : $pct50X"."%\n";
		print DETAILED "%Reference Covered at 200X : $pct200X"."%\n";
		
		print DETAILED "\n* Coverage Data By Gene (vs Reference) *\n";
		foreach my $gene(@genes)
		{
			my @genecov;
			undef @genecov;
			
			for(my $pos = $$gene{start}; $pos <= $$gene{stop}; $pos++)
			{
				push(@genecov, $ntfreq{$reftype}{$pos}{coverage}{count});
			}
			my $curGeneCov = sprintf("%.3f",&average(\@genecov));
			push(@allGeneCov, $curGeneCov);
			print DETAILED $$gene{name}."\t".$curGeneCov."\n";
		}

		print DETAILED "\n* Coverage Data By Amplicon (vs Reference, non-overlapping region only) *\n";
		foreach my $amp(@uniqueamps)
		{
			my @ampcov;
			undef @ampcov;
			
			for(my $pos = $$amp{start}; $pos <= $$amp{stop}; $pos++)
			{
				push(@ampcov, $ntfreq{$reftype}{$pos}{coverage}{count});
			}
			my $curAmpCov = sprintf("%.3f",&average(\@ampcov));
			push(@allAmpCov, $curAmpCov);
			print DETAILED $$amp{name}."\t".$curAmpCov."\n";
		}
	}elsif($reftype eq 'assem')
	{
		print DETAILED "\n* Coverage Data vs Assembly *\n";
		$avgCovAss = sprintf("%.3f",&average(\@coverages));
		$stdCovAss = sprintf("%.3f",&stdev(\@coverages));
		print DETAILED "Average Coverage Assembly: $avgCovAss\n";
		print DETAILED "Standard Deviation Coverage: $stdCovAss\n";
		
		
		unless($minorvar{mismatch}{count}){$minorvar{mismatch}{count} = 0;}
		unless($minorvar{del}{count}){$minorvar{del}{count} = 0;}
		unless($minorvar{insertion}{count}){$minorvar{insertion}{count} = 0;}
		$nonDomBase = $minorvar{mismatch}{count}."/".scalar(@coverages);
		$nonDomDel = $minorvar{del}{count}."/".scalar(@coverages);
		$nonDomIns = $minorvar{insertion}{count}."/".scalar(@coverages);
		print DETAILED "% Positions with non-dominant base call :\t".$minorvar{mismatch}{count}."/".scalar(@coverages)."\n";
		print DETAILED "% Positions with non-dominant deletion call :\t".$minorvar{del}{count}."/".scalar(@coverages)."\n";
		print DETAILED "% Positions with non-dominant insertion call :\t".$minorvar{insertion}{count}."/".scalar(@coverages)."\n";

	}
}

sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}
sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}
