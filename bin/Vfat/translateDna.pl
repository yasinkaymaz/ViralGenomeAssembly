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
	start   	=> 0,
	stop		=> 0,
	genelist	=> '',
);
GetOptions(
  "start=i"		=> \$option{start},
  "stop=i"		=> \$option{stop},
  "genelist=s"	=> \$option{genelist},
) || die("Problem processing command-line options: $!\n");


my $input = shift || die("Usage : perl translateDna.pl input.fa output.fa\nOptions:\n-start : start position for the translation\n-stop : stop position for the translation\n-genelist : supply a genelist. All genes will be translated separately and put in an output folder instead of an output file\n");
my $output = shift || die("Usage : perl translateDna.pl input.fa output.fa\nOptions:\n-start : start position for the translation\n-stop : stop position for the translation\n-genelist : supply a genelist. All genes will be translated separately and put in an output folder instead of an output file\n");

open(INPUT, $input) || die("Unable to open $input\n");

my %geneticcode = (
'TTT'=>'F','TTC'=>'F',
'TTA'=>'L','TTG'=>'L','CTT'=>'L','CTC'=>'L','CTA'=>'L','CTG'=>'L', 'CTN'=>'L',
'ATT'=>'I','ATC'=>'I','ATA'=>'I',
'ATG'=>'M',	
'TCT'=>'S','TCC'=>'S','TCA'=>'S','TCG'=>'S','AGT'=>'S','AGC'=>'S', 'TCN'=>'S',
'CCT'=>'P','CCC'=>'P','CCA'=>'P','CCG'=>'P', 'CCN'=>'P',
'ACT'=>'T','ACC'=>'T','ACA'=>'T','ACG'=>'T', 'ACN'=>'T',
'GCT'=>'A','GCC'=>'A','GCA'=>'A','GCG'=>'A', 'GCN'=>'A',
'TAT'=>'Y','TAC'=>'Y',
'TAA'=>'*','TAG'=>'*','TGA'=>'*',
'CAT'=>'H','CAC'=>'H',
'CAA'=>'Q','CAG'=>'Q',
'AAT'=>'N','AAC'=>'N',
'AAA'=>'K','AAG'=>'K',
'GAT'=>'D','GAC'=>'D',
'GAA'=>'E','GAG'=>'E',
'TGT'=>'C','TGC'=>'C',
'TGG'=>'W',
'AGA'=>'R',	'AGG'=>'R', 'CGT'=>'R',	'CGC'=>'R','CGA'=>'R','CGG'=>'R', 'CGN'=>'R',
'GGT'=>'G','GGC'=>'G','GGA'=>'G','GGG'=>'G', 'GGN'=>'G',
'GTT'=>'V','GTC'=>'V','GTA'=>'V','GTG'=>'V', 'GTN'=>'V');

my @seqs;
my $seqflag = -1;
my %genes;

if($option{genelist})
{
	open(GENES, $option{genelist});
	while(my $line = <GENES>)
	{
		chomp $line;
		if($line =~ /(.+?)\t(.+?)\t(.+)/)
		{
			$genes{$1}{start} = $2;
			$genes{$1}{stop} = $3;
		}
	}
}

while(my $line = <INPUT>)
{
	chomp $line;
	
	if($line =~ />(.*)/)
	{
		$seqflag++;
		$seqs[$seqflag]{'name'} = $1;
		next;
	}else{
		$seqs[$seqflag]{'seq'} .= $line;
	}
}

foreach my $curseq (@seqs)
{
	if($option{genelist})
	{
		system("mkdir $output");
#		open(MAINOUTPUT, ">$output/allPeptides.fa");
		foreach my $gene (sort keys %genes)
		{
			open(CURGENE, ">$output/$gene"."_pep.fa");
			my $geneseq = substr($$curseq{'seq'}, ($genes{$gene}{start}-1), ($genes{$gene}{stop}-$genes{$gene}{start}+1));
			my $aaseq = translateDnaString($geneseq);
#			print MAINOUTPUT ">".$gene."\n$aaseq\n";
			print CURGENE ">".$gene."\n$aaseq\n";
			close CURGENE;
		}
	}else{
		open(OUTPUT, ">$output");
		if($option{start})
		{
			$$curseq{'seq'} = substr($$curseq{'seq'}, ($option{start}-1), ($option{stop}-$option{start}+1));
		}
		my $aaseq = translateDnaString($$curseq{'seq'});
		print OUTPUT ">".$$curseq{'name'}."\n$aaseq\n";
	}
}




sub translateDnaString
{
	my $dnastring = shift;
	my $aastring = "";

	while(length($dnastring) > 2)
	{
		my $curcod = substr($dnastring, 0, 3);
		my $curaa;
		if($geneticcode{$curcod}){$curaa = $geneticcode{$curcod};}
		elsif($curcod =~ /\-/){$curaa = "-";}
		else{$curaa = "X";}
		
		$aastring .= $curaa;
		if(length($dnastring) == 3){$dnastring = "";}
		else{$dnastring = substr($dnastring, 3);}
	}
	return $aastring;
}
