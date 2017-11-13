#!/share/pkg/perl/5.18.1/bin/ perl -w

# Copyright © 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

# Convert a fastq to a fasta/qual combo using BioPerl, with some Linux commands

my $input = shift || die("Usage : perl fastq2fasta.pl input.fq output\nWill create the files output.fa and output.qual\n");
my $output = shift || die("Usage : perl fastq2fasta.pl input.fq output\nWill create the files output.fa and output.qual\n");

open(INPUT, "$input");
open(FASTA, ">$output.fa");
open(QUAL, ">$output.qual");

my $count = 0;
while(my $line = <INPUT>)
{
	chomp $line;
	if($line =~ /\@(.+)/ && $count == 0)
	{
		$count = 1;
		my $head = $1;
		$head =~ s/:/_/g;
		print FASTA ">$head\n";
		print QUAL ">$head\n";
	}elsif($count == 1)
	{
		print FASTA $line."\n";
		$count = 2;
	}elsif($count == 2)
	{
		$count = 3;
	}elsif($count == 3)
	{
		for(my $i = 0; $i < length($line); $i++)
		{
			if($i == 0)
			{
#				print QUAL (ord(substr($line,$i,1) - 33));
				print QUAL (ord(substr($line,$i,1)) - 33);
			}else{
				print QUAL " ".(ord(substr($line,$i,1)) - 33);
			}
		}
		print QUAL "\n";
		$count = 0;
	}
}
