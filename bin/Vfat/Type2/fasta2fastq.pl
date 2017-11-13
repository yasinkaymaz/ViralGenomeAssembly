#!/share/pkg/perl/5.18.1/bin/ perl -w

# Copyright © 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

use warnings;
use strict;
use File::Basename;

my $inFasta = $ARGV[0]  || die ("Usage : perl fasta2fastq.pl inputReads\nThe script will assume that the files are named inputReads.fa/.fasta/.fna, inputReads.qual and inputReads.fastq\n");
my $baseName = basename($inFasta, qw/.fasta .fna .fa/);
my $inQual = $baseName . ".qual";
my $outFastq = $baseName . ".fastq";

my %seqs;

$/ = ">";

open (FASTA, "<$inFasta");
my $junk = (<FASTA>);

while (my $frecord = <FASTA>) {
	chomp $frecord;
	my ($fdef, @seqLines) = split /\n/, $frecord;
	my $seq = join '', @seqLines;
	$seqs{$fdef} = $seq;
}

close FASTA;

open (QUAL, "<$inQual");
$junk = <QUAL>;
open (FASTQ, ">$outFastq");

while (my $qrecord = <QUAL>) {
	chomp $qrecord;
	my ($qdef, @qualLines) = split /\n/, $qrecord;
	my $qualString = join ' ', @qualLines;
	my @quals = split / /, $qualString;
	print FASTQ "@","$qdef\n";
	print FASTQ "$seqs{$qdef}\n";
	print FASTQ "+\n";
	foreach my $qual (@quals) {
		print FASTQ chr($qual + 33);
	}
	print FASTQ "\n";
}

close QUAL;
close FASTQ;
