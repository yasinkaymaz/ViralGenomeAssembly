#!/usr/local/bin/perl
#
# genbank2gff -- convert a genbank flatfile to a GFF file with all
# 		labelled features.
#
# Copyright Neomorphic Software, 1998.  All Rights Reserved.
#
# Description: A genbank flatfile is read and the features fields extracted.
#
# Usage: genbank2gff [-d [dir]] < genbank-file
#
# Options:
#
#	-d [dir]	write output in separate file in the directory
#			'dir', one file per sequence, 
#			using the name of the LOCUS with an extension
#			of '.gff'.  If no directory is specified, then 
#			use the current working directory.
#			If -d is not specified, then output is sent to stdout.
#
#	-h		print this help information
#
# Example:
#
#	genbank2gff -d /tmp < file.genbank
#

use File::Basename;
use lib dirname($0);
use strict;
use FileHandle;
use Getopt::Long;
use Gff;
$Getopt::Long::ignorecase = 0;

my $outfile;
my $printUsage;
GetOptions("d:s" => \$outfile,
	   "h!" => \$printUsage);

usage() if ($printUsage);

my $gff;
my $fh;
my $name;
while (<>) {
    if (/^LOCUS\s+(\w+)\s+(\d+)/) {
	$name = $1;
	$gff->Print($fh) if (defined($gff));

	if ($fh) { $fh->close(); }
	if (defined($outfile)) {
	    my $fname;
	    $fname = "$outfile/" if $outfile;
	    $fname .= $name . ".gff";
	    $fh = new FileHandle(">$fname") || die "Can't open $fname for writing: $!\n";
	}
	else {
	    $fh = new FileHandle(">&STDOUT");
	}

	$gff = new Gff;
    }
    elsif (/^FEATURES/) {
	chop;
	my $line = <>;
	while ($line !~ /^\S/) {
	    my ($feat,$desc) = ($line =~ /     (\S+)\s+(.*)/);
	    my ($doFrame,$frame) = (0,'.');

	    # calc frame for CDS.  all others are 0
	    ($doFrame,$frame) = (1,0) if ($feat =~ /^CDS$/);

	    while (<>) {
		chop;
		if (/=\"/) {
		    ReadQuoted($fh,$_);
		    next;
		}
		next if (/^\s+\//);
		last if (/^     \w|^\S+/);
		$desc .= substr($_,21);
	    }
	    $line = $_;		# save for next loop through

	    my $strand = '+';
	    if ($desc =~ /complement/) {
		$strand = '-';
		$desc =~ s/complement\((.*)\)/$1/;
	    }
	    if ($desc =~ /join/) {
		$desc =~ s/join\((.*)\)/$1/;
	    }

	    # skip entries with '<' or other stuff.
	    next if ($desc =~ /[^,\s\.\d+]/);

	    my @pieces = split(/,/,$desc);
	    foreach (@pieces) {
		my ($start,$stop) = split(/\.\./,$_);
		$stop = $start if (!defined($stop));

		# calc new frame for bottom strand
		if ($doFrame && $strand eq '-') {
		    $frame = ($frame + ($stop - $start + 1)) % 3;
		}

		new GffLine($name,'genbank',$feat,$start-1,$stop-1,0,
			    $strand,
			    # do GFF style frame
			    ($frame eq '.' || $frame == 0)?$frame:$frame^3)
		  ->Print($fh);

		# calc new frame for top strand
		if ($doFrame && $strand eq '+') {
		    $frame = ($frame + ($stop - $start + 1)) % 3;
		}
	    }

	    if ($doFrame && $frame != 0) {
		die "Bad frame calc for $name.";
	    }
	}
    }
    elsif (/^ORIGIN/) {
	my $seq;
	while (<>) {
	    chop;
	    last if (/^\/\//);
	    s/\s|\d//g;
	    $seq .= $_;
	}
	$gff->SetSeq($name,$seq);
    }
}

$gff->Print($fh) if defined($gff);

# open this file and print the header.  cheapo. 
sub usage {
  open(ME,$0) || die;
  while (<ME>) {
    if (/^[^\#]/) {
      exit;
    }
    s/^.//;
    print STDERR unless /^!/;
  }
}


sub ReadQuoted {
    my ($fh,$line) = @_;
    my ($feat,$rest) = split(/=\"/,$line);
    $feat =~ s/\s+\///;
    if ($rest =~ /\"/) {
	chop $rest;
	print $fh "# $feat: $rest\n";
	return;
    } else {
	print $fh "# $feat: $rest\n";
    }
    while (<>) {
	s/^ {21}//;
	chop;
	if (/\"/) {
	    chop;
	    print $fh "# $_\n";
	    return;
	}
	else {
	    print $fh "# $_\n";
	}
    }
}
< 1 min to Spreed
