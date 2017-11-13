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
	ref			=> '',
	fa			=> '',
	o			=> '',
	genelist	=> '',
	qlx			=> '',
	minhomosize	=> 3,
	readwindow	=> 5,
	nofix		=> '',
	fixhomo		=> '',
	forcefix	=> '',
	gap3window	=> 10,
	h			=> '',
);
GetOptions(
  "ref=s"			=> \$option{ref},
  "fa=s"			=> \$option{fa},
  "o=s"				=> \$option{o},
  "qlx=s"			=> \$option{qlx},
  "genelist=s"		=> \$option{genelist},
  "minhomosize=i"	=> \$option{minhomosize},
  "readwindow=i"	=> \$option{readwindow},
  "nofix"			=> \$option{nofix},
  "fixhomo"			=> \$option{fixhomo},
  "forcefix=f"		=> \$option{forcefix},
  "gap3window=i"	=> \$option{gap3window},
  "h"				=> \$option{h},
) || die("Problem processing command-line options: $!\n");

if($option{h})
{
	print "Usage : perl fixFrameshifts.pl -fa assembly.fa -ref reference.fa -o outputBaseName -genelist reference_genelist.txt\n";
#	print "fixFrameshifts.pl requires MUSCLE to run\n";
	print "fixFrameshifts.pl requires CLUSTALW to run\n";
	print "Options:\n";
	print "-qlx\tRead alignment in Qlx format\n";
	print "-minhomosize(3)\tMinimum length to consider a string of bases an homopolymer. Recommend 2 for 454, 3 for Illumina\n";
	print "-readwindow(5)\tLength of sequence on each side of a frameshift to look at in reads to find a non-frameshifted correction\n";
	print "-nofix\tUsing -nofix will still look for frameshifts and write the log, but it will not automatically correct them in the assembly\n";
	print "-fixhomo\tForces correction of frameshifts in homopolymer regions. Will work without reads supplied, but requires at least 1 read support to fix if reads are supplied\n";
	print "-forcefix\tWill fix frameshifts as long as one of the correct length windows is present above # (supplied, i.e. 0.25) fraction of the reads. Off by default\n";
	exit();
}

my $inputFasta = $option{fa};
my $output = $option{o};
my %consseqs;
my $consseq;

#my $musclepath = "/project/umw_jeffrey_bailey/share/bin_sync/Muscle/";
my $clustalpath = "/project/umw_jeffrey_bailey/share/bin_sync/clustalw-2.1/bin/";

# Truncate Name...
open(FA, $inputFasta) || die("Unable to open $inputFasta\n");
my $count = 1;
my $curid = '';
while(my $line = <FA>)
{
	chomp $line;
	if($line =~ />.+/)
	{
		$curid = 'Seq'.$count;
		$count++;
	}else{
		$consseqs{$curid} .= $line;
	}
}

if(scalar(keys %consseqs) == 1){$consseq = $consseqs{'Seq1'};}

# Get Geneorder
open(GENELIST, $option{genelist}) || die("Unable to open ".$option{genelist}."\n");
my @geneorder;
my %genes;

while(my $line = <GENELIST>)
{
	chomp $line;
	if($line =~ /(.+?)\t(.+?)\t(.+)/)
	{
		push(@geneorder, $1);
		$genes{$1}{start} = $2;
		$genes{$1}{stop} = $3;
	}
}
my $refseq = readSeq($option{ref});
(my $alignedRef, my $alignedCons) = alignPair($refseq, $consseq, $output."_alignPair");

my %readSeqs;
my %readWindows;
my %potentialFs;
my @fsWindows;

my $isFs = isFrameshifted($alignedCons, $alignedRef);
open(FSLOG, ">".$output."_frameshifts_log.txt");
my $newAssembly = $consseq;

open(ASSEMBLYOUT, ">".$output."_fixed_assembly.fa");
print ASSEMBLYOUT ">".$output."\n";

### FIX FRAMESHIFT WITHOUT READS
unless($option{qlx})
{
	print "No .qlx file supplied, fixFrameshift running without read support...\n";
	if($isFs)
	{
		getAlignWindows();
		my $posmodifier = 0;
		my %fixingwindows;
		foreach my $fswindow (@fsWindows)
		{
			print FSLOG "Potential FS Window starting at ".$$fswindow{start}." in gene ".$$fswindow{gene}." with gap of length ".$$fswindow{len}."\n";
			print FSLOG "It contains ".scalar(keys %{$$fswindow{allfs}})." frameshifts:\n";
			foreach my $curfs (sort keys %{$$fswindow{allfs}})
			{
				print FSLOG "$curfs ";
				if($$fswindow{allfs}{$curfs}{homo}){
					print FSLOG "in homopolymer\n";
				}else{
					print FSLOG "out of homopolymer\n";
				}
			}
			print FSLOG "\nReference Window: ".$$fswindow{alnref}."\n";
			print FSLOG "Assembly Window : ".$$fswindow{alncons}."\n";

			if($option{fixhomo})
			{
				$newAssembly = '';
				foreach my $fspos (sort {$b <=> $a} keys %{$$fswindow{allfs}})
				{
					my $fsposgaplen = $$fswindow{allfs}{$fspos}{len};
					if($$fswindow{allfs}{$fspos}{homo})
					{
						print FSLOG "Homopolymer Error : potential fix, ";
						
						if($fsposgaplen < 0){print FSLOG "add ";}else{print FSLOG "remove ";}
						print FSLOG abs($fsposgaplen)." ".$$fswindow{allfs}{$fspos}{homo}." ";
						print FSLOG "at position ".$fspos."\n";
						
						my $conswindow = $$fswindow{alncons};
						if($fixingwindows{$$fswindow{start}}){$conswindow = $fixingwindows{$$fswindow{start}}{seq};}
						
						my $fixedwindow = '';
						my $curpos = $$fswindow{start} - 1;
						for(my $alnpos = 0; $alnpos < length($conswindow); $alnpos++)
						{
							my $curcons = substr($conswindow, $alnpos, 1);
							if($curcons ne '-'){$curpos++;}
							if($curpos < $fspos){
								$fixedwindow .= $curcons;
							}elsif($curpos == $fspos)
							{
								if($$fswindow{len} < 0)
								{
									$fixedwindow .= $curcons;
									$fixedwindow .= $$fswindow{allfs}{$fspos}{homo} x abs($fsposgaplen);
									$alnpos += abs($fsposgaplen);
								}else{
									for(my $i = 0; $i < abs($fsposgaplen) - 1; $i++)
									{
										$alnpos++;
										$curpos++;
									}
								}
							}else{
								$fixedwindow .= $curcons;
							}
						}
						print FSLOG "Fixed Window : ".$fixedwindow."\n";
						$fixingwindows{$$fswindow{start}}{seq} = $fixedwindow;
						$fixingwindows{$$fswindow{start}}{stop} = $$fswindow{stop};
						print FSLOG "FIXED\n";
					}
				}
			}else{
				print FSLOG "UNFIXED\n";
			}
		
			print FSLOG "\n**********************************\n";
		}
		
		if($option{fixhomo})
		{
			for(my $conspos = 1; $conspos <= length($consseq); $conspos++)
			{
				if($fixingwindows{$conspos})
				{
					$fixingwindows{$consseq}{seq} =~ s/\-//;
#					print "$conspos\t".$fixingwindows{$conspos}{seq}."\t".substr($consseq, $conspos - 1, ($fixingwindows{$conspos}{stop} - $conspos) + 1)."\n";
					$newAssembly .= $fixingwindows{$conspos}{seq};
					$conspos += ($fixingwindows{$conspos}{stop} - $conspos);
				}else{
					$newAssembly .= substr($consseq, $conspos - 1, 1);
				}
			}
		}
	}
	print ASSEMBLYOUT $newAssembly."\n";
}else{

### FIX FRAMESHIFT WITH READS
	getReadSeqs($option{qlx});
	if($isFs)
	{
		my %fsFixes;
		foreach my $fswindow (@fsWindows)
		{
			print FSLOG "Potential FS Window starting at ".$$fswindow{start}." in gene ".$$fswindow{gene}." with gap of length ".$$fswindow{len}."\n";
			print FSLOG "It contains ".scalar(keys %{$$fswindow{allfs}})." frameshifts:\n";
			foreach my $curfs (sort keys %{$$fswindow{allfs}})
			{
				print FSLOG "$curfs ";
				if($$fswindow{allfs}{$curfs}{homo}){
					print FSLOG "in homopolymer\n";
				}else{
					print FSLOG "out of homopolymer\n";
				}
			}
			
			my $bestWindow = '';
			my $secondary = 0;
	
			print FSLOG "\nAssembly Window :\n".$readWindows{$$fswindow{start}}{ref}."\n";
			
			print FSLOG "\nTop 3 Read Windows:\n";
			
			my $count = 0;
			foreach my $fullwindow (sort {$readWindows{$$fswindow{start}}{full}{$b} <=> $readWindows{$$fswindow{start}}{full}{$a}} keys %{$readWindows{$$fswindow{start}}{full}})
			{
				if($fullwindow eq 'count'){next;}
				if($count < 3){
					print FSLOG $fullwindow."\t".$readWindows{$$fswindow{start}}{full}{$fullwindow}."\n";
				}
				if(length($fullwindow) == ($$fswindow{stop} - $$fswindow{start} + 1 - $$fswindow{len}))
				{
					unless($bestWindow)
					{
						$bestWindow = $fullwindow;
					}
				}else{
					unless($bestWindow)
					{
						$secondary = 1;
					}
				}
				$count++;
			}
			if($bestWindow){
				print FSLOG "\nBest Replacement Window : $bestWindow (".$readWindows{$$fswindow{start}}{full}{$bestWindow}."/".$readWindows{$$fswindow{start}}{full}{count}.")\n";
				my $accepted = 1;
				if($secondary){
					print FSLOG "* Secondary Window *\n";
					$accepted = 0;
					if($option{forcefix} && ($readWindows{$$fswindow{start}}{full}{$bestWindow}/$readWindows{$$fswindow{start}}{full}{count}) > $option{forcefix})
					{
						$accepted = 1;
					}
				}
				
				if($accepted == 0)
				{
					foreach my $curfs (sort keys %{$$fswindow{allfs}})
					{
						if($$fswindow{allfs}{$curfs}{homo} && $option{fixhomo}){$accepted = 1;}
						else{$accepted = 0; last;}
					}
				}
				if($option{nofix}){$accepted = 0;}
				if($accepted)
				{
					$fsFixes{$$fswindow{start}}{length} = length($readWindows{$$fswindow{start}}{ref});
					$fsFixes{$$fswindow{start}}{window} = $bestWindow;
					$newAssembly = fixAssembly($newAssembly, $$fswindow{start}, $fsFixes{$$fswindow{start}}{length}, $fsFixes{$$fswindow{start}}{window});
					print FSLOG "FIXED\n";
				}else{
					print FSLOG "UNFIXED\n";
				}
			}else{
					print FSLOG "UNFIXED\n";
			}
			print FSLOG "\n**********************************\n";
		}
	}
	print ASSEMBLYOUT $newAssembly."\n";
}
	
($alignedRef, $alignedCons) = alignPair($refseq, $newAssembly, $output."_alignPair_fixed");


undef @fsWindows;
$isFs = isFrameshifted($alignedCons, $alignedRef);
if($isFs){
	print "Assembly still contains Frameshifts (see log for details)\n";}
else{print "All Frameshifts have been fixed (see log for details)\n";}



sub getAlignWindows
{
	foreach my $fswindow (@fsWindows)
	{
		my $curconsstr = '';
		my $currefstr = '';
		my $conspos = 0;
		
		for(my $alnpos = 0; $alnpos < length($alignedCons); $alnpos++)
		{
			my $curcons = substr($alignedCons, $alnpos, 1);
			my $curref = substr($alignedRef, $alnpos, 1);
			
			if($curcons ne '-'){$conspos++;}
			if($conspos >= $$fswindow{start} && $conspos <= $$fswindow{stop})
			{
				$curconsstr .= $curcons;
				$currefstr .= $curref;
			}
		}
#		print $$fswindow{start}."\t".$curconsstr."\t".$currefstr."\n";
		$$fswindow{alncons} = $curconsstr;
		$$fswindow{alnref} = $currefstr;
	}
}



sub fixAssembly
{
	my $assembly = shift;
	my $startPos = shift;
	my $windowLength = shift;
	my $fixWindow = shift;
	
	return substr($assembly, 0, $startPos - 1).$fixWindow.substr($assembly, $startPos + $windowLength - 1);
}

sub mergeFs
{
	my $prevfs = 100000000;
	my $fswinflag = -1;
	
	foreach my $fspos (sort {$b <=> $a} keys %potentialFs)
	{
		if(($prevfs - $fspos) > ($option{readwindow} * 2))
		{
			$fswinflag++;
			$fsWindows[$fswinflag]{start} = $fspos - $option{readwindow};
			if($fsWindows[$fswinflag]{start} < 0){$fsWindows[$fswinflag]{start} = 0;}
			$fsWindows[$fswinflag]{stop} = $fspos + $potentialFs{$fspos}{len} + $option{readwindow};
			if($fsWindows[$fswinflag]{stop} >= length($consseq)){$fsWindows[$fswinflag]{stop} = (length($consseq) - 1);}
			$fsWindows[$fswinflag]{allfs}{$fspos}{homo} = $potentialFs{$fspos}{homo};
			$fsWindows[$fswinflag]{allfs}{$fspos}{len} = $potentialFs{$fspos}{len};
			$fsWindows[$fswinflag]{gene} = $potentialFs{$fspos}{gene};
			$fsWindows[$fswinflag]{len} = $potentialFs{$fspos}{len};
		}else{
			$fsWindows[$fswinflag]{start} = $fspos - $option{readwindow};
			if($fsWindows[$fswinflag]{start} < 0){$fsWindows[$fswinflag]{start} = 0;}
			$fsWindows[$fswinflag]{allfs}{$fspos}{homo} = $potentialFs{$fspos}{homo};
			$fsWindows[$fswinflag]{allfs}{$fspos}{len} = $potentialFs{$fspos}{len};
			if($potentialFs{$fspos}{gene} ne $fsWindows[$fswinflag]{gene})
			{
				$fsWindows[$fswinflag]{gene} .= '-'.$potentialFs{$fspos}{gene};
			}
			$fsWindows[$fswinflag]{len} += $potentialFs{$fspos}{len};
		}
		$prevfs = $fspos;
	}

	my $arrayindex = 0;
	foreach my $window (@fsWindows)
	{
#		print "Checking ".$$window{start}."...\n";
		if($$window{len} == 0)
		{
#			print "Removing\n";
			splice (@fsWindows, $arrayindex, 1);
		}else{
#			print "Keeping\n";
			$arrayindex++;
		}
	}
#	print "Returning : ".scalar(@fsWindows)."\n";
	return scalar(@fsWindows);
}

sub getReadSeqs
{
	my $qlxfile = shift;
	my $readid;
	my $readstart;
	my $readend;
	my $readlen;
	my $strand;
	my $ref;
	my $refstart;
	my $refend;
	
	my $refstr = '';
	my $readstr = '';
	my $qualstr = '';
	my $mode;
	
	open(QLX, $qlxfile);
	
	while (my $line = <QLX>)
	{
		chomp $line;
	#	print $line;
		if($line =~ />Read/)
		{
			my @linedata = split(/ /,$line);
			$readid = $linedata[1];
			$readstart = $linedata[2];
			$readend = $linedata[3];
			$readlen = $linedata[4];
			$strand = $linedata[5];
			$ref = $linedata[6];
			$refstart = $linedata[7];
			$refend = $linedata[8];
			$mode = 1;
		}elsif($mode == 1)
		{
			$refstr = $line;
			$mode = 2;
		}elsif($mode == 2)
		{
			$readstr = $line;
			foreach my $fswindow(@fsWindows)
			{
				if($refstart <= $$fswindow{start} && $refend >= $$fswindow{stop})
				{
					getReadWindow($readstr, $refstr, $refstart, $refend, $$fswindow{start}, $$fswindow{stop});
				}
			}
			$mode = 3;
		}elsif($mode == 3)
		{
			$mode = 4;
		}elsif($mode == 4)
		{
			$qualstr = $line;
			$mode = 0;
		}
	}
}

sub getReadWindow
{
	my $string = shift;
	my $refstring = shift;
	my $readpos = shift;
	my $readend = shift;
	my $winstart = shift;
	my $winstop = shift;
	
	my $window = '';
	my $refwindow = '';
	
	my $type = 'full';
	if($readpos > $winstart || $readend < $winstop)
	{
		$type = 'partial';
	}
	
	
	for(my $alnpos = 0; $alnpos < length($string); $alnpos++)
	{
		my $curres = substr($string, $alnpos, 1);
		my $curref = substr($refstring, $alnpos, 1);
		
		if($readpos >= $winstart && $readpos <= $winstop)
		{
			unless($curres eq '-')
			{
				$window .= $curres;
			}
			if($type eq 'full' && $curref ne '-')
			{
				$refwindow .= $curref;
			}	
		}
		if($curref ne '-'){
			$readpos++;
		}
	}
	$readWindows{$winstart}{$type}{$window}++;
	$readWindows{$winstart}{$type}{'count'}++;
	unless($readWindows{$winstart}{ref})
	{
		$readWindows{$winstart}{ref} = $refwindow;
	}
}

sub isFrameshifted
{
	my $fsstr = shift;
	my $refstr = shift;
	
	my $fspos = 0;
	my $refpos = 0;
	undef %potentialFs;
	
#	print ">FSSTR\n".$fsstr."\n\n>REFSTR\n".$refstr."\n\n";
	
	for(my $alnpos = 0; $alnpos < length($fsstr); $alnpos++)
	{
		my $curfs = substr($fsstr, $alnpos, 1);
		my $curref = substr($refstr, $alnpos, 1);
		my $inhomo = 0;
		if($curfs ne '-'){$fspos++;}
		if($curref ne '-'){$refpos++;}

	#	print $inhomo."\n";
#		print $fspos."\t".$refpos."\n";
		if($fspos == 0 || $refpos == 0){next;}
		if($curfs eq '-' || $curref eq '-')
		{
			$inhomo = inHomo($alnpos, $fsstr, $refstr);
			unless($inhomo)
			{
				$inhomo = inHomo($alnpos, $refstr, $fsstr);
			}
#			print "INHOMO? $inhomo\n";
			my $curgene = inGene($refpos);
			unless($curgene){next;}
			my $gapdiff =  gapArea($refstr, $alnpos) - gapArea($fsstr, $alnpos);
#			print $alnpos."\t".$gapdiff."\n";
			if($gapdiff % 3 != 0)
			{
			#	$frameshifted = 1;
				$potentialFs{$fspos}{len} = $gapdiff;
				$potentialFs{$fspos}{homo} = $inhomo;
				$potentialFs{$fspos}{gene} = $curgene;
				$potentialFs{$fspos}{alnpos} = $alnpos;
			}
			if($curfs eq '-'){
				my $gaplen = gapLength($fsstr, $alnpos);
				$alnpos += $gaplen - 1;
				$refpos += $gaplen - 1;
			}elsif($curref eq '-'){
				my $gaplen = gapLength($refstr, $alnpos);
				$alnpos += $gaplen - 1;
				$fspos += $gaplen - 1;
			}
		}
	}
	my $frameshifted = mergeFs();
	return $frameshifted;
}


sub inGene
{
	my $pos = shift;
	
	foreach my $gene (keys %genes)
	{
		if($pos >= $genes{$gene}{start} && $pos <= $genes{$gene}{stop})
		{
			return $gene;
		}
	}
	return 0;
}

sub alignPair
{
	my $seq1 = shift;
	my $seq2 = shift;
	my $alignPairOut = shift;
	
	open(MFA, ">$alignPairOut.mfa");
	print MFA ">Ref\n$seq1\n>Assembly\n$seq2\n";
	close MFA;
	
#	system($musclepath."muscle -in $alignPairOut.mfa -out $alignPairOut.afa -quiet");
	system($clustalpath."clustalw2 -INFILE=$alignPairOut.mfa -OUTFILE=$alignPairOut.afa -ALIGN -TYPE=DNA -OUTPUT=FASTA");
	
	open(AFA, "$alignPairOut.afa");
	my $curid;
	my $alignseq1;
	my $alignseq2;
	while(my $line = <AFA>)
	{
		chomp $line;
		if($line =~ />(.+)/)
		{
			$curid = $1;
		}elsif($curid eq 'Ref')
		{
			$alignseq1 .= $line;
		}elsif($curid eq 'Assembly')
		{
			$alignseq2 .= $line;
		}
	}

	return $alignseq1, $alignseq2;
}

sub readSeq
{
	my $file = shift;
	my %seqs;
	my $curid = '';
	open(FILE, $file);
	while(my $line = <FILE>)
	{
		chomp $line;
		if($line =~ />(.+)/)
		{
			$curid = $1;
			next;
		}else{
			$seqs{$curid} .= $line;
		}
	}
	if(scalar(keys %seqs) == 1)
	{
		return $seqs{$curid};
	}else{
		return \%seqs;
	}
}

sub inHomo
{
	my $postocheck = shift;
	my $dnastring = shift;
	my $compstring = shift;
	my $homostrres = 1;
	my $homores = substr($dnastring, $postocheck, 1);

	my $previousmatch = 0;

	for(my $homopos = $postocheck - 1; $homopos >= $postocheck - $option{minhomosize}; $homopos--)
	{
		if($homopos < 0){last;}

		my $resnow = substr($dnastring, $homopos, 1);
		#print $resnow."\n";
		my $compres = substr($compstring, $homopos, 1);

		if($resnow eq $homores)
		{
			if($compres eq $resnow){$previousmatch = 1;}
			$homostrres++;
		}else{last;}

		if($homostrres == $option{minhomosize} && $previousmatch)
		{
			return $homores;
		}
	}
	
	for(my $homopos = $postocheck + 1; $homopos <= $postocheck + $option{minhomosize}; $homopos++)
	{
		if($homopos >= length($dnastring)){last;}
		
		my $resnow = substr($dnastring, $homopos, 1);
		my $compres = substr($compstring, $homopos, 1);

		if($resnow eq $homores)
		{
			if($compres eq $resnow){$previousmatch = 1;}
			$homostrres++;
		}else{last;}
		
		if($homostrres == $option{minhomosize} && $previousmatch)
		{
			return $homores;
		}
	}
	return 0;
}

### CALCULATE THE LENGTH OF A GAP
sub gapLength
{
	my $dnastring = shift;
	my $postocheck = shift;

	my $resnow = "-";
	my $gaplength = 1;

	my $flag = 1;
	while($resnow eq "-")
	{
		$resnow = substr($dnastring, ($postocheck - $flag), 1);
		if($resnow eq "-")
		{
			$gaplength++;
			$flag++;
		}
	}

	$flag = 1;
	$resnow = "-";
	while($resnow eq "-")
	{
		$resnow = substr($dnastring, ($postocheck + $flag), 1);
		if($resnow eq "-")
		{
			$gaplength++;
			$flag++;
		}
	}
	return $gaplength;
}

sub gapArea
{
	my $dnastring = shift;
	my $postocheck = shift;
	my $totalgap = 0;
	my $areasize = $option{gap3window};

	my $curpos = $postocheck;
	
	if(terminalGap($dnastring, $postocheck)){return 0;}
	
	my $initialgaplen = 0;
	my $initialposlow;
	my $initialposhigh;
	
	while(my $resnow = substr($dnastring, $curpos, 1))
	{
		if($resnow eq "-")
		{
			$initialgaplen++;
		}else{
			$initialposlow = $curpos;
			last;
		}
		$curpos--;
		if($curpos < 0)
		{
			return 0;
		}
	}
	
	$curpos = $postocheck + 1;
	while(my $resnow = substr($dnastring, $curpos, 1))
	{
		if($resnow eq "-")
		{
			$initialgaplen++;
		}else{
			$initialposhigh = $curpos;
			last;
		}
		$curpos++;
		if($curpos > length($dnastring))
		{
			return 0;
		}
	}

	if($initialgaplen % 3 != 0)
	{
		$totalgap += $initialgaplen;
	}
	
	my $curgaplen = 0;

	$curpos = $initialposlow;
	
	while(1)
	{
		if(substr($dnastring, $curpos, 1) eq "-")
		{
			$curgaplen++;
		}else{
			if($curgaplen)
			{
				if($curgaplen % 3 != 0)
				{
					unless(terminalGap($dnastring, $curpos + 1))
					{
						$totalgap += $curgaplen;
					}
				}
				$curgaplen = 0;
			}
			if(abs($postocheck - $curpos) > $areasize)
			{
				last;
			}
		}
		$curpos--;
		if($curpos < 0)
		{
			last;
		}
	}

	$curpos = $initialposhigh;
	$curgaplen = 0;
	while(1)
	{
		if(substr($dnastring, $curpos, 1) eq "-")
		{
			$curgaplen++;
		}else{
			if($curgaplen)
			{
				if($curgaplen % 3 != 0)
				{
					$totalgap += $curgaplen;
				}
				$curgaplen = 0;
			}
			if(abs($postocheck - $curpos) > $areasize)
			{
				last;
			}
		}
		$curpos++;
		if($curpos > length($dnastring))
		{
			last;
		}	
	}

	return $totalgap;
}

sub terminalGap
{
	my $dnastring = shift;
	my $postocheck = shift;
#	print "Terminal $postocheck\t".length($dnastring)."\n";
	my $resnow = "-";
	my $flag = 1;
	while($resnow eq "-")
	{
		$resnow = substr($dnastring, ($postocheck - $flag), 1);
#		print "$resnow ".($postocheck - $flag)."\n";
		$flag++;
		if($postocheck - $flag < 0){return 1;}
	}
	
	$flag = 1;
	$resnow = "-";
	while($resnow eq "-")
	{
		$resnow = substr($dnastring, ($postocheck + $flag), 1);
#		print "$resnow ".($postocheck + $flag)."\n";
		$flag++;
		if(($postocheck + $flag) >= (length($dnastring) - 1)){return 1;}
	}
	
	return 0;
}
