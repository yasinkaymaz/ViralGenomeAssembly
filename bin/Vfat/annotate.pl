#!/share/pkg/perl/5.18.1/bin/ perl -w

# Copyright � 2012 The Broad Institute, Inc.
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
	pepfolder	=> '',
	o			=> '',
	align		=> '',
	genelist	=> '',
	maxannogaplen	=> 50,
	minpctid	=> 0.4,
	minpctlen	=> 0.2,
	sampname	=> '',
	h			=> '',
);
GetOptions(
  "ref=s"		=> \$option{ref},
  "fa=s"		=> \$option{fa},
  "pepfolder=s"	=> \$option{pepfolder},
  "o=s"			=> \$option{o},
  "align=s"		=> \$option{align},
  "genelist=s"	=> \$option{genelist},
  "maxannogaplen=i"	=> \$option{maxannogaplen},
  "minpctid=f"	=> \$option{minpctid},
  "minpctlen=f"	=> \$option{minpctlen},
  "sampname=s"	=> \$option{sampname},
  "h"			=> \$option{h},
) || die("Problem processing command-line options: $!\n");

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

#### MUST DO 'setenv WISECONFIGDIR /seq/annotation/bio_tools/GeneWise/wise2.2.0/wisecfg/' BEFORE SCRIPT ####
my $genewisepath = "/project/umw_jeffrey_bailey/share/bin_sync/wise2.2.0/src/bin/";
my $genewisecfgpath = "/project/umw_jeffrey_bailey/share/bin_sync/wise2.2.0/wisecfg/";
#my $musclepath = "/project/umw_jeffrey_bailey/share/bin_sync/Muscle/";
my $clustalpath = "/project/umw_jeffrey_bailey/share/bin_sync/clustalw-2.1/bin/";
my $mafftpath = "/project/umw_jeffrey_bailey/share/bin_sync/mafft-7.215-with-extensions/scripts/";

$ENV{WISECONFIGDIR} = $genewisecfgpath;

if($option{h})
{
	print "Usage : perl annotate.pl -fa input.fa -ref reference.fa -o outputBaseName -pepfolder peptideFolder -genelist reference_genelist.txt\n";
	print "annotate.pl requires	GeneWise (wise2.2) to run\n";
	print "Options:\n";
	print "-maxannogaplen(50)\tMaximum gap in the alignment before splitting gene\n";
	print "-minpctlen(0.2)\t\tMinimum % length of a gene aligning to consider it partial instead of missing\n";
	print "-minpctid(0.4)\t\tMinimum % identity of a gene aligning to consider it partial instead of missing\n";
	print "-sampname()\t\tSpecifies the sample name to write in NCBI submission files. If none is specified, uses the output name (-o)\n";
	exit();
}

my $inputFasta = $option{fa};
my $output = $option{o};
my $inputFolder = $option{pepfolder};
my $inputAlignment = $option{align};

unless($option{sampname}){$option{sampname} = $output;}

if(substr($inputFolder, length($inputFolder) - 1, 1) eq '/'){chop $inputFolder;}
my %consseqs;
my $consseq;

# Truncate Name...
open(FA, $inputFasta) || die ("Unable to open $inputFasta\n");
my $tmpFa = "$output"."_tmpFasta.fa";
open(FATMP, ">$tmpFa");
my $count = 1;
my $curid = '';
open(NCBIFA, ">$output"."_assembly_ncbi.fa");
my $seqname = '';
while(my $line = <FA>)
{
	chomp $line;
	if($line =~ />(.+)/)
	{
		print FATMP ">Seq"."$count\n";
		$curid = 'Seq'.$count;
		$consseqs{$curid}{name} .= $1;
		$count++;
	}else{
		print FATMP $line."\n";
		$consseqs{$curid}{seq} .= $line;
	}
}
if(scalar(keys %consseqs) == 1){$consseq = $consseqs{'Seq1'}{seq};}else{
	print "annotate.pl is currently designed to handle a single sequence as a reference only\n";
	exit();
}

print NCBIFA ">".$option{sampname}."\n";
for(my $pos = 0; $pos < length($consseq) - 1; $pos++)
{
	my $ntpair = substr($consseq, $pos, 2);
	my $curnt = substr($consseq, $pos, 1);

	if($ntpair eq 'NN')
	{
		print NCBIFA "\n>?";
		my $Nlen = 1;
		while(substr($consseq, $pos, 1) eq 'N')
		{
			$pos++;
			$Nlen++;
		}
		$pos--;
		$Nlen--;
		print NCBIFA $Nlen."\n";
	}else{
		print NCBIFA $curnt;
	}
}

open(DETAILSOUT, ">$output"."_details.txt");
open(SUMMARY, ">$output"."_summary.txt");
open(NCBI, ">$output"."_ncbi.txt");

# Get Ref Gene data
open(GENELIST, $option{genelist}) || die("Unable to open ".$option{genelist}."\n");
my @geneorder;
my %genelist;
while(my $line = <GENELIST>)
{
	chomp $line;
	if($line =~ /(.+?)\t(.+?)\t(.+)/)
	{
		my $genename = $1;
		my $startpos = $2;
		my $stoppos = $3;
		my $exonno = 1;

		if($genename =~ /(.+)\_exon(.+)/)
		{
			$genename = $1;
			$exonno = $2;
		}

		if($exonno == 1){push(@geneorder, $1)};
		$genelist{$genename}{nbexon}++;
		$genelist{$genename}{'ex'.$exonno}{start} = $startpos;
		$genelist{$genename}{'ex'.$exonno}{stop} = $stoppos;
	}
}

# Fill Geneorder
my %refGeneorder;
foreach my $gene (keys %genelist)
{
	foreach my $exon (keys %{$genelist{$gene}})
	{
		if($exon eq 'nbexon'){next;}
		foreach my $compgene (keys %genelist)
		{
			foreach my $compexon (keys %{$genelist{$compgene}})
			{
				if($compexon eq 'nbexon'){next;}
				if($compgene eq $gene && $compexon eq $exon){next;}

				$refGeneorder{$gene}{$exon}{$compgene}{$compexon}{startvstart} = compValue($genelist{$gene}{$exon}{start}, $genelist{$compgene}{$compexon}{start});
				$refGeneorder{$gene}{$exon}{$compgene}{$compexon}{startvstop} = compValue($genelist{$gene}{$exon}{start}, $genelist{$compgene}{$compexon}{stop});
				$refGeneorder{$gene}{$exon}{$compgene}{$compexon}{stopvstart} = compValue($genelist{$gene}{$exon}{stop}, $genelist{$compgene}{$compexon}{start});
				$refGeneorder{$gene}{$exon}{$compgene}{$compexon}{stopvstop} = compValue($genelist{$gene}{$exon}{stop}, $genelist{$compgene}{$compexon}{stop});
#				print "$gene\t$exon\t$compgene\t$compexon\t".$refGeneorder{$gene}{$exon}{$compgene}{$compexon}{startvstart}."\t".$refGeneorder{$gene}{$exon}{$compgene}{$compexon}{startvstop}."\t".$refGeneorder{$gene}{$exon}{$compgene}{$compexon}{stopvstart}."\t".$refGeneorder{$gene}{$exon}{$compgene}{$compexon}{stopvstop}."\n";
			}
		}
	}
}



#die;

# Get reference sequence
my $refseq = readSeq($option{ref});

# Get alignment sequences

my $alignref;
my $alignass;
if($option{align})
{
	($alignref, $alignass) = readAlign($option{align});
}else{
	($alignref, $alignass) = alignPair($refseq, $consseq, "$output"."_aligned");
}

my %geneOrderPositions;
my %alignGenelist;
genelistFromAlignment($alignref, $alignass);


# Run GeneWise
unless(-d $output)
{
	system("mkdir $output");
}
opendir(AAS, $inputFolder) || die("Unable to open folder $inputFolder\n");
my %genefeat;
my %finalGeneStrings;

my $cdsstat = "Multi";
while(my $file = readdir(AAS))
{
	if($file =~ /(.+)\_pep\.fa/)
	{
		my $gene = $1;
		open(GENE, $inputFolder."/".$file);

		my $aastr = '';
		while(my $line = <GENE>)
		{
			chomp $line;
			if($line =~ />/)
			{
				next;
			}else{
				$aastr .= $line;
			}
		}
		$genelist{$1}{size} = length($aastr);
		$genelist{$1}{seq} = $aastr;

#		system($genewisepath."genewise $inputFolder/$file ".$option{fa}." > $output/$gene"."_genewise_genes.txt -genes");
#		system($genewisepath."genewise $inputFolder/$file $tmpFa > $output/$gene"."_genewise.txt");
		system($genewisepath."genewise $inputFolder/$file $inputFasta -both -pep > $output/$gene"."_genewisePEP.txt");
		system($genewisepath."genewise $inputFolder/$file $inputFasta -both -cdna > $output/$gene"."_genewiseDNA.txt");

	}elsif($file =~ /Peptides_Features.txt/)
	{
		open(AAFEAT, $inputFolder."/".$file);
		while(my $line = <AAFEAT>)
		{
			chomp $line;
			if($line =~ /CDSStatus\t(.+)/)
			{
				if(lc($1) eq 'single'){
					$cdsstat = 'Single';
				}
			}elsif($line =~ /(.+?)\t(.+?)\t(.+?)\t(.+)/)
			{
				if($2 eq 'Start'){next;}
				$genefeat{$1}{start} = $2;
				$genefeat{$1}{stop} = $3;
				$genefeat{$1}{product} = $4;
			}
		}
	}
}

my %genes;
my %geneFlags;

open(OUTPUT, ">$output"."_genelist.txt");
my $newgoflag = 0;

foreach my $gene (@geneorder)
{
	my $nbexon = readGene("$output/$gene"."_genewise_genes.txt", $gene, $genelist{$gene}{nbexon});
#Temp!!!
	print $gene."\t".$genes{$gene}{nbexon}." exon(s)\n";
#	for(my $i = 1; $i <= $genes{$gene}{nbexon}; $i++)
#	{
#Temp!!!
#		print "Exon $i\t".$genes{$gene}{'ex'.$i}{start}."\t".$genes{$gene}{'ex'.$i}{stop}."\t".$genes{$gene}{frameshifted}."\n";
#	}
#Temp!!!Here is giving error: illegal division by zero
	my $completegene = readGenewise("$output/$gene"."_genewise.txt", $gene);
	my $flagstatus = checkGeneFlags($gene, $completegene);

#	print "Flags\t".$geneFlags{$gene}{size}."\t".$geneFlags{$gene}{startstop}."\t".$geneFlags{$gene}{completeness}."\t".$geneFlags{$gene}{exons}."\t".$geneFlags{$gene}{frameshift}."\n\n";
	$genes{$gene}{flags} = $flagstatus;
}

my $okgo = checkGeneorder();

unless($okgo){
#	print "Bad final geneorder!\n";
}

system("rm $tmpFa");


### print OUTPUTs

print SUMMARY "Gene\tNbExon\tCompletion\tStart\tStop\tFlag Size\tFlag Start/Stop\tFlag Frameshift\tFlag Internal N String\tFlag Exons\tAll Flags Ok\n";
foreach my $gene (@geneorder)
{
	my $complete = 'Full';
	if($geneFlags{$gene}{completeness} == 0){$complete = 'Partial';}
	elsif($geneFlags{$gene}{completeness} == -1){$complete = 'Missing';}

	if($complete eq 'Missing')
	{
		print SUMMARY $gene."\t0\t".$complete."\tN/A\tN/A\t0\t0\t0\t0\t0\t0\n";
	}else{

		my $start = 'N/A';
		my $stop = 'N/A';

		if($genes{$gene}{ex1}{start}){$start = $genes{$gene}{ex1}{start};}
			if($genes{$gene}{'ex'.$genes{$gene}{nbexon}}{stop}){$stop = $genes{$gene}{ex1}{stop};}

		print SUMMARY $gene."\t".$genes{$gene}{nbexon}."\t".$complete."\t".$start."\t".$stop."\t".$geneFlags{$gene}{size}."\t".$geneFlags{$gene}{startstop}."\t".$geneFlags{$gene}{frameshift}."\t".$geneFlags{$gene}{Nstr}."\t".$geneFlags{$gene}{exons}."\t".$genes{$gene}{flags}."\n";
	}
}

foreach my $gene (@geneorder)
{
	if($geneFlags{$gene}{completeness} == -1)
	{
		print DETAILSOUT "$gene\t0\t0\tMissing\n";
	}elsif($genes{$gene}{nbexon} == 1)
	{
		print DETAILSOUT $gene."\t".$genes{$gene}{ex1}{start}."\t".$genes{$gene}{ex1}{stop}."\t";
		print OUTPUT $gene."\t".$genes{$gene}{ex1}{start}."\t".$genes{$gene}{ex1}{stop}."\n";
		if($geneFlags{$gene}{completeness} == 1){
			print DETAILSOUT "Full\n";
		}else{
			print DETAILSOUT "Partial\n";
		}
	}else{
		for(my $i = 1; $i <= $genes{$gene}{nbexon}; $i++)
		{
			print DETAILSOUT $gene."_ex".$i."\t".$genes{$gene}{'ex'.$i}{start}."\t".$genes{$gene}{'ex'.$i}{stop}."\t";
			print OUTPUT $gene."_ex".$i."\t".$genes{$gene}{'ex'.$i}{start}."\t".$genes{$gene}{'ex'.$i}{stop}."\n";
			if($geneFlags{$gene}{completeness} == 1){
				print DETAILSOUT "Full\n";
			}else{
				print DETAILSOUT "Partial\n";
			}
		}
	}
}

### Print NCBI

print NCBI ">Feature\t".$option{sampname}."\n";
if($geneFlags{$geneorder[0]}{completeness} == 1 && $genes{$geneorder[0]}{ex1}{start} > 1)
{
	print NCBI "1\t".($genes{$geneorder[0]}{ex1}{start} - 1)."\t5'UTR\n";
	print NCBI "\t\t\tnote\tindels in UTR have not been validated\n";
}

if($cdsstat eq 'Single')
{
	my $firstgene = 1;
	my $firststart = 0;
	my $laststop = 0;
	foreach my $gene (@geneorder)
	{
		if($geneFlags{$gene}{completeness} == -1)
		{
			$firstgene = 0;
			next;
		}elsif($genes{$gene}{nbexon} == 1)
		{
			my $ncbistart = $genes{$gene}{ex1}{start};
			my $ncbistop = $genes{$gene}{ex1}{stop};
			if($genes{$gene}{partstart} == 1 || $firstgene == 0){$ncbistart = "<".$ncbistart;}
			if($genes{$gene}{partstop} == 1){$ncbistop = ">".$ncbistop;}
			unless($firststart){
				print NCBI $ncbistart."\t";
				$firststart = 1;
			}
			$laststop = $ncbistop;
		}else{
			my $ncbistart = $genes{$gene}{ex1}{start};
			my $ncbistop = $genes{$gene}{'ex'.$genes{$gene}{nbexon}}{stop};
			if($genes{$gene}{partstart} == 1 || $firstgene == 0){$ncbistart = "<".$ncbistart;}
			if($genes{$gene}{partstop} == 1){$ncbistop = ">".$ncbistop;}
			unless($firststart){
				print NCBI $ncbistart."\t";
				$firststart = 1;
			}
			$laststop = $ncbistop;
		}
	}
	print NCBI $laststop."\tCDS\n";
	print NCBI "\t\t\tproduct\tpolyprotein\n";
	foreach my $gene (@geneorder)
	{
		if($geneFlags{$gene}{completeness} == -1)
		{
			next;
		}elsif($genes{$gene}{nbexon} == 1)
		{
			my $ncbistart = $genes{$gene}{ex1}{start};
			my $ncbistop = $genes{$gene}{ex1}{stop};
			if($genes{$gene}{partstart} == 1){$ncbistart = "<".$ncbistart;}
			if($genes{$gene}{partstop} == 1){$ncbistop = ">".$ncbistop;}
			print NCBI $ncbistart."\t".$ncbistop."\tmat_peptide\n";
		}else{
			my $ncbistart = $genes{$gene}{ex1}{start};
			my $ncbistop = $genes{$gene}{'ex'.$genes{$gene}{nbexon}}{stop};
			if($genes{$gene}{partstart} == 1){$ncbistart = "<".$ncbistart;}
			if($genes{$gene}{partstop} == 1){$ncbistop = ">".$ncbistop;}
			print NCBI $ncbistart."\t".$ncbistop."\tmat_peptide\n";
		}
		print NCBI "\t\t\tproduct\t".$genefeat{$gene}{product}."\n";
	}
}else{
	foreach my $gene (@geneorder)
	{
		if($geneFlags{$gene}{completeness} == -1)
		{
			next;
		}elsif($genes{$gene}{nbexon} == 1)
		{
			my $ncbistart = $genes{$gene}{ex1}{start};
			my $ncbistop = $genes{$gene}{ex1}{stop};
			if($genes{$gene}{partstart} == 1){$ncbistart = "<".$ncbistart;}
			if($genes{$gene}{partstop} == 1){$ncbistop = ">".$ncbistop;}
			print NCBI $ncbistart."\t".$ncbistop."\tgene\n";
			print NCBI "\t\t\tgene\t".$gene."\n";
			print NCBI $ncbistart."\t".$ncbistop."\tCDS\n";
			print NCBI "\t\t\tproduct\t".$genefeat{$gene}{product}."\n";
		}else{
			my $ncbistart = $genes{$gene}{ex1}{start};
			my $ncbistop = $genes{$gene}{'ex'.$genes{$gene}{nbexon}}{stop};
			if($genes{$gene}{partstart} == 1){$ncbistart = "<".$ncbistart;}
			if($genes{$gene}{partstop} == 1){$ncbistop = ">".$ncbistop;}
			print NCBI $ncbistart."\t".$ncbistop."\tgene\n";
			print NCBI "\t\t\tgene\t".$gene."\n";

			foreach my $exon (sort keys %{$genes{$gene}})
			{
				unless($exon =~ /^ex[0-9]+/){next;}
				my $ncbistart = $genes{$gene}{$exon}{start};
				my $ncbistop = $genes{$gene}{$exon}{stop};
				if($exon eq 'ex1' && $genes{$gene}{partstart} == 1){$ncbistart = "<".$ncbistart;}
				elsif($exon eq 'ex'.$genes{$gene}{nbexon} && $genes{$gene}{partstop} == 1){$ncbistop = ">".$ncbistop;}
				if($exon eq 'ex1')
				{
					print NCBI $ncbistart."\t".$ncbistop."\tCDS\n";
				}else{
					print NCBI $ncbistart."\t".$ncbistop."\n";
				}
			}
			print NCBI "\t\t\tproduct\t".$genefeat{$gene}{product}."\n";
		}
	}
}

if($geneFlags{$geneorder[scalar(@geneorder) - 1]}{completeness} == 1 && $genes{$geneorder[scalar(@geneorder) - 1]}{'ex'.$genes{$geneorder[scalar(@geneorder) - 1]}{nbexon}}{stop} < length($consseq))
{
	print NCBI "".($genes{$geneorder[scalar(@geneorder) - 1]}{'ex'.$genes{$geneorder[scalar(@geneorder) - 1]}{nbexon}}{stop} + 1)."\t".length($consseq)."\t3'UTR\n";
	print NCBI "\t\t\tnote\tindels in UTR have not been validated\n";
}




################### SUBS ###################

sub checkGeneorder
{
	my $goodorder = 1;

	foreach my $gene (keys %genes)
	{
		if($genes{$gene}{flags} == 0){next;}
		foreach my $exon (keys %{$genes{$gene}})
		{
			unless($exon =~ /^ex/){next};
			foreach my $compgene (keys %genes)
			{
				if($genes{$compgene}{flags} == 0){next;}
				foreach my $compexon (keys %{$genes{$compgene}})
				{
					unless($compexon =~ /^ex/){next};
					if($gene eq $compgene && $exon eq $compexon){next;}
#					print $gene."\t".$exon."\t".$compgene."\t".$compexon."\n";
					if(compValue($genes{$gene}{$exon}{start}, $genes{$compgene}{$compexon}{start}) != $refGeneorder{$gene}{$exon}{$compgene}{$compexon}{startvstart})
					{
						$goodorder = 0;
					}
					if(compValue($genes{$gene}{$exon}{start}, $genes{$compgene}{$compexon}{stop}) != $refGeneorder{$gene}{$exon}{$compgene}{$compexon}{startvstop})
					{
						$goodorder = 0;
					}
					if(compValue($genes{$gene}{$exon}{stop}, $genes{$compgene}{$compexon}{start}) != $refGeneorder{$gene}{$exon}{$compgene}{$compexon}{stopvstart})
					{
						$goodorder = 0;
					}
					if(compValue($genes{$gene}{$exon}{stop}, $genes{$compgene}{$compexon}{stop}) != $refGeneorder{$gene}{$exon}{$compgene}{$compexon}{stopvstop})
					{
						$goodorder = 0;
					}
				}
			}
		}
	}
	return $goodorder;
}

sub checkGeneFlags
{
	my $gene = shift;
	my $completegene = shift;

	my $allgood = 1;

	# Check if size is within 10% of original gene
	$geneFlags{$gene}{size} = checkGeneSize($genelist{$gene}{seq}, $finalGeneStrings{$gene}{aa});
	if($geneFlags{$gene}{size} == 0){$allgood = 0;}

	$geneFlags{$gene}{Nstr} = checkNString($gene);
	if($geneFlags{$gene}{Nstr} == 0){$allgood = 0;}

	# Check if gene has start and stop codon (if required)
	$geneFlags{$gene}{startstop} = checkStartStop($gene, $finalGeneStrings{$gene}{aa});
	if($geneFlags{$gene}{startstop} == 0){$allgood = 0;}

	# Check if gene has full coverage
	$geneFlags{$gene}{completeness} = $completegene;
	if($geneFlags{$gene}{completeness} != 1){$allgood = 0;}

	# Check if gene has good number of exons
	if($genes{$gene}{nbexon} == $genelist{$gene}{nbexon})
	{
		$geneFlags{$gene}{exons} = 1;
	}else{
		$geneFlags{$gene}{exons} = 0;
	}
	if($geneFlags{$gene}{exons} == 0){$allgood = 0;}

	# Check if gene has frameshift; 1 = no, 0 = yes
	$geneFlags{$gene}{frameshift} = abs(1-$genes{$gene}{frameshifted});
	if($geneFlags{$gene}{frameshift} == 0){$allgood = 0;}

	return $allgood;
}


sub readGene
{
	my $genefile = shift;
	my $gene = shift;
	my $expectedExons = shift;

	open(GENEOUT, $genefile);
	my @exons;
	my $nbExon = 0;

	my $nbGenes = 0;
	while(my $line = <GENEOUT>)
	{
		chomp $line;
		if($line =~ /Exon (.+?) (.+?) .+? (.+)/)
		{
			$exons[$nbExon]{start} = $1;
			$exons[$nbExon]{stop} = $2;
			$nbExon++;
		}elsif($line =~ /^Gene ([0-9]+)$/)
		{
			$nbGenes = $1;
		}
	}
	close GENEOUT;

	$genes{$gene}{frameshifted} = 0;

	if($nbGenes > 1){
		$genes{$gene}{frameshifted} = 1;
	}

	if($nbExon == 1 && $expectedExons == 1)
	{
		$genes{$gene}{ex1}{start} = $exons[0]{start};
		$genes{$gene}{ex1}{stop} = $exons[0]{stop};
	}else{
		my $prevstop = 0;
		my @longExons;
		my $longExonFlag = 0;

		foreach my $curexon (@exons)
		{
			if($prevstop){
				if($$curexon{start} <= ($prevstop + $option{maxannogaplen}))
				{
					$longExons[$longExonFlag]{stop} = $$curexon{stop};
				}else{
					$longExonFlag++;
					$longExons[$longExonFlag]{start} = $$curexon{start};
					$longExons[$longExonFlag]{stop} = $$curexon{stop};
				}
			}else{
				$longExons[$longExonFlag]{start} = $$curexon{start};
				$longExons[$longExonFlag]{stop} = $$curexon{stop};
			}
			$prevstop = $$curexon{stop};
		}

		if(scalar(@longExons) == 1)
		{
			$genes{$gene}{ex1}{start} = $longExons[0]{start};
			$genes{$gene}{ex1}{stop} = $longExons[0]{stop};
		}else{
			$longExonFlag = 0;
			foreach my $longExon (@longExons)
			{
				$genes{$gene}{'ex'.($longExonFlag+1)}{start} = $longExons[$longExonFlag]{start};
				$genes{$gene}{'ex'.($longExonFlag+1)}{stop} = $longExons[$longExonFlag]{stop};
				$longExonFlag++;
			}
		}
		$nbExon = scalar(@longExons);
	}

	$genes{$gene}{nbexon} = $nbExon;
}

sub readGenewise
{
	my $inputfile = shift;
	my $gene = shift;

	my $minpctid = $option{minpctid};
	my $minpctlen = $option{minpctlen};
	my $diffForCompletion = 10;

	open(GENEWISE, $inputfile);

	my $firstpos = 0;
	my $lastpos = 0;
	my $fullgene = 1;
	my $matchstr = '';
	my $firstnt = $genes{$gene}{ex1}{start};
	my $lastnt = $genes{$gene}{'ex'.$genes{$gene}{nbexon}}{stop};

	$genes{$gene}{partstart} = 0;
	$genes{$gene}{partstop} = 0;

	while(my $line = <GENEWISE>)
	{
		chomp $line;
		if($line =~ /^($gene\s+)(.+?)\s(.+)/)
		{
#Temp!!!
#			$line = <GENEWISE>;

			chomp $line;
			my $pos = $2;
			my $aastr = $3;
			$aastr =~ s/\s+$//;
			if($aastr =~ /^\s+/)
			{
				$pos++;
				$aastr =~ s/^\s+//;
			}
#Temp
#			print $aastr."\n";
			unless($firstpos){
				$firstpos = $pos;
			}
			$aastr =~ s/-//g;
			$lastpos = $pos + length($aastr) - 1;
#Temp
#			print $lastpos."\n";
			$line = <GENEWISE>;

			chomp $line;
			$line =~ /\s+(.+?)\s+$/;
			$matchstr .= $1;
#Temp
#			print length($matchstr)."\n";
#			print $matchstr."\n";
			$line = <GENEWISE>;
		}
	}
#Temp!!
#	print $gene."\n";
#	print $matchstr."\n";
	### Add one position if gene has stop codon
	if($lastpos == ($genelist{$gene}{size} - 1) && $genefeat{$gene}{stop} == 1 && $geneticcode{substr($consseq, $lastnt, 3)} eq '*')
	{
		$lastpos++;
		$lastnt += 3;
		$genes{$gene}{'ex'.$genes{$gene}{nbexon}}{stop} += 3;
	}


	### Check if you can fill the ends in case the nt-based alignment disagree
#	print "alignStart : ".$alignGenelist{$gene}{'ex1'}{start}."\t$firstnt\n";

	if($firstpos > 1 && $firstpos < $diffForCompletion && $alignGenelist{$gene}{'ex1'}{start} < $firstnt)
	{
#		print "TRYING TO EXTEND START ".$gene."\t".$firstpos."\n";

		my $assemstr = '';
		my $refstr = '';
		my $asspos = 0;
		my $posdiff = 0;
		my $aaposdiff = 0;
		for(my $alnpos = 0; $alnpos < length($alignass); $alnpos++)
		{
			my $curass = substr($alignass, $alnpos, 1);
			my $curref = substr($alignref, $alnpos, 1);

			if($curass ne '-'){$asspos++;}
			if($asspos >= $alignGenelist{$gene}{'ex1'}{start} && $asspos <= $firstnt)
			{
				$assemstr .= $curass;
				$refstr .= $curref;
				if($curass eq '-'){
					$posdiff++;
				}
				if($curref eq '-')
				{
					$posdiff--;
				}
			}
		}
		$aaposdiff = $posdiff/3;

#		print $assemstr."\n";
#		print $refstr."\n";
#		print $posdiff."\n";
#		print $aaposdiff."\n";

		my $tmpfirstpos = $firstpos;
		my $tmpfirstnt = $firstnt;
		while($tmpfirstpos > 1 && $tmpfirstnt > 0)
		{
			$tmpfirstpos -= 1;
			$tmpfirstnt -= 3;
#			print $tmpfirstnt."\t".$geneticcode{substr($consseq, $tmpfirstnt - 1, 3)}."\n";
			if($tmpfirstnt == $alignGenelist{$gene}{'ex1'}{start})
			{
#				print "FIXED\n";
				$tmpfirstnt -= $posdiff;
				$tmpfirstpos -= $aaposdiff;
				$firstnt = $tmpfirstnt;
				$firstpos = $tmpfirstpos;
				$genes{$gene}{'ex1'}{start} = $firstnt;
				last;
			}
		}
	}

	if($lastpos < $genelist{$gene}{size} && $lastpos > ($genelist{$gene}{size} - $diffForCompletion) && $alignGenelist{$gene}{'ex'.$genes{$gene}{nbexon}}{stop} > $lastnt)
	{
#		print "TRYING TO EXTEND ".$gene."\t".$lastpos."\n";

		my $assemstr = '';
		my $refstr = '';
		my $asspos = 0;
		my $posdiff = 0;
		my $aaposdiff = 0;
		for(my $alnpos = 0; $alnpos < length($alignass); $alnpos++)
		{
			my $curass = substr($alignass, $alnpos, 1);
			my $curref = substr($alignref, $alnpos, 1);

			if($curass ne '-'){$asspos++;}
			if($asspos >= $lastnt && $asspos <= $alignGenelist{$gene}{'ex'.$genes{$gene}{nbexon}}{stop})
			{
				$assemstr .= $curass;
				$refstr .= $curref;
				if($curass eq '-'){
					$posdiff++;
				}
				if($curref eq '-')
				{
					$posdiff--;
				}
			}
		}
		$aaposdiff = $posdiff/3;

		my $tmplastpos = $lastpos;
		my $tmplastnt = $lastnt;
		while($tmplastpos < $genelist{$gene}{size} && $tmplastnt < length($consseq))
		{
			$tmplastpos++;
			$tmplastnt += 3;
#			print $tmplastpos."\t".$genelist{$gene}{size}."\t".$geneticcode{substr($consseq, $tmplastnt - 3, 3)}."\n";
			if($geneticcode{substr($consseq, $tmplastnt - 3, 3)} eq '*')
			{
				last;
			}
		}
#		print $tmplastnt."\t".$alignGenelist{$gene}{'ex'.$genes{$gene}{nbexon}}{stop}."\n";
		if($tmplastnt == $alignGenelist{$gene}{'ex'.$genes{$gene}{nbexon}}{stop})
		{
#			print "FIXED\n";
			$tmplastpos -= $aaposdiff;
			$lastnt = $tmplastnt;
			$lastpos = $tmplastpos;
			$genes{$gene}{'ex'.$genes{$gene}{nbexon}}{stop} = $lastnt;
		}
	}

	for(my $i = 1; $i <= $genes{$gene}{nbexon}; $i++)
	{
		$finalGeneStrings{$gene}{nt} .= substr($consseq, $genes{$gene}{'ex'.$i}{start}-1, $genes{$gene}{'ex'.$i}{stop}-$genes{$gene}{'ex'.$i}{start}+1);
	}


	$finalGeneStrings{$gene}{aa} = translateDnaString($finalGeneStrings{$gene}{nt});
#	print "\t".$firstnt."\t".$lastnt."\t".$firstpos."\t".$lastpos."\t".$genelist{$gene}{size}."\n";#.$genefeat{$gene}{stop}."\t".substr($consseq, $lastnt, 3)."\n";

	if($firstpos == 1 && $lastpos == $genelist{$gene}{size}){
		return $fullgene;
	}

	### Gene is partial or missing
	$fullgene = 0;
	my $newmatchstr = $matchstr;
	$newmatchstr =~ s/ //g;
	$newmatchstr =~ s/\+//g;
	my $pctid = length($newmatchstr)/length($matchstr);
	my $pctlen = ($lastpos - $firstpos + 1)/$genelist{$gene}{size};

	if($pctid < $minpctid || $pctlen < $minpctlen)
	{
		$pctid = int($pctid*10000)/100;
		$pctlen = int($pctlen*10000)/100;
#		print "$gene is Missing (genewise %id = $pctid and %length = $pctlen)\n";
		$fullgene = -1;
	}else{
		$pctid = int($pctid*10000)/100;
		$pctlen = int($pctlen*10000)/100;
#		print "$gene is Partial (positions $firstpos-$lastpos instead of 1-".$genelist{$gene}{size}.", %id of $pctid)\n";
		if($firstpos > 1){$genes{$gene}{partstart} = 1;}
		if($lastpos < $genelist{$gene}{size}){$genes{$gene}{partstop} = 1;}
	}
	return $fullgene;
}





sub genelistFromAlignment
{
	my $refseq = shift;
	my $otherseq = shift;

	my %posconversion;

	my $refpos = 0;
	my $otherpos = 0;

	for(my $seqpos = 0; $seqpos < length($refseq); $seqpos++)
	{
		my $curref = substr($refseq, $seqpos, 1);
		my $curother = substr($otherseq, $seqpos, 1);

		if($curref eq "-")
		{
			$otherpos++;
		}else{
			if($curother ne "-")
			{
				$otherpos++;
				$refpos++;
				$posconversion{$refpos} = $otherpos;
			}else{
				$refpos++;
				$posconversion{$refpos} = $otherpos."+";
			}
		}
	}

	my $lastotherpos = $otherpos;

#	open(OUTPUT, ">$output");

	foreach my $gene (keys %genelist)
	{
		foreach my $exon (keys %{$genelist{$gene}})
		{
			unless($exon =~ /^ex.+/){next;}
			my $start = $posconversion{$genelist{$gene}{$exon}{start}};
			my $end = $posconversion{$genelist{$gene}{$exon}{stop}};

			my $exact = 1;
			my $startplus = 0;
			if($start =~ /(.+)\+/)
			{
				$start = $1;
				$startplus = 1;
				$exact = 0;
			}
			my $endplus = 0;
#Temp!!!
#			print $gene."\t".$exon."\t".$end."\n";
			if($end =~ /(.+)\+/)
			{
				$end = $1;
				$endplus = 1;
				$exact = 0;
			}

			my $full = 1;

			if($start == 0)
			{
				$start = "BeforeStart";
				$full = 0;
			}elsif($start eq $lastotherpos && $startplus)
			{
				$start = "AfterEnd";
				$full = 0;
			}

			if($end == 0)
			{
				$end = "BeforeStart";
				$full = 0;
			}elsif($end eq $lastotherpos && $endplus)
			{
				$end = "AfterEnd";
					$full = 0;
			}

			$alignGenelist{$gene}{$exon}{start} = $start;
			$alignGenelist{$gene}{$exon}{stop} = $end;
			$alignGenelist{$gene}{$exon}{full} = $full;
			$alignGenelist{$gene}{$exon}{exact} = $exact;
#			print "AlignGeneList\t$gene\t$exon\tStart\t".$alignGenelist{$gene}{$exon}{start}."\n";
#			print "AlignGeneList\t$gene\t$exon\tStop\t".$alignGenelist{$gene}{$exon}{stop}."\n";
		}
	}
}




sub checkGeneSize
{
	my $refaa = shift;
	my $assaa = shift;

	my $lengthFactor = 0.1;

	if(abs(length($refaa) - length($assaa))/max(length($refaa), length($assaa)) <= $lengthFactor)
	{
		return 1;
	}
	return 0;
}

sub checkNString
{
	my $gene = shift;
	for(my $exon = 1; $exon <= $genes{$gene}{nbexon}; $exon++)
	{
		my $dnastring = substr($consseq, $genes{$gene}{'ex'.$exon}{start} - 1, $genes{$gene}{'ex'.$exon}{stop} - $genes{$gene}{'ex'.$exon}{start} + 1);

		for(my $pos = 0; $pos < length($dnastring) - 1; $pos++)
		{
			if(substr($dnastring, $pos, 1) eq 'N' && substr($dnastring, $pos + 1, 1) eq 'N')
			{
				return 0;
			}
		}
	}
	return 1;
}

sub checkStartStop
{
	my $gene = shift;
	my $aaseq = shift;

	# Check if starts with M
	if($genefeat{$gene}{start} == 1 && substr($aaseq, 0, 1) ne 'M')
	{
		return 0;
	}

	# Check if ends with stop codon
	if($genefeat{$gene}{stop} == 1 && substr($aaseq, length($aaseq) - 1, 1) ne '*')
	{
		return 0;
	}

	return 1;
}


sub alignPair
{
	my $seq1 = shift;
	my $seq2 = shift;
	my $alignPairOut = shift;

	open(MFA, ">$alignPairOut.mfa");
	print MFA ">Ref\n$seq1\n>Assembly\n$seq2\n";
	close MFA;

#	system($musclepath."muscle -in $alignPairOut.mfa -out $alignPairOut.afa  -quiet");
#	system($clustalpath."clustalw2 -INFILE=$alignPairOut.mfa -OUTFILE=$alignPairOut.afa -ALIGN -TYPE=DNA -OUTPUT=FASTA");
	system($mafftpath."mafft --thread 4 --clustalout --auto $alignPairOut.mfa > $alignPairOut.afa");
	
	(my $alignseq1, my $alignseq2) = readAlign("$alignPairOut.afa");
	return $alignseq1, $alignseq2;
}

sub readAlign
{
	my $afa = shift;
	open(AFA, $afa);
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

sub max
{
	my $nb1 = shift;
	my $nb2 = shift;
	if($nb1 >= $nb2){return $nb1;}
	return $nb2;
}

sub compValue
{
	my $val1 = shift;
	my $val2 = shift;

	if($val1 > $val2){return 1;}
	if($val1 == $val2){return 0;}
	if($val1 < $val2){return -1;}
}
