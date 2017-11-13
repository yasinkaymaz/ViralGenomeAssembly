#!/usr/bin/perl -w

## estMOI: estimating multiplicity of infection using parasite deep sequencing data
## Version 1.03
## contact: samuel.assefa at lshtm.ac.uk OR mark.preston at lshtm.ac.uk

use strict;
use Carp;
use Getopt::Long;
use POSIX qw/strftime/;
my $version =1.03;
my %commands;
my %options;
my %files;

$commands{zless}    = "zless";
$commands{grep}     = "grep";
$commands{awk}      = "awk";
$commands{head}     = "head";
$commands{samtools} = "samtools";

$options{readl}   = 76;
$options{maxsnp}  = 3;
$options{mincov}  = 10;
$options{minhap}  = 3;
$options{maxfact} = 90;
$options{mindis}  = 10;
$options{maxdis}  = 500;
$options{flank}   = 500;
$options{snpq}    = 30;
$options{exclude} = "";
$options{out}     = "estMoi";
$options{tmp}     = 0;
$options{debug}   = 0;

if (@ARGV < 3) {
  print STDERR "Usage: $0 <bam> <vcf.gz> <fasta> [OPTIONS]
    \t--readl=int\tMinimum read length [default $options{readl} ]
    \t--maxsnp=int\tMaximum number of SNPs [default $options{mincov} ]
    \t--mincov=int\tMininum number of reads over a SNP [default $options{mincov}]
    \t--minhap=int\tMinimum frequency of haplotypes to consider [default $options{minhap}]
    \t--maxfact=int\tPercentile cutoff for adjusting genomewide MOI estimate [default $options{maxfact}]
    \t--mindis=int\tMinimum distance between any two SNPs [default $options{mindis}]
    \t--maxdis=int\tMaximum distance between the first and last SNP [default $options{maxdis}]
    \t--flank=int\tFlanking size for excluded regions [default $options{flank}] 
    \t--out=string\tOutput file prefix [default $options{out}]
    \t--debug \tDebug output [default $options{debug}]
    Note:
    \tFor paired end libraries, maxdis should not exceed the insert/fragment size of library
    \t--minhap is used to disregard haplotypes that are less frequent than minhap\n
    \tOptimisation on your data could be done by changing --minhap and --maxfact parameters\n\n";
  exit;
}

GetOptions(
           "readl=i"      => \$options{readl},
           "maxsnp=i"     => \$options{maxsnp},
           "mincov=i"     => \$options{mincov},
           "minhap=i"     => \$options{minhap},
           "maxfact"      => \$options{maxfact},
           "maxdis=i"     => \$options{maxdis},
           "exclude=s"    => \$options{exclude},
           "out=s"        => \$options{out},
           "tmp"          => \$options{tmp},
           "debug"        => \$options{debug},
           "rerun=i"      => \$options{rerun});

($files{bam}, $files{vcf}, $files{ref}) = @ARGV;
$files{fai} = $files{ref} . ".fai";
$options{out} = $files{bam} if ($options{out} eq "");

die("Can't find bam file: $files{bam}\n")             unless (-e $files{bam});
die("Can't find vcf.gz file: $files{vcf}\n")          unless (-e $files{vcf});
die("Can't find reference index file: $files{fai}\n") unless (-e $files{fai});

#  Preprocessing
my ($r1, $r2) = hash_fai();
my %lhChromosomes = %{$r1};
my @laChromosomes = @{$r2};

my $dir        = $options{out} . ".moi";
my $base       = $options{out} . ".moi.$options{maxsnp}.$options{mindis}.$options{maxdis}.$options{minhap}";
$files{log}    = $base . ".log";
$files{output} = $base . ".txt";
my $giMaxfact  = $options{maxfact}; #defining a global variable for adjusting the maximun MOI per chromosome.
print STDERR "#\tRUNNING estMOI version $version\n";
mkdir $dir;
debug("## Started");

open(OUT, ">$files{output}") or die ("Unable to open output file: $files{output}\n");
open(LOG, ">$files{log}") or die ("Unable to open log file: $files{log}\n");

my @MOI;
my %finalMOI;
print STDERR "#\t";
foreach my $chr (@laChromosomes) {
  print STDERR ".";
  debug("# Working on Chromosome. $chr");
  $files{snps}   = "$dir/$chr.snps";
  $files{sam}    = "$dir/$chr.sam";
  $files{pileup} = "$dir/$chr.pileup";

  debug("#\tGet samfiles for $chr");
  `$commands{samtools} view -f 0x0002 -q 10 $files{bam} $chr | $commands{awk} '\$6!~/[IDS]/' > $files{sam}` unless (-e $files{sam});

  debug("#\tGet filtered SNPs for $chr");
  `$commands{zless} $files{vcf} | $commands{grep} $chr | $commands{awk} '{ if (\$6 >= $options{snpq}) { print } }' | $commands{grep} -v INDEL > $files{snps}` unless (-e $files{snps} && -s $files{snps} > 0);

  debug("#\tGet pileup file for $chr");
  samToPileup($chr) unless (-e $files{pileup});

  debug("#\tEstimating MOI using paired-reads");
  my ($refMoi, $log) = pairedMOI($chr);
  push @MOI, @{$refMoi};
  print LOG "$log";
  if($options{tmp} == 0) {`rm -f $files{snps} $files{sam} $files{pileup}`;}
}
print STDERR "\n#\tPRINT ESTIMATES...\n";
my $liSize = scalar(@MOI);
my %lhCount;
map { $lhCount{$_}++ } @MOI;
my $total=0;
my $foundMOI=0;
my $moi=1;

print OUT "#MOI\tCount\t%Total\n";
print STDERR "#\tMOI\tCount\t%Total\n";
foreach (sort(keys(%lhCount)))
{
  $total+=$lhCount{$_};
  my $prop= sprintf ("%.2f", ($total/$liSize)*100);
  if ($prop >=$giMaxfact && $foundMOI==0) {
    $moi=$_; $foundMOI=1;
    print OUT "$_\t$lhCount{$_}\t$prop\tMOI-estimate\n";
    print STDERR "#\t$_\t$lhCount{$_}\t$prop\tMOI-estimate\n";
  }
  else{
    print OUT "$_\t$lhCount{$_}\t$prop\n";
    print STDERR "#\t$_\t$lhCount{$_}\t$prop\n";
    }
}

print STDERR "\t#DONE MOI-estimate using $files{bam}\n";
close OUT;
close LOG;
debug("Estimates_of_MOI $files{bam}  at a cutoff of $giMaxfact is $moi");
debug("DONE MOI-estimate using $files{bam}");
if($options{tmp} == 0) {`rm -rf $dir`;}

exit;

##############################################################################
#
#
#
sub samToPileup {
  my $psChromosome = shift;
  my $liExclusions = 0;
  my @laStart;
  my @laStop;
  my %lhSnps;
  my %lhSam;
  

  if ($options{exclude} ne "") {
    debug("## Hash-Excluded-regions");
    if (-e $options{exclude}) {
      open(EXCLUDE, $options{exclude}) or die("Cannot open exclude file: $options{exclude}\n");
      while (<EXCLUDE>) {
        my ($lsChromosome, $start, $stop) = split /\s+/, $_;
        if ($lsChromosome eq $psChromosome) {
           $liExclusions++;
           push @laStart, $start - $options{flank};
           push @laStop, $stop + $options{flank};
        }
      }
      close EXCLUDE;
      debug("## Hash-Excluded-regions:done");
    } else {
      print STDERR "WARNING: couldn't open exclude-file $options{exclude}...use absolute file name? \n";
    }
  }

  open(SNPS, $files{snps}) or die("Unable to open snps file: $files{snps}");
  while (<SNPS>) {
    my $pos = (split /\s+/, $_)[1];
    for (my $i = $pos - $options{readl}; $i < $pos + $options{readl}; $i++) {
      my $lbKeep = 1;
      for (my $j = 0; $j < $liExclusions; $j++) {
        if ($laStart[$j] <= $i && $i <= $laStop[$j]) {
          $lbKeep = 0;
          last;
        }
      }
      $lhSnps{$i} = 1 if ($lbKeep);
    }
  }
  close SNPS;

  open SAM, $files{sam} or die("Unable to open sam file: $files{sam}");
  while (<SAM>) {
    my @laLine = split /\s+/, $_;
    next if (@laLine < $options{mincov});
    my $start = $laLine[3];
    next if (! defined $lhSnps{$start});
    my @bases = split //, $laLine[9];
    for (my $i = 0; $i <= $#bases; $i++) {
      if (exists($lhSam{$start + $i})) {
        $lhSam{$start + $i} .= "%$laLine[0],$bases[$i]";
      } else {
        $lhSam{$start + $i} = "$laLine[0],$bases[$i]";
      }
    }
  }
  close SAM;

  open PILEUP, ">$files{pileup}" or die("Unable to open pileup file for writing: $files{pileup}");
  foreach (keys %lhSam) {
    print PILEUP "$_\t" . $lhSam{$_} . "\n";
  }
  close PILEUP
}

sub pairedMOI {
  my ($chr) = @_;
  my (@snp_Positions, %counts, %snps, %sam_positions, %lhLocalPil, @moiArray);
  my $chrMoi = 1;
  my $log    = "";

  open(SNPS, $files{snps}) or die ("Unable to open snps file: $files{snps}\n");
  while (<SNPS>) {
    chomp;
    my @sp = split /\s+/, $_;
    $snps{$sp[1]} = $_;
    push @snp_Positions, $sp[1];
  }
  close SNPS;

  open PILEUP, $files{pileup} or die("Unable to open pileup file: $files{pileup}");
  while (<PILEUP>) {
    chomp;
    my ($pos, $line) = split /\s+/, $_;
    $sam_positions{$pos} = $line;
  }
  close PILEUP;

  my $excluded = 0;
  my $countS   = scalar(@snp_Positions);
  debug("#\t", $countS, " snps found in file");
  my ($pass, $fail, $comboFail) = (0, 0, 0);

  for (my $i = 0 ; $i < ($countS - $options{maxsnp}) ; $i++) {    #take three snps.. HRD codede for three snps for now,,
    my $pos1 = $snp_Positions[$i];

    my $id1 = "$chr:$pos1";
    my @laPos;
    for (my $j = $i + 1 ; $j < ($countS - $options{maxsnp}) ; $j += 1) {
      my $pos2 = $snp_Positions[$j];
      if ($pos2 - $pos1 < $options{maxdis}) {
        push @laPos, $pos2;
      } else {
        last;
      }
    }
    my @snpCombs;
    for (my $k = 0 ; $k < @laPos ; $k += 1) {
      my $p1 = $laPos[$k];
      for (my $n = $k + 1 ; $n < @laPos ; $n += 1) {
        my $p2 = $laPos[$n];
        my @snp = ($pos1, $p1, $p2);
        push @snpCombs, \@snp;
      }
    }

    foreach (@snpCombs) {
      my @snps = @{$_};
      my @distances;
      my $liexcluded  = 0;
      my $onBothReads = 0;
      for (my $j = 1 ; $j < @snps ; $j += 1) {
        my $p = $snps[$j];
        my $d = $snps[$j] - $snps[$j - 1];
        push @distances, $d;
        if ($d < $options{mindis}) { $liexcluded += 1; }
        if ($d >= $options{readl}) {
          $onBothReads += 1;
        }
      }
      if ($liexcluded > 0 or $onBothReads eq 0) {
        $comboFail += 1;
        next;
      } else {
        $pass += 1;
      }
      my %readHapHash;
      my %haploHash;
      my @haplos;
      for (my $k = 0 ; $k < @snps ; $k += 1) {
        my $p1 = $snps[$k];
        if (!defined $sam_positions{$p1}) { next; }

        my @sp = split /%/, $sam_positions{$p1};
        if (@sp < $options{mincov}) {
          next;
        }

        foreach (@sp) {
          my ($readID, $call) = split /,/, $_;
          if (defined $readHapHash{$readID}) {
            $readHapHash{$readID} .= "%$k:$call";
          } else {
            $readHapHash{$readID} = "$k:$call";
          }

        }
      }
      foreach (keys %readHapHash) {
        my @sp = split /%/, $readHapHash{$_};
        my %tmpH;

        if (@sp == $options{maxsnp}) {
          foreach (@sp) {
            my ($posi, $call) = split /:/, $_;
            $tmpH{$posi} = $call;
          }
        }
        my $hap = "";
        my $inc = 0;
        for (my $i = 0 ; $i < $options{maxsnp} ; $i += 1) {
          my $posi = $_;
          if (defined $tmpH{$i}) {
            $hap .= " " . $tmpH{$i};
          } else {
            $hap .= " -";
            $inc = 1;
          }
        }
        if ($inc == 0) {
          if (defined $haploHash{$hap}) {
            $haploHash{$hap} += 1;
          } else {
            $haploHash{$hap} = 1;
            push @haplos, $hap;
          }
        }
      }
      my $numHap = 0;
      my $res    = "";
      if (@haplos > 0) {

        foreach (@haplos) {
          my $haploFrequency = $haploHash{$_};
          if ($haploFrequency >= $options{minhap}) {
#            $res .= "\t# $files{bam} $chr @snps Hapotype:\t$_\t$haploHash{$_}\n";
            $res .= "\t# $files{bam} $chr @snps Haplotype:\t$_\t$haploHash{$_}\t$haploFrequency\n";

            $numHap += 1;
          }
        }
      }
      if ($numHap > 0) {
        push @moiArray, $numHap;
        $log .= "$files{bam}\t$chr\t@snps\t$numHap\n$res";
        debug("$files{bam}\t$chr\t@snps\t$numHap\n$res");
      }
    }
  }
  debug("#PASSedSNPs=$pass\tFAILEDsnps=$fail\tFAILedSNPCombinations=$comboFail");

  return (\@moiArray, $log);
}

sub hash_fai {
  open(F, $files{fai}) or die "Cant open fai file: $files{fai}\n";
  my %rhLengths;
  my @raChromosomes;
  while (<F>) {
    chomp;
    my @laLine = split /\s+/, $_;
    my $lsChromosome = $laLine[0];
    push @raChromosomes, $lsChromosome;
    $rhLengths{$lsChromosome} = $laLine[1];
  }
  close(F);
  return (\%rhLengths, \@raChromosomes);
}

sub trim {
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}

sub debug {
  my $lsMessage = shift;
  if ($options{"debug"}) {
    print STDERR strftime('%F %T', localtime) . ": $lsMessage\n";
  }
}

