#! perl -w

# Copyright © 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.
#
# This code was developed by Patrick Charlebois <patrickc@broadinstitute.org>

use strict;

my $configfile = shift || die("Usage : perl configPaths.pl configfile");

open(CONFIG, $configfile) || die("Unable to open $configfile");

my $scriptpath = '';
my $musclepath = '';
my $mosaikpath = '';
my $mosaiknetworkpath = '';
my $samtoolspath = '';
my $perlpath = '';
my $Rpath = '';
my $refDataPath = '';
my $genewisepath = '';
my $genewisecfgpath = '';

while(my $line = <CONFIG>)
{
	if($line =~ /scriptpath \= \'(.*)\'.*/)
	{
		$scriptpath = $1;
		if($scriptpath && substr($scriptpath, length($scriptpath) - 1) ne '/'){$scriptpath .= '/';}
	}elsif($line =~ /musclepath \= \'(.*)\'.*/)
	{
		$musclepath = $1;
		if($musclepath && substr($musclepath, length($musclepath) - 1) ne '/'){$musclepath .= '/';}
	}elsif($line =~ /samtoolspath \= \'(.*)\'.*/)
	{
		$samtoolspath = $1;
		if($samtoolspath && substr($samtoolspath, length($samtoolspath) - 1) ne '/'){$samtoolspath .= '/';}
	}elsif($line =~ /mosaikpath \= \'(.*)\'.*/)
	{
		$mosaikpath = $1;
		if($mosaikpath && substr($mosaikpath, length($mosaikpath) - 1) ne '/'){$mosaikpath .= '/';}
	}elsif($line =~ /mosaiknetworkpath \= \'(.*)\'.*/)
	{
		$mosaiknetworkpath = $1;
		if($mosaiknetworkpath && substr($mosaiknetworkpath, length($mosaiknetworkpath) - 1) ne '/'){$mosaiknetworkpath .= '/';}
	}elsif($line =~ /perlpath \= \'(.*)\'.*/)
	{
		$perlpath = $1;
	}elsif($line =~ /Rpath \= \'(.*)\'.*/)
	{
		$Rpath = $1;
		if($Rpath && substr($Rpath, length($Rpath) - 1) ne '/'){$Rpath .= '/';}
	}elsif($line =~ /genewisepath \= \'(.*)\'.*/)
	{
		$genewisepath = $1;
		if($genewisepath && substr($genewisepath, length($genewisepath) - 1) ne '/'){$genewisepath .= '/';}
	}elsif($line =~ /genewisecfgpath \= \'(.*)\'.*/)
	{
		$genewisecfgpath = $1;
		if($genewisecfgpath && substr($genewisecfgpath, length($genewisecfgpath) - 1) ne '/'){$genewisecfgpath .= '/';}
	}elsif($line =~ /refDataPath \= \'(.*)\'.*/)
	{
		$refDataPath = $1;
		if($refDataPath && substr($refDataPath, length($refDataPath) - 1) ne '/'){$refDataPath .= '/';}
	}
}

close(CONFIG);

open(ANNOT, "$scriptpath"."annotate.pl");
open(ANNOTTMP, ">$scriptpath"."annotatetmp.pl");

my $linecount = 0;
while(my $line = <ANNOT>)
{
	if($linecount == 0 && $line =~ /^\#\!.+/)
	{
		print ANNOTTMP "#!".$perlpath." perl -w\n";
	}elsif($line =~ /my \$genewisepath \= .+/)
	{
		print ANNOTTMP 'my $genewisepath = "'.$genewisepath."\";\n";
	}elsif($line =~ /my \$genewisecfgpath \= .+/)
	{
		print ANNOTTMP 'my $genewisecfgpath = "'.$genewisecfgpath."\";\n";
	}elsif($line =~ /my \$musclepath \= .+/)
	{
		print ANNOTTMP 'my $musclepath = "'.$musclepath."\";\n";
	}else{
		print ANNOTTMP $line;
	}
	$linecount++;
}

close(ANNOTTMP);
close(ANNOT);

open(RUNMOSAIK, "$scriptpath"."runMosaik2.pl");
open(RUNMOSAIKTMP, ">$scriptpath"."runMosaiktmp.pl");

$linecount = 0;
while(my $line = <RUNMOSAIK>)
{
	if($linecount == 0 && $line =~ /^\#\!.+/)
	{
		print RUNMOSAIKTMP "#!".$perlpath." perl -w\n";
	}elsif($line =~ /my \$mosaiknetworkpath \= .+/)
	{
		print RUNMOSAIKTMP 'my $mosaiknetworkpath = "'.$mosaiknetworkpath."\";\n";
	}elsif($line =~ /my \$mosaikpath \= .+/)
	{
		print RUNMOSAIKTMP 'my $mosaikpath = "'.$mosaikpath."\";\n";
	}elsif($line =~ /my \$samtoolspath \= .+/)
	{
		print RUNMOSAIKTMP 'my $samtoolspath = "'.$samtoolspath."\";\n";
	}else{
		print RUNMOSAIKTMP $line;
	}
	$linecount++;
}

close(RUNMOSAIKTMP);
close(RUNMOSAIK);

open(MERGER, "$scriptpath"."contigMerger.pl");
open(MERGERTMP, ">$scriptpath"."contigMergertmp.pl");

$linecount = 0;
while(my $line = <MERGER>)
{
	if($linecount == 0 && $line =~ /^\#\!.+/)
	{
		print MERGERTMP "#!".$perlpath." perl -w\n";
	}elsif($line =~ /my \$musclepath \= .+/)
	{
		print MERGERTMP 'my $musclepath = "'.$musclepath."\";\n";
	}elsif($line =~ /my \$scriptpath \= .+/)
	{
		print MERGERTMP 'my $scriptpath = "'.$scriptpath."\";\n";
	}elsif($line =~ /my \$Rpath \= .+/)
	{
		print MERGERTMP 'my $Rpath = "'.$Rpath."\";\n";
	}else{
		print MERGERTMP $line;
	}
	$linecount++;
}

close(MERGERTMP);
close(MERGER);

open(CONT2ASS, "$scriptpath"."vfat.pl");
open(CONT2ASSTMP, ">$scriptpath"."vfattmp.pl");

$linecount = 0;
while(my $line = <CONT2ASS>)
{
	if($linecount == 0 && $line =~ /^\#\!.+/)
	{
		print CONT2ASSTMP "#!".$perlpath." perl -w\n";
	}elsif($line =~ /my \$refDataPath \= .+/)
	{
		print CONT2ASSTMP 'my $refDataPath = "'.$refDataPath."\";\n";
	}elsif($line =~ /my \$scriptpath \= .+/)
	{
		print CONT2ASSTMP 'my $scriptpath = "'.$scriptpath."\";\n";
	}else{
		print CONT2ASSTMP $line;
	}
	$linecount++;
}

close(CONT2ASSTMP);
close(CONT2ASS);

open(FIXFS, "$scriptpath"."fixFrameshifts.pl");
open(FIXFSTMP, ">$scriptpath"."fixFrameshiftstmp.pl");

$linecount = 0;
while(my $line = <FIXFS>)
{
	if($linecount == 0 && $line =~ /^\#\!.+/)
	{
		print FIXFSTMP "#!".$perlpath." perl -w\n";
	}elsif($line =~ /my \$musclepath \= .+/)
	{
		print FIXFSTMP 'my $musclepath = "'.$musclepath."\";\n";
	}else{
		print FIXFSTMP $line;
	}
	$linecount++;
}

close(FIXFSTMP);
close(FIXFS);

open(ORI, "$scriptpath"."orientContig.pl");
open(ORITMP, ">$scriptpath"."orientContigtmp.pl");

$linecount = 0;
while(my $line = <ORI>)
{
	if($linecount == 0 && $line =~ /^\#\!.+/)
	{
		print ORITMP "#!".$perlpath." perl -w\n";
	}elsif($line =~ /my \$musclepath \= .+/)
	{
		print ORITMP 'my $musclepath = "'.$musclepath."\";\n";
	}else{
		print ORITMP $line;
	}
	$linecount++;
}

close(ORITMP);
close(ORI);


open(QA, "$scriptpath"."QA_stats.pl");
open(QATMP, ">$scriptpath"."QA_statstmp.pl");

$linecount = 0;
while(my $line = <QA>)
{
	if($linecount == 0 && $line =~ /^\#\!.+/)
	{
		print QATMP "#!".$perlpath." perl -w\n";
	}elsif($line =~ /my \$Rpath \= .+/)
	{
		print QATMP 'my $Rpath = "'.$Rpath."\";\n";
	}else{
		print QATMP $line;
	}
	$linecount++;
}

close(QATMP);
close(QA);


open(QLXTOSAM, "$scriptpath"."qlxToSam.pl");
open(QLXTOSAMTMP, ">$scriptpath"."qlxToSamtmp.pl");

$linecount = 0;
while(my $line = <QLXTOSAM>)
{
        if($linecount == 0 && $line =~ /^\#\!.+/)
        {
                print QLXTOSAMTMP "#!".$perlpath." perl -w\n";
        }elsif($line =~ /my \$samtoolspath \= .+/)
        {
                print QLXTOSAMTMP 'my $samtoolspath = "'.$samtoolspath."\";\n";
        }else{
                print QLXTOSAMTMP $line;
        }
        $linecount++;
}

close(QLXTOSAMTMP);
close(QLXTOSAM);

open(SAMTOQLX, "$scriptpath"."samToQlx.pl");
open(SAMTOQLXTMP, ">$scriptpath"."samToQlxtmp.pl");

$linecount = 0;
while(my $line = <SAMTOQLX>)
{
        if($linecount == 0 && $line =~ /^\#\!.+/)
        {
                print SAMTOQLXTMP "#!".$perlpath." perl -w\n";
        }elsif($line =~ /my \$samtoolspath \= .+/)
        {
                print SAMTOQLXTMP 'my $samtoolspath = "'.$samtoolspath."\";\n";
        }else{
                print SAMTOQLXTMP $line;
        }
        $linecount++;
}

close(SAMTOQLXTMP);
close(SAMTOQLX);


open(FASTA2FASTQ, "$scriptpath"."fasta2fastq.pl");
open(FASTA2FASTQTMP, ">$scriptpath"."fasta2fastqtmp.pl");

$linecount = 0;
while(my $line = <FASTA2FASTQ>)
{
        if($linecount == 0 && $line =~ /^\#\!.+/)
        {
                print FASTA2FASTQTMP "#!".$perlpath." perl -w\n";
	}else{
                print FASTA2FASTQTMP $line;
        }
        $linecount++;
}

close(FASTA2FASTQTMP);
close(FASTA2FASTQ);

open(FASTQ2FASTA, "$scriptpath"."fastq2fasta.pl");
open(FASTQ2FASTATMP, ">$scriptpath"."fastq2fastatmp.pl");

$linecount = 0;
while(my $line = <FASTQ2FASTA>)
{
        if($linecount == 0 && $line =~ /^\#\!.+/)
        {
                print FASTQ2FASTATMP "#!".$perlpath." perl -w\n";
	}else{
                print FASTQ2FASTATMP $line;
        }
        $linecount++;
}

close(FASTQ2FASTATMP);
close(FASTQ2FASTA);

open(FQPAIR2FASTA, "$scriptpath"."fqpair2fasta.pl");
open(FQPAIR2FASTATMP, ">$scriptpath"."fqpair2fastatmp.pl");

$linecount = 0;
while(my $line = <FQPAIR2FASTA>)
{
        if($linecount == 0 && $line =~ /^\#\!.+/)
        {
                print FQPAIR2FASTATMP "#!".$perlpath." perl -w\n";
	}else{
                print FQPAIR2FASTATMP $line;
        }
        $linecount++;
}

close(FQPAIR2FASTATMP);
close(FQPAIR2FASTA);

open(TRANSDNA, "$scriptpath"."translateDna.pl");
open(TRANSDNATMP, ">$scriptpath"."translateDnatmp.pl");

$linecount = 0;
while(my $line = <TRANSDNA>)
{
        if($linecount == 0 && $line =~ /^\#\!.+/)
        {
                print TRANSDNATMP "#!".$perlpath." perl -w\n";
	}else{
                print TRANSDNATMP $line;
        }
        $linecount++;
}

close(TRANSDNATMP);
close(TRANSDNA);


system("mv $scriptpath"."annotatetmp.pl $scriptpath"."annotate.pl");
system("mv $scriptpath"."runMosaiktmp.pl $scriptpath"."runMosaik2.pl");
system("mv $scriptpath"."contigMergertmp.pl $scriptpath"."contigMerger.pl");
system("mv $scriptpath"."vfattmp.pl $scriptpath"."vfat.pl");
system("mv $scriptpath"."fixFrameshiftstmp.pl $scriptpath"."fixFrameshifts.pl");
system("mv $scriptpath"."orientContigtmp.pl $scriptpath"."orientContig.pl");
system("mv $scriptpath"."QA_statstmp.pl $scriptpath"."QA_stats.pl");
system("mv $scriptpath"."qlxToSamtmp.pl $scriptpath"."qlxToSam.pl");
system("mv $scriptpath"."samToQlxtmp.pl $scriptpath"."samToQlx.pl");
system("mv $scriptpath"."fasta2fastqtmp.pl $scriptpath"."fasta2fastq.pl");
system("mv $scriptpath"."fastq2fastatmp.pl $scriptpath"."fastq2fasta.pl");
system("mv $scriptpath"."fqpair2fastatmp.pl $scriptpath"."fqpair2fasta.pl");
system("mv $scriptpath"."translateDnatmp.pl $scriptpath"."translateDna.pl");

print "Config done\n";
