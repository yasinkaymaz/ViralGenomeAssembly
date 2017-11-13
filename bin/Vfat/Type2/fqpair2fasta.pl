#######!/share/pkg/perl/5.18.1/bin/ perl -w

# Copyright © 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

use strict;

my $fastq1 = shift || die("Usage : perl fqpair2fasta.pl input1.fq input2.fq output\nWill create the files output.fa and output.qual\n");
my $fastq2 = shift || die("Usage : perl fqpair2fasta.pl input1.fq input2.fq output\nWill create the files output.fa and output.qual\n");
my $output = shift || die("Usage : perl fqpair2fasta.pl input1.fq input2.fq output\nWill create the files output.fa and output.qual\n");

my $scriptpath = "/project/umw_jeffrey_bailey/share/bin_sync/Vfat/";

system("perl $scriptpath/fastq2fasta.pl $fastq1 $output"."_1");
system("perl $scriptpath/fastq2fasta.pl $fastq2 $output"."_2");
system("cat $output"."_1.fa $output"."_2.fa > $output.fa");
system("cat $output"."_1.qual $output"."_2.qual > $output.qual");
system('perl -p -i -e \'s/>(.+?) .+/>$1/g\' '.$output.'.fa');
system('perl -p -i -e \'s/>(.+?) .+/>$1/g\' '.$output.'.qual');

system("rm $output"."_1.fa");
system("rm $output"."_1.qual");
system("rm $output"."_2.fa");
system("rm $output"."_2.qual");
