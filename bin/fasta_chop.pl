#!/usr/bin/perl

use strict;
use warnings;

my $len = 5000;
my $over = 0;
my ($seq_id, $seq);

while (<>) {
    chomp;
    if (m/^>/) { $seq_id = $_; } else { $seq .= $_; }
}

for (my $i = 1; $i <= length $seq; $i += ($len - $over)) {
    my $s = substr ($seq, $i - 1, $len);
    print "$seq_id:$i-", $i + (length $s) - 1, "\n$s\n";
}


