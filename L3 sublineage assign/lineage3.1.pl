#!/usr/bin/perl
use strict;
use warnings;

my %sites;

open my $F1, "<", $ARGV[0] or die $!;
while (<$F1>) {
    chomp;
    next if /^$/;
    my @a = split /\t/;
    $sites{$a[0]} = 1;
}
close $F1;

# L3.1.1 sublineage diagnostic SNPs
my %snp = (
    2118846 => "L3.1.1.1",
    1914217 => "L3.1.1.2",
    2812520 => "L3.1.1.3",
    3350236 => "L3.1.1.4",
    1134143 => "L3.1.1.5",
    3592529 => "L3.1.1.6",
);

my $assigned = 0;

foreach my $pos (keys %sites) {
    if (exists $snp{$pos}) {
        print "$snp{$pos}_$ARGV[0]\n";
        $assigned = 1;
        last;
    }
}

if ($assigned == 0) {
    print "???_$ARGV[0]\n";
}