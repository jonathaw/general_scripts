#!/usr/bin/perl

# Written by David La
# Updated Tue May 29 03:07:47 PDT 2012

# Description: Fetches the first DNA sequnece out of the DNAworks output!

my $usage = "getSeq.pl <dnaworks.output.txt>\n";

use strict;

my $file = $ARGV[0] or die $usage;
my $five_prime = $ARGV[1];
my $three_prime = $ARGV[2];

my $seq;

open(FILE,"$file") or die "Cannot open file $file\n";

my $start1;
my $start2;
while (<FILE>) {
	chomp;
	
	if (/ The DNA sequence #\s+1/) {
		$start1++;
	}
	elsif ($start1) {
		$start1 = 0;
		$start2++ if $_ =~ /^ -------/;
	}
	elsif ($start2 and $_ =~ /^ -------/) {
		$start2 = 0;
	}
	elsif ($start2) {
		$_ =~ /\d+\s+(\w+)/;
		$seq .= $1;
	}
	
}


print "$five_prime$seq$three_prime\n";

close(FILE);
