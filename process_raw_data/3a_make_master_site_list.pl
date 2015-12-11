#!/bin/perl
use warnings;
#use strict;

#This script takes a list of tab separated snp tables, and merges all the sites but sorts the chromosomes in the same as the reference.
#It doesn't print groupXIX because it's the sex chromosome.

my $in = $ARGV[0]; #A list of snp tables. One file per line
#Stickleback chrom order because in the reference the chromosomes are not in an easy to sort order.
my @chromlist = qw(groupIV groupI groupVII groupII groupIX groupXIII groupXX groupVIII groupXII groupXVI groupVI groupIII groupXI groupXVIII groupXV groupX groupXIV groupXVII groupV groupXXI);

print "CHROM\tPOS";
foreach my $j (0..$#chromlist){
	my %sites;
	print STDERR "Collecting sites for $chromlist[$j] from...\n";
	open IN, $in;
	while(<IN>){
		chomp;
		print STDERR "$_\n";
		open (my $sample, '<', $_) or die "Could not open $_";;
		while (<$sample>){
			chomp;
			my @a = split(/\t/,$_);
			my $chr = $a[0];
			my $pos = $a[1];
			if ($chromlist[$j+1]){
				if ($chr eq $chromlist[$j+1]){
					goto NEXTFILE;
				}
			}
			if ($chr ne $chromlist[$j]){next;}
			$sites{$pos}++;
		}
		NEXTFILE:
		close $sample;
	}
	for my $pos (sort {$a <=> $b} keys %sites){
		print "\n$chromlist[$j]\t$pos";
	}
	close IN;
}
