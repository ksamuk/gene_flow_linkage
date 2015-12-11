#!/bin/perl
use warnings;
use strict;

#This script takes a list of sites, that may have snps with it or may not. It loads those lines and also a new set of samples. It puts the samples in or it fills it with NN. The list and the new set must have sites sorted in the same manner and the list should have all possible sites.

my $new = $ARGV[0];

my $counter = 0;
my $first_row = 0;
my $current_chr;
my %data_hash;
my %non_empty_pos;
my $header;
my $last_pos;
my $header_printed;
my $new_col_number;
while(<STDIN>){
	chomp;
	if ($. == 1){
		$header = "$_";
		next;
	}
	my $line = $_;
	my @a = split(/\t/,$_);
	my $chr = $a[0];
	my $pos = $a[1];
	unless ($current_chr){
		$current_chr = $chr;
	}
	if ((($counter > 0) and($counter % 10000000 == 0)) or (eof) or ($current_chr ne $chr)){
		print STDERR "At line $chr\t$pos. The previous chr was $current_chr\n";
		&pull_data;
		unless($header_printed){
			print "$header";
			$header_printed++;
		}
		foreach my $key (sort {$a <=> $b}  keys %data_hash){
			if ($non_empty_pos{$key}){
				print "\n$data_hash{$key}";
			}else{
				print "\n$data_hash{$key}";
				foreach my $i (2..$new_col_number){
					print "\tNN";
				}
			}
			unless ($data_hash{$key}){
				print STDERR "\n$key";
			}
		}
		undef(%data_hash);
		undef(%non_empty_pos);
		$counter = 0;
	}
	$data_hash{$pos} = $line;
	$current_chr = $chr;
	$last_pos = $pos;
	$counter++;
}

#This subroutine opens the new data. It goes through each row and asks if it is the first row in the old data. If it isn't then it shifts off rows from the old data and prints them, adding NN for the new data until it matches.
sub pull_data{
	open NEW, $new;
	while(<NEW>){
		chomp;
		my @a = split(/\t/,$_);
		my $chr = $a[0];
		my $pos = $a[1];
		unless ($new_col_number){
			$new_col_number = $#a;
		}
		if ($first_row == 0){
			foreach my $i (2..$#a){
				$header .= "\t$a[$i]";
			}
			$first_row++;
			next;
		}elsif($chr eq "CHROM"){
			next;
		}
		if ($chr ne $current_chr){
			next;
		}
		if ($pos > $last_pos){
			return;
		}
		if ($data_hash{$pos}){
			$non_empty_pos{$pos}++;
			foreach my $i (2..$#a){
				$data_hash{$pos} .= "\t$a[$i]";
			}
		}
	}
	close NEW;
}
