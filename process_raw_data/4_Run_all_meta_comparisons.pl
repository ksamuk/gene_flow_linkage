#!bin/perl
use warnings;
#use strict;
#This script takes the popfile, groupfile and snpfile and produces pbs files for each population combination. The pbs will run the SNPtable2fstdxy and sliding_window scripts. 
my $popfile = "samplelist.txt";
my $groupfile = "grouplist.txt";
my $SNPfile = "snptable.tab";

my %pop;
my @poplist;
open POP, $popfile;
while(<POP>){
	chomp;
	my @a = split(/\t/,$_);
	$pop{$a[0]} = $a[1];
	push(@poplist, $a[1]);
}
close POP;

my @uniq_pops = do { my %seen; grep { !$seen{$_}++ } @poplist };

my %location;
my %type;
open GROUP, $groupfile;
while(<GROUP>){
	chomp;
	my @a = split(/\t/,$_);
	$location{$a[0]} = $a[1];
	$type{$a[0]} = $a[2];
}
my %counter;
my $counts;
foreach my $pop1 (0..($#uniq_pops-1)){
	foreach my $pop2 (($pop1+1)..$#uniq_pops){
		$counts++;
		my $type_comp;
		unless ($type{$uniq_pops[$pop2]}){
			print STDERR "$uniq_pops[$pop2]\n";
		}
		if ($type{$uniq_pops[$pop1]} eq $type{$uniq_pops[$pop2]}){
			$type_comp = "S";
		}else{
			$type_comp = "D"
		}
		my $location_comp;
		if ($location{$uniq_pops[$pop1]} eq $location{$uniq_pops[$pop2]}){
			$location_comp = "para";
		}else{
			$location_comp = "allo"
		}
		if (-e "$uniq_pops[$pop1].$uniq_pops[$pop2].$location_comp.$type_comp.slidingwindow.150000.txt"){
			open SW, '<', "$uniq_pops[$pop1].$uniq_pops[$pop2].$location_comp.$type_comp.slidingwindow.150000.txt";
			chomp (my @rows = <SW>);
			close SW;
			if ($rows[-1] =~ /^groupXXI/){
				my @tmpdata = split(/\t/, $rows[-1]);
				if ($tmpdata[1] > 11400000){
					print STDERR "Skipping $uniq_pops[$pop1].$uniq_pops[$pop2].$location_comp.$type_comp because it was already run.\n";
					next;
				}
			}
		}
		if (-e "$uniq_pops[$pop1].$uniq_pops[$pop2].$location_comp.$type_comp.stats.txt.gz"){
			next;
		}
#		print "Comparison of $uniq_pops[$pop1].$uniq_pops[$pop2].$location_comp.$type_comp\n";
		$counter{$location_comp.$type_comp}++;
		open (my $groupfile,">", "$uniq_pops[$pop1].$uniq_pops[$pop2].$location_comp.$type_comp.grouplist.txt");
		print $groupfile "$uniq_pops[$pop1]\t1\n$uniq_pops[$pop2]\t2";
		close $groupfile;
		open(my $filename,">","$uniq_pops[$pop1].$uniq_pops[$pop2].$location_comp.$type_comp.pbs");
		print $filename qq[

#!/bin/bash
#PBS -S /bin/bash
#PBS -l walltime=360:00
#PBS -l mem=2GB
#PBS -l nodes=1:ppn=1
#PBS -N "$uniq_pops[$pop1].$uniq_pops[$pop2].$location_comp.$type_comp".comp
#PBS -M EMAIL\@address.com
#PBS -m bea

comparisons='comp'

cd \$PBS_O_WORKDIR
JOBINFO="$uniq_pops[$pop1].$uniq_pops[$pop2].$location_comp.$type_comp"_\${PBS_JOBID}
echo \"Starting run at: \`date\`\" >> \$JOBINFO

#Run analysis
cat $SNPfile | perl 4a_SNPtable2Fstdxy_filtered_speedy.pl $popfile $uniq_pops[$pop1].$uniq_pops[$pop2].$location_comp.$type_comp.grouplist.txt > $uniq_pops[$pop1].$uniq_pops[$pop2].$location_comp.$type_comp.stats.txt
cat $uniq_pops[$pop1].$uniq_pops[$pop2].$location_comp.$type_comp.stats.txt | perl 4b_SlidingWindow_v5.pl 150000 > $uniq_pops[$pop1].$uniq_pops[$pop2].$location_comp.$type_comp.slidingwindow.150000.txt
cat $uniq_pops[$pop1].$uniq_pops[$pop2].$location_comp.$type_comp.stats.txt | perl 4b_SlidingWindow_v5.pl 75000 > $uniq_pops[$pop1].$uniq_pops[$pop2].$location_comp.$type_comp.slidingwindow.75000.txt
gzip $uniq_pops[$pop1].$uniq_pops[$pop2].$location_comp.$type_comp.stats.txt

echo \"Program finished with exit code \$? at: \`date\`\" >> \$JOBINFO"

];
		my $qsub = "qsub -c s $uniq_pops[$pop1].$uniq_pops[$pop2].$location_comp.$type_comp.pbs";
#		system($qsub);
		

		#THIS MAKES THE SCRIPT END AFTER SUBMITTING 3 JOBS. COMMENT OUT THE FOLLOWING THREE LINES TO MAKE IT RUN ALL
		if ($counts > 3){
#			exit;
		}
	}
}


#foreach my $key (keys %counter){
#	print "$key = $counter{$key}\n";
#}
