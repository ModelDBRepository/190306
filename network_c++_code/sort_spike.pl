#!/usr/bin/perl

#
# CountAP.pl
# Script to count the number of APs in dat file
# Author: Ernest Ho
# Date: September 10th 2009
#

use strict;
#print "Please enter input file:\n";
#my $infile = <STDIN>;
#chomp($infile);

#print "Please enter output file:\n";
#my $outfile = <STDIN>;
#chomp($outfile);

my $infile=$ARGV[0];
# my $numNeurons=$ARGV[1] || 20;
my @splitinfile=split /\./, $infile;
my $infilesuffix = pop @splitinfile;
my $outfile=join("_",@splitinfile, "sorted");
$outfile=join(".", $outfile, $infilesuffix);

#open(OUTPUT, "> ./$outfile") or die "Couldn\'t open $infile for writing\n";
#open(INPUT, "< ./$infile") or die "Couldn\'t open file for reading\n";
open(INPUT, "< $infile") or die "Couldn\'t open file for reading\n";
open(OUTPUT, "> $outfile") or die "Couldn\'t open file for writing\n";

my @thisline;
my @lines;
my @times;
while (<INPUT>){
	push @lines, $_;
}

foreach my $element(sort my_sort @lines){
#foreach my $element(@lines){
#print OUTPUT $element;
my @thisline = split /\s/, $element;
print OUTPUT $thisline[0]."\t".$thisline[1]."\t".$thisline[2]."\t".$thisline[3]."\t".$thisline[4]."\t".$thisline[5]."\t".$thisline[6]."\n";
}


close(OUTPUT);
close(INPUT);

sub my_sort{
	my @line_1 = split /\s/, $a;
	my @line_2 = split /\s/, $b;
	return $line_1[0]<=>$line_2[0];
}
