#!/usr/bin/perl

#
# rearrange_rows.pl
# rearrange rows according to the first entry
# Author: Ernest Ho
# Date: May 10th, 2013.
#

use strict;
my $infile=$ARGV[0];
my @g_array;
#my $fcol=$ARGV[1];
#my $scol=$ARGV[2];

open(INPUT, "< $infile") or die "Couldn\'t open file for reading\n";
while (<INPUT>){
    next if ($_=~m/^\#/);
    my @thisline = split /\s/, $_;
		push @g_array, \@thisline;
#$_=~s/\n/\\n/g;
#   print $_;

}

@g_array = sort {@{$a}[0] <=> @{$b}[0]} @g_array;
#@g_array = sort {@{$a}[1] <=> @{$b}[1]} @g_array;

foreach my $entries(@g_array){
print @{$entries}[0]."\t".@{$entries}[1]."\t".@{$entries}[2]."\n";
	
}
print "Programme ends"."\n";
