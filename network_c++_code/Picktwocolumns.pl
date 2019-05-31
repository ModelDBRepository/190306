#!/usr/bin/perl

#
# Picktwocolumns.pl
# Pick two columns from a multi column file
# Author: Ernest Ho
# Date: April 30th, 2010
#

use strict;
my $infile=$ARGV[0];
my $fcol=$ARGV[1];
my $scol=$ARGV[2];

open(INPUT, "< $infile") or die "Couldn\'t open file for reading\n";
while (<INPUT>){
    next if ($_=~m/^\#/);
    my @thisline = split /\s/, $_;
    #print "thisline[1] is ".$thisline[1]."\n";
    print $thisline[$fcol]."\t".$thisline[$scol]."\n";

}

