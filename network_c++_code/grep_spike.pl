#!/usr/bin/perl
use strict;
use warnings;

my $dir = $ENV{'PWD'};
#print $dir;
opendir(DIR, $dir) or die $!;
#my @allfiles = readdir(DIR);
my %rs_freq = ();
my %fs_freq = ();
my %ib_freq = ();
open OUTRS, ">fi_rs.dat" or die $!;
open OUTFS, ">fi_fs.dat" or die $!;
open OUTIB, ">fi_ib.dat" or die $!;


foreach my $file(readdir(DIR)){
	next if ($file!~/spikeout/i);
	my @sfile = split /\_/, $file;
	print $file."\n";
	my $num_rs = `grep '^0\t' $file | wc -l`; # grep the total number of spikes in 10 second interval for each type of neuron
	my $num_ib = `grep '^5\t' $file | wc -l`;
	my $num_fs = `grep '^13\t' $file | wc -l`;

	my $current = `echo $sfile[1] | sed s/d/./g`;
	$current=~s/\n//g;
	$num_rs=~s/\n//g;
	$num_ib=~s/\n//g;
	$num_fs=~s/\n//g;
#print "current is $current\n";
	$rs_freq{$current} = $num_rs/10;
	$ib_freq{$current} = $num_ib/10;
	$fs_freq{$current} = $num_fs/10;
		
}
	
	foreach my $key(sort {$a <=> $b} keys %rs_freq){
#print " key is $key\n";
		print OUTRS $key."\t".$rs_freq{$key}."\n";
print OUTIB $key."\t".$ib_freq{$key}."\n";
		print OUTFS $key."\t".$fs_freq{$key}."\n";
	}

close(OUTRS);
close(OUTIB);
close(OUTFS);

print "Programme ends\n";
