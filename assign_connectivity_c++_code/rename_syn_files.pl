#!/usr/bin/perl
use strict;

my $dir=$ARGV[0];
opendir(DIR, $dir) || die("Cannot open directory");

my @files=readdir(DIR);
print @files;

foreach my $file(@files){
	if ($file=~/^syn_sorted_proc_/){
		next if ($file=~/template/);
		my @split=split("\_", $file);
		my $rearranged="syn_sorted_proc_template_".$split[3];
		my $copycommand="cp $dir/$file $dir/$rearranged ";
		`$copycommand`;
	}
}
