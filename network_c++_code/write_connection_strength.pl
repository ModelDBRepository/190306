#!/usr/bin/perl
use strict;

my $directory = $ARGV[0];
my @directory_path = split "\/", $directory;
my $last_path = pop @directory_path;
print "last path is ".$last_path."\n";
my @path = split "_", $last_path;
my $numproc = pop @path;
print "Last argument is ".$ARGV[5]."\n";


for (my $i=0; $i<$numproc; $i++){
	print "Num proc is ".$numproc."\n";
	my $file = $directory."/"."syn_sorted_proc_template_".$i.".DAT";
	my $outfile =$directory."/"."syn_sorted_".$ARGV[5]."_proc_".$i.".DAT";
	next if (-s $outfile); # File exists has non-zero size
	open (OUT, "> $outfile") or die "Cannot open file\n";
	my $out0="sed -e s/EtoE/$ARGV[1]/g -e s/EtoI/$ARGV[2]/g -e s/ItoE/$ARGV[3]/g -e s/ItoI/$ARGV[4]/g $file  ";
		my $outfinal=$out0;
#$out4=`$out3`;
#print $outfinal;
print OUT `$outfinal`;
#print $out0;
	close(OUT);
}

print "Programme ends\n";


#open(INPUT, "< test_par_init.dat");
#print OUT $cort->calPoissonIndex(0)."\t"."-15\t0.0\n\n";
#$cort->printParInit(\*OUT);
#$cort->printPar("test_par.dat");
#$cort->printInit("test_init.dat");
#$cort->printConnectionMatrix(\*OUT);
#close(OUT);



#$cort->assignConnect;
#$cort->printConnectionMatrix();


=pod
print "L23PYR to L23PYR total connection is ".$cort->{"connectionCount"}->{"L23PYR"}->{"L23PYR"}."\n";
print "L23PYR to L5PYR total connection is ".$cort->{"connectionCount"}->{"L23PYR"}->{"L5PYR"}."\n";
print "L23PYR to L4ST total connection is ".$cort->{"connectionCount"}->{"L23PYR"}->{"L4ST"}."\n";
print "L5PYR to L23PYR total connection is ".$cort->{"connectionCount"}->{"L5PYR"}->{"L23PYR"}."\n";
print "L5PYR to L5PYR total connection is ".$cort->{"connectionCount"}->{"L5PYR"}->{"L5PYR"}."\n";
print "L5PYR to L6PYR total connection is ".$cort->{"connectionCount"}->{"L5PYR"}->{"L6PYR"}."\n";

=cut


#print $cort->{"cellDensity"}->{"DOUB"};
#print split(" ", keys %{$cort->{"distanceHash"}});
#my $count = 0;
#print "$count"."\n";
