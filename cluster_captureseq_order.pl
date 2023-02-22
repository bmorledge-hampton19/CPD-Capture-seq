#!/usr/bin/perl

use strict;
use warnings;

# get file with order of capture regions
my $order_file = "../Compiled_capturearray_simple_order.txt";
open( ORDER, $order_file ) || die "Couldn't open file: $order_file\n";

my %genelist;
while ( <ORDER> )
{
	chomp $_;
	my @fields = split /\t/, $_;
	$genelist{$fields[0]} = $fields[1];
}

close( ORDER );

print STDERR "Please enter filename for data matrix\n";
my $filename = <STDIN>;
chomp $filename;

my $fileout = $filename;
$fileout =~ s/\.txt/_order\.txt/ || die "filename of data matrix is wrong type!\n";
open ( OUT, ">$fileout") || die "couldn't open file\n";

open( FILE, $filename ) || die "Couldn't open file: $filename\n";

if ( $filename =~ /_centered/ || $filename =~ /_center/ || $filename =~ /_norm/ )
{
	;
}
else
{
	# remove top header because (normally removed by logtransform_center.pl program)
	my $topline = <FILE>;
}

my $header = <FILE>;
my %matrix;
while ( my $line = <FILE> )
{
	chomp $line;
	my @fields = split /\t/, $line;
	$matrix{$fields[0]} = $line;
}

#print less list
print OUT $header;
foreach my $acc ( sort {$genelist{$a} <=> $genelist{$b} or $a cmp $b } keys %genelist )
{
        if ( exists $matrix{$acc} )
        {
		print OUT $matrix{$acc} . "\n";
	}
}

