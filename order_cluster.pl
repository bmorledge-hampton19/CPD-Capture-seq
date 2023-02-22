#!/usr/bin/perl

use strict;
use warnings;

print STDERR "Please enter filename sorting values\n";

my $greater_file = <STDIN>;
chomp $greater_file;
open( GREATER, $greater_file ) || die "Couldn't open file: $greater_file\n";
my $headline = <GREATER>;
print STDERR "Header of sort file: $headline";
my %genelist;
while ( <GREATER> )
{
	chomp $_;
	my @fields = split /\t/, $_;
	$genelist{$fields[0]} = $fields[1];
}

close( GREATER );
print STDERR "Please enter filename for data matrix\n";
my $filename = <STDIN>;
chomp $filename;

my $fileout = $filename;
$fileout =~ s/\.txt/_filesort\.txt/ || die "filename of data matrix is wrong type!\n";
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
foreach my $acc ( sort {$genelist{$b} <=> $genelist{$a} or $a cmp $b } keys %genelist )
{
        if ( exists $matrix{$acc} )
        {
		print OUT $matrix{$acc} . "\n";
	}
}

