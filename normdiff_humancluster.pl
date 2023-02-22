#!/usr/bin/perl

use strict;
use warnings;

print "Enter name of 0hr cluster file\n";
my $mutclusterfile = <STDIN>;
chomp $mutclusterfile;

my $diffclusterout = $mutclusterfile;

print "Enter name of Naked DNA cluster file\n";
my $wtclusterfile = <STDIN>;
chomp $wtclusterfile;

print "Enter normalization factor for naked DNA (multipy)\n";
my $norm = <STDIN>;
chomp $norm;

$diffclusterout =~ s/\.txt/_normdiff\.txt/ || die "file not .txt file\n";
open ( MUT, "$mutclusterfile" ) || die "Couldn't open file $mutclusterfile\n";
open ( WT, "$wtclusterfile" ) || die "Couldn't open file $wtclusterfile\n";
open ( OUT, ">$diffclusterout" ) || die "Couldn't open file $diffclusterout\n";

if ( $mutclusterfile =~ /_centered/ || $mutclusterfile =~ /_center/ )
{
        ;
}
else
{
        # remove top header because (normally removed by logtransform_center.pl program)
        my $topline = <MUT>;
}

if ( $wtclusterfile =~ /_centered/ || $wtclusterfile =~ /_center/ ) 
{
        ;
}
else
{
        # remove top header because (normally removed by logtransform_center.pl program)
        my $top = <WT>;
}

my $header = <MUT>;
my $otherhead = <WT>;
if ( $header eq $otherhead )
{
	print OUT $header;
}
else
{
	die "mismatched headers\n";
}

my @diffdata = ();
while ( my $line = <MUT> )
{
	chomp $line;

	# to keep empty last tab, see https://stackoverflow.com/questions/3711649/perl-split-with-empty-text-before-after-delimiters
	my @mutvals = split /\t/, $line, -1;

	my $wt = <WT> || die "error in number of WT lines\n";
	chomp $wt;
	my @wtvals = split /\t/, $wt, -1;

	if ( $wtvals[0] eq $mutvals[0] )
	{
		print OUT $mutvals[0];
	}
	else
	{
		die "Mismatched accs for lines: $line\t$wt\n";
	}

        if ( $wtvals[1] eq $mutvals[1] )
        {
                print OUT "\t$mutvals[1]";
        }
        else
        {
                die "Mismatched accs for lines: $line\t$wt\n";
        }
	if ( scalar @wtvals != scalar @mutvals )
	{
		die "mismatched entries for lines: $line\t$wt\n";
	}
	for ( my $i = 2; $i < scalar @mutvals; $i++ )
	{
		my $diffval = "";
		if ( $mutvals[$i] ne "" && $wtvals[$i] ne "" )
		{
			$diffval = $mutvals[$i] - ( $norm * $wtvals[$i] );
			push @diffdata, $diffval;

		}
		print OUT "\t$diffval";
	}
	print OUT "\n";
}

# calculate median, based off of code from: https://stackoverflow.com/questions/5119034/using-perl-to-find-median-mode-standard-deviation

my $median;
my $midindex = int ( (scalar @diffdata) / 2 );
my @sorted = sort {$a <=> $b} @diffdata;
if ( (scalar @diffdata) % 2 == 1 )
{
	$median = 1.0 * $sorted[$midindex];
}
else
{
	$median = ( $sorted[$midindex] + $sorted[$midindex - 1])/ 2.0;
}


print STDERR "Median of differences is $median\n";

