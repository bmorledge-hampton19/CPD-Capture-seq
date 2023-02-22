#!/usr/bin/perl

use strict;
use warnings;

my $count_filtered = 0;
while( <STDIN> )
{
	if ( $_ =~ /^pUC19/ )
	{	$count_filtered++;	}
	else
	{	print $_;	}	

}
print STDERR "Filtered lines: $count_filtered\n";
