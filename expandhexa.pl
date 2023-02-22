use strict;
use warnings;

while ( my $line = <STDIN> )
{
	chomp $line;
	my @fields = split /\t/, $line;
	$fields[2] += 2;
	$fields[1] -= 2;

	print "$fields[0]";
	for ( my $i = 1; $i < scalar @fields; $i++ )
	{
		print "\t$fields[$i]";
	}	
	print "\n";
}
