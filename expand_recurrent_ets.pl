use strict;
use warnings;

my $exp = 7;
while ( my $line = <STDIN> )
{
	chomp $line;
	my @fields = split /\t/, $line;
	my $str = $fields[1] - $exp;
	my $end = $fields[2] + $exp;
	print "$fields[0]\t$str\t$end\t$fields[3]\t$fields[4]\t$fields[5]\n"; 
}
