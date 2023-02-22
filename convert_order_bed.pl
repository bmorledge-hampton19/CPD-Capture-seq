use strict;
use warnings;

my $win = 10;
while ( my $line = <STDIN> )
{
	chomp $line;
	my @fields = split /\t/, $line;
	my $chr;
	my $mid;
	if ($fields[0] =~ /^(chr[0-9XY]+):([0-9]+)$)/ )
	{
		$chr = $1;
		$mid = $2;
	}
	else
	{
		die "Misformatted line: $line\n";
	}
	my $str = $mid - $win;
	$str--; #0-based
	my $end = $mid + $win;
	print "$chr\t$str\t$end\t$fields[0]\t.\t.\n";
	#wont work don't know strand
}
