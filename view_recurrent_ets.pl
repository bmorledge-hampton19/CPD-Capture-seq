use strict;
use warnings;

my $count = 0;
my $weakcount = 0;
my $pos = "";
while ( my $line = <STDIN> )
{
	chomp $line;

	if ( $line =~ /^>/ )
	{
		if ( $line =~ /^>(chr[0-9XY]+):([0-9]+)-([0-9]+)\(/ )
		{
			my $chr = $1;
			my $str = $2 + 7;
			my $end = $3 - 7;
			$pos = "${chr}:$str-$end";
		}
		else
		{
			die "Misformatted line\n";
		}
	}
	else 
	{
		if ( length( $line ) == 15 )
		{
			my $mid = substr $line, 7, 1;
			if ( $mid =~ /[AG]/ )
			{
				$line = reverse $line;
				$line =~ tr/ACGT/TGCA/;
			}

			if ( $line =~ /^[ACGT]{6,7}[CT][CT]TTCCGG/ )
			{
				print "Super match for position: $pos\n$line\n";
				$count++;
			}
			elsif ( $line =~ /^[ACGT]{6,7}[CT][CT]TTCC/ )
			{
				print "Partial super match for position: $pos\n$line\n";
                                $count++;
			}
			elsif ( $line =~ /^[ACGT]{4,5}TTCCGG/ )
			{
				print "ETS match for position: $pos\n$line\n";
                                $count++;
			}
                        elsif ( $line =~ /^[ACGT]{4,5}TTCC/ )
                        {
                                #print "Weak ETS match for position: $pos\n$line\n";
                                $weakcount++;
                        }
		}			
	}

}
print "ETS count = $count\nWeak ETS count = $weakcount\n";
