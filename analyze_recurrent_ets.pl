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
			my $str = $2;
			my $end = $3;
			$pos = "${chr}:$str-$end";
		}
		else
		{
			die "Misformatted line\n";
		}
	}
	else 
	{
		if ( length( $line ) == 15 || length ( $line ) == 16 )
		{
			my $mid = substr $line, 7, 1;
			if ( $mid =~ /[AG]/ )
			{
				$line = reverse $line;
				$line =~ tr/ACGT/TGCA/;
			}

			if ( $line =~ /^[ACGT]{6,7}[CT][CT]TTCCG/ )
			{
				print "$pos\n";
				$count++;
			}
                        elsif ( $line =~ /^[ACGT]{6,7}[CT]CTTCC/ )
                        {       
                                print "$pos\n";
                                $count++;
                        }
			elsif ( $line =~ /^[ACGT]{4,5}TTCCG/ )
			{
				print "$pos\n";
                                $count++;
			}
                        elsif ( $line =~ /^[ACGT]{3,4}CTTCC/ )
                        {
				print "$pos\n";
                                $count++;
                        }
		}			
	}

}
print STDERR "ETS count = $count\n";
