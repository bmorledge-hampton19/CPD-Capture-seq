#!/usr/bin/perl 
# intersect_TFBS.pl

use strict;
use warnings;

print STDERR "Enter filename of wig to intersect with TF coordinates\n";
my $tfbs_file = <STDIN>;
chomp $tfbs_file;
open ( TFBS, $tfbs_file ) || die "Couldn't open file: $tfbs_file\n";

# open DHS coordinates
print STDERR "Enter filename of TF coordinate file:\n";
my $dhs_file = <STDIN>;
chomp $dhs_file;
open (DHS, $dhs_file) || die "Couldn't open promoter file: $dhs_file\n";
print STDERR "processing TF coordinates...\n";

# process promoters into lookup hash (1 => promoter, undef/not exist ==> not promoter)
# I am assuming these coordinates are 1 based
my %DHS_coord;
my $chr = "";
while ( my $line = <DHS> )
{
	chomp $line;
	my @fields = split /\t/, $line;
	if ( $fields[0] =~ /(chr[0-9XYM_]+)/ )
	{
		my $temp = $1;
		if ( $chr ne $temp )
		{
			print STDERR "Starting to process $temp\n";
			$chr = $temp;
		}

		push @{$DHS_coord{$chr}}, [$fields[1], $fields[2], $fields[3]];
	}
	else
	{
		die "Misformatted line: $line\n";
	}
}
close (DHS);

# intersect TF coordinates with DHS --> any overlap means active DHS TF

$tfbs_file =~ s/\.wig//;

my $InDHS_file = $tfbs_file . "_ets34dipy.txt";

#my $notInDHS_file = $tfbs_file . "_notDHS.bed";

open (DHS, ">$InDHS_file") || die "Couldn't open InDHS file: $InDHS_file\n";
#open (NOTDHS, ">$notInDHS_file") || die "Couldn't open NOTInDHS file: $notInDHS_file\n";

print STDERR "Starting processing TFBS coordinates...\n";
$chr = "";
my %DHS_lookup;

my %wigcount = ();
while ( my $line = <TFBS>)
{
        chomp $line;
        $line =~ s/^M//g;
        if ( $line eq "" )
        {
                next;
        }
        elsif( $line =~ /chrom=(chr[0-9A-Za-z_]+) / )
        {
                $chr = $1;
                print STDERR "Starting to process $chr\n";
                # empty hash
                %DHS_lookup = ();
                foreach my $dhs ( @{$DHS_coord{$chr}} )
                {
                        my $start = $dhs->[0];
                        # start coordinate is 0-based, so convert to 1-based for consistency
                        $start++;
                        my $end = $dhs->[1];
                        for ( my $i = $start; $i <= $end; $i++ )
                        {
                                # Define this nucleotide as a promoter
                                $DHS_lookup{$i} = $dhs->[2];
                        }
                }

        }
        else
        {
                my @fields = split /\t/, $line;

                # hash of hashes
		my $pos = $fields[0];
		my $count = $fields[1];
		if ( exists $DHS_lookup{$pos} )
		{
			my $site = $DHS_lookup{$pos};
			if ( exists $wigcount{$site} )
			{
				$wigcount{$site} .= "\t";
			}
			$wigcount{$site} .= "$pos\t$count";
				
		}

	}
}


foreach my $tf ( sort keys %wigcount )
{
	print DHS "$tf\t$wigcount{$tf}\n";
}
