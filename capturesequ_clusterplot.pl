#!/usr/bin/perl
# ctcf_mutation_unique_indivtfbs.pl
use strict;
use warnings;

use lib '../';

use ClusterPosCoord;
use CPDHumanReads;

print STDERR "Enter filename for output .txt file\n";
my $outfile = <STDIN>;

print STDERR "Enter length of window (in bp)\n";
my $window = <STDIN>;

chomp $window;
chomp $outfile;

$outfile =~ s/\.txt/_${window}bp\.txt/ || die "wrong file format\n";

my $bed_file = $outfile;
$bed_file =~ s/\.txt/\.bed/ || die "output file must be a .txt file\n";

open ( OUT, ">$outfile" ) || die "couldn't open output file\n";
open ( BED, ">$bed_file" ) || die "couldn't open bed file\n";

# ask for probe filename to analyze
print STDERR "Enter filename of plus strand reads\n";
my $plusfile = <STDIN>;
chomp($plusfile);

print STDERR "Enter filename of minus strand reads\n";
my $minusfile = <STDIN>;
chomp($minusfile);

# add filename here for probe values 
#my $filename = "h2bwtdata.txt";

print STDERR "Loading TF coordinates\n";
my $tfsites = ClusterPosCoord->new();
print STDERR "Loading Probe Values\n";
my $reads = CPDHumanReads->new($plusfile, $minusfile);

=pod
print STDERR "Loading ETS motif offsets\n";
my $offsetfile = "../nfy_motifs_offsets.txt";
open( OFFSET, $offsetfile ) || die "couldn't open ETS offset file\n";
my %ets_offset = ();
while ( my $line = <OFFSET> )
{
	chomp $line;
	my @fields = split /\t/, $line;
	# file gives offset of first G; set to offset of second G by adding 1
	my $offset = $fields[1] + 1;
	$ets_offset{$fields[0]} = $offset;
}
=cut

my %tf = $tfsites->get_tf_boundaries();
my $tf_filename = $tfsites->get_tf_filename();

# flip strand for calculating other strand of motif 
my %flipstrand = ( "+" => "-", "-" => "+" );

#print header
print OUT "Gene position binding site data from: $tf_filename";
print OUT "\tSequencing data from file: $plusfile\t$minusfile\n";
print OUT "Genome Region midpoint (1-based)\tRegion Midpoint (strand}";
for (my $i = (2 * $window ); $i >= 0; $i-- )
{
        my $pos = $window - $i;
        print OUT "\t$pos";
}
print OUT "\n";

foreach my $chr (sort keys %tf)
{
	print STDERR "Starting $chr\n";
	
	# locate TFBS coordinates throughout chromosome to avoid impingement during overlap
	my @tfpos = @{$tf{$chr}};

	my %plusreads = $reads->get_plus_reads_for_chromosome($chr);
	my $num_plusreads = scalar keys %plusreads;
	my %minusreads = $reads->get_minus_reads_for_chromosome($chr);
	my $num_minusreads = scalar keys %minusreads;
	print STDERR "$chr reads: $num_plusreads plus reads and $num_minusreads minus reads\n";

	# previous tfmidpoints;
	my %prevtfmids = ();
        for ( my $i = 0; $i < scalar @tfpos; $i++ )
        {
                my $tfstart = $tfpos[$i]->[0];
                my $tfend = $tfpos[$i]->[1];
		my $tfmotif = $tfpos[$i]->[2];
		my $tfstrand = $tfpos[$i]->[3];
		my $bedstart = $tfstart;
		# increase start position by 1, since wig files (data files) are one-based, and bed files (position files) are zero-based (only start, not end coordinate)
		$tfstart++;

		my @plusmutationval = ();
		my @minusmutationval = ();
		my @bothval = ();
		my $tfmidpoint;
		if ( $tfstrand eq "+" )
		{
	        	$tfmidpoint = ( $tfstart + $tfend ) / 2.0;

	        	if ( $tfmidpoint - int($tfmidpoint) != 0.0 )
        		{
	                	$tfmidpoint = int ( $tfmidpoint + 1);
			}

			# don't repeat analysis of TF with same midpoint as previous TF
			if ( exists $prevtfmids{$tfmidpoint} && $prevtfmids{$tfmidpoint} == 1 )
			{
				next;
			}
	
			$prevtfmids{$tfmidpoint} = 1;
	
			# expand size of window
			my $winstart = $tfmidpoint - $window;
			my $winend = $tfmidpoint + $window;

			for ( my $j = $winstart; $j <= $winend; $j++ )
			{	
				my $relpos = $j - $winstart;
	
				if ( exists $plusreads{$j} )
				{
					$plusmutationval[$relpos] += $plusreads{$j};
					$bothval[$relpos] += $plusreads{$j};	
				}
			
	        	        if ( exists $minusreads{$j} )
	                	{
	                                $minusmutationval[$relpos] += $minusreads{$j};
					$bothval[$relpos] += $minusreads{$j};
	                        }
			}
		}
                elsif ( $tfstrand eq "-" )
                {
                        $tfmidpoint = ( $tfstart + $tfend ) / 2.0;

                        if ( $tfmidpoint - int($tfmidpoint) != 0.0 )
                        {
				# rounding down on minus strand is same as rounding up on plus strand 
                                $tfmidpoint = int ( $tfmidpoint);
                        }

                        # don't repeat analysis of TF with same midpoint as previous TF
                        if ( exists $prevtfmids{$tfmidpoint} && $prevtfmids{$tfmidpoint} == 1 )
                        {
                                next;
                        }

                        $prevtfmids{$tfmidpoint} = 1;

                        # expand size of window
                        my $winstart = $tfmidpoint - $window;
                        my $winend = $tfmidpoint + $window;

                        for ( my $j = $winend; $j >= $winstart; $j-- )
                        {
                                my $relpos = $winend - $j;

				# switch strands because motif is on minus strand
                                if ( exists $plusreads{$j} )
                                {
                                        $minusmutationval[$relpos] += $plusreads{$j};
					$bothval[$relpos] += $plusreads{$j};
                                }

                                if ( exists $minusreads{$j} )
                                {
                                        $plusmutationval[$relpos] += $minusreads{$j};
					$bothval[$relpos] += $minusreads{$j};
                                }
                        }
		}
		else
		{
			die "No strand information for $tfmotif\n";
		}

		#print OUT "$chr:$bedstart-$tfend\t${tfmotif}_${tfmidpoint}_${tfstrand}\tBothStrands";
		print BED "$chr\t$bedstart\t$tfend\t$chr:$bedstart-$tfend\t$tfmotif\t$tfstrand\n";
		print OUT "$chr:$tfmidpoint\t$tfmotif";

		# calculate median, based off of code from: https://stackoverflow.com/questions/5119034/using-perl-to-find-median-mode-standard-deviation
		my @data = ();
		my $data_string = "";
                for (my $i = 0; $i <= ( 2 * $window ); $i++ )
                {
                        my $val = 0;
                        if ( exists $bothval[$i] )
                        {
                                $val += $bothval[$i];
				push @data, $val;
                        }
				$data_string .= "\t$val";
			}
=pod			
               	my $median;
                my $midindex = int ( (scalar @data) / 2 );
                my @sorted = sort {$a <=> $b} @data;
                if ( (scalar @data) % 2 == 1 )
                {       
                        $median = 1.0 * $sorted[$midindex];
                }
                else    
                {       
                        $median = ( $sorted[$midindex] + $sorted[$midindex - 1])/ 2.0;
                }

		my $snr = "--";
		if ( $median > 0  )
		{
			if ( exists $bothval[$window] )
			{
				$snr = 1.0 * $bothval[$window] / $median;
			}
			else
			{
				$snr = 0;
			}
		} 
=cut
		
		print OUT "$data_string\n";
	}
}
