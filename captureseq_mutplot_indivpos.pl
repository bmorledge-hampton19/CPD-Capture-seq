#!/usr/bin/perl
# ctcf_mutation_unique_indivtfbs.pl
use strict;
use warnings;

use lib '../';

use GenomePosCoord;
use CPDHumanReads;

print STDERR "Enter filename for output .txt file\n";
my $outfile = <STDIN>;

chomp $outfile;

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
my $tfsites = GenomePosCoord->new();
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
print OUT "Mutation location (bed)\tMutation type/count\tCPD count\n";

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
		my $tfname = $tfpos[$i]->[4];
		my $bedstart = $tfstart;
		# increase start position by 1, since wig files (data files) are one-based, and bed files (position files) are zero-based (only start, not end coordinate)
		$tfstart++;
=pod
		# only keep single base mutations -- exclude DNVs etc.
		if ( $tfstart != $tfend )
		{
			next;
		}
=cut
		my $plusmutationval = 0;
		my $minusmutationval = 0;
		my $bothval = 0;
		my $tfmidpoint;
		if ( $tfstrand eq "+" )
		{
	        	$tfmidpoint = $chr . ":" . $tfstart . "_" . $tfend;

			# don't repeat analysis of TF with same midpoint as previous TF
			if ( exists $prevtfmids{$tfmidpoint} && $prevtfmids{$tfmidpoint} == 1 )
			{
				next;
			}
	
			$prevtfmids{$tfmidpoint} = 1;
	
			for ( my $j = $tfstart; $j <= $tfend; $j++ )
			{	
				if ( exists $plusreads{$j} )
				{
					$plusmutationval += $plusreads{$j};
					$bothval += $plusreads{$j};	
				}
			
	        	        if ( exists $minusreads{$j} )
	                	{
	                                $minusmutationval += $minusreads{$j};
					$bothval += $minusreads{$j};
	                        }
			}
		}
                elsif ( $tfstrand eq "-" )
                {
                        $tfmidpoint = $chr . ":" . $tfstart . "_" . $tfend;

                        # don't repeat analysis of TF with same midpoint as previous TF
                        if ( exists $prevtfmids{$tfmidpoint} && $prevtfmids{$tfmidpoint} == 1 )
                        {
                                next;
                        }

                        $prevtfmids{$tfmidpoint} = 1;

                        for ( my $j = $tfend; $j >= $tfstart; $j-- )
                        {
				# switch strands because motif is on minus strand
                                if ( exists $plusreads{$j} )
                                {
                                        $minusmutationval += $plusreads{$j};
					$bothval += $plusreads{$j};
                                }

                                if ( exists $minusreads{$j} )
                                {
                                        $plusmutationval += $minusreads{$j};
					$bothval += $minusreads{$j};
                                }
                       }
		}
		else
		{
			die "No strand information for $tfmotif\n";
		}

		print OUT "$chr:$bedstart-$tfend\t$tfname\t$tfmotif\t$bothval\n";
		print BED "$chr\t$bedstart\t$tfend\t$chr:$bedstart-$tfend\t$tfmotif\t$tfstrand\n";
	}
}
