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
print OUT "Relative Position to TF binding site midpoint";
for (my $i = (2 * $window ); $i >= 0; $i-- )
{
        my $pos = $window - $i;
        print OUT "\t$pos";
}
print OUT "\n";

my @pluscpdval = ();
my @minuscpdval = ();
my @bothval = ();
my @pluscpdcount = ();
my @minuscpdcount = ();

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
		# increase start position by 1, since wig files (data files) are one-based, and bed files (position files) are zero-based (only start, not end coordinate)
		$tfstart++;
		my $bedstart;
		my $bedend;

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
			$bedstart = $winstart - 1;
                        $bedend = $winend;
			for ( my $j = $winstart; $j <= $winend; $j++ )
			{	
				my $relpos = $j - $winstart;
				$pluscpdcount[$relpos]++;
				$minuscpdcount[$relpos]++;
				if ( exists $plusreads{$j} )
				{
					$pluscpdval[$relpos] += $plusreads{$j};
					$bothval[$relpos] += $plusreads{$j};	
				}
			
	        	        if ( exists $minusreads{$j} )
	                	{
	                                $minuscpdval[$relpos] += $minusreads{$j};
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
			$bedstart = $winstart - 1;
			$bedend = $winend;
                        for ( my $j = $winend; $j >= $winstart; $j-- )
                        {
                                my $relpos = $winend - $j;

                                $pluscpdcount[$relpos]++;
                                $minuscpdcount[$relpos]++;

				# switch strands because motif is on minus strand
                                if ( exists $plusreads{$j} )
                                {
                                        $minuscpdval[$relpos] += $plusreads{$j};
					$bothval[$relpos] += $plusreads{$j};
                                }

                                if ( exists $minusreads{$j} )
                                {
                                        $pluscpdval[$relpos] += $minusreads{$j};
					$bothval[$relpos] += $minusreads{$j};
                                }
                        }
		}
		else
		{
			die "No strand information for $tfmotif\n";
		}

		#print OUT "$chr:$bedstart-$bedend\t${tfmotif}_${tfmidpoint}_${tfstrand}\tBothStrands";
		print BED "$chr\t$bedstart\t$bedend\t$chr:$bedstart-$bedend\t$tfmotif\t$tfstrand\n";

	}
}

#print header

print OUT "CPD read counts";
print OUT "\nPlus Strand";
for (my $i = 0; $i < scalar @pluscpdval; $i++ )
{
        print OUT "\t$pluscpdval[$i]";
}
print OUT "\nMinus Strand";
for (my $i = 0; $i < scalar @minuscpdval; $i++ )
{
        print OUT "\t$minuscpdval[$i]";
}
print OUT "\n\nTotal nucleotide counts";
print OUT "\nPlus strand:";
for (my $i = 0; $i < scalar @pluscpdcount; $i++ )
{
        print OUT "\t$pluscpdcount[$i]";
}
print OUT "\nMinus strand:";
for (my $i = 0; $i < scalar @minuscpdcount; $i++ )
{
        print OUT "\t$minuscpdcount[$i]";
}

print OUT "\n\nAverage CPDs:";
print OUT "\nPlus strand:";
for (my $i = 0; $i < scalar @pluscpdcount; $i++ )
{
        if ( $pluscpdcount[$i] != 0 )
        {
                my $mean = 1.0 * $pluscpdval[$i] / $pluscpdcount[$i];
                print OUT "\t$mean";
        }
        else
        {
                print OUT "\t ";
        }
}
print OUT "\nMinus strand:";
for (my $i = 0; $i < scalar @minuscpdcount; $i++ )
{
        if ( $minuscpdcount[$i] != 0 )
        {
                my $mean = 1.0 * $minuscpdval[$i] / $minuscpdcount[$i];
                print OUT "\t$mean";
        }
        else
        {
                print OUT "\t ";
        }
}

print OUT "\n\nStrand Average:";
print OUT "\nCPD read counts";
for (my $i = 0; $i < scalar @pluscpdval || $i < scalar @minuscpdval; $i++ )
{
        my $tempsum = $pluscpdval[$i] + $minuscpdval[$i];
	if ( $tempsum != $bothval[$i] )
	{	die "mismatched counts\n";	}
        print OUT "\t$tempsum";
}
print OUT "\nTotal nucleotide counts:";
for (my $i = 0; $i < scalar @pluscpdcount; $i++ )
{
        if ($pluscpdcount[$i] != $minuscpdcount[$i] )
        {
                die "count of nucleotides on plus and minus strand don't match\n";
        }
        else
        {
                print OUT "\t$pluscpdcount[$i]";
        }
}
print OUT "\nAverage CPDs";
for (my $i = 0; $i < scalar @pluscpdval || $i < scalar @minuscpdval; $i++ )
{
        my $tempcpdsum = $pluscpdval[$i] + $minuscpdval[$i];
        my $tempnucsum;
        if ($pluscpdcount[$i] != $minuscpdcount[$i] )
        {
                die "count of nucleotides on plus and minus strand don't match\n";
        }
        else
        {
                $tempnucsum = $pluscpdcount[$i];
        }
        my $mean = 1.0 * $tempcpdsum / $tempnucsum;
        print OUT "\t$mean";
}
print OUT "\n\nPosition Relative to Motif Midpoint";
for (my $i = 0; $i < scalar @pluscpdcount; $i++ )
{
        my $pos = $i - $window;
        print OUT "\t$pos";
}
print OUT "\n\n";

