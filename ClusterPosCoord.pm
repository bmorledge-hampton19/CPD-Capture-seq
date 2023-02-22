#!/usr/bin/perl

use strict;
use warnings;

package ClusterPosCoord;

sub new
{
	my ($class) = @_;
	
	my $self = bless {}, $class;

	#my $tf_file_name = "../hayward_exons_recurrentmuts_mela.bed";
	#my $tf_file_name = "../ENCODE.tf.bound.union_known_sorted_proximal_DHS_etsmotif_GGatpos23.bed";
	#my $tf_file_name = "../captureseqprobes_with600plus_sorted.bed";
	#my $tf_file_name = "../Compiled_capturearray_simple.bed";
	my $tf_file_name = "../Unique_ets_known_34dipy_DHS_revshr_Capture.bed";
	#my $tf_file_name = "../Unique_ets_known_NOT34dipy_DHS_revshr_Capture.bed";

	open( GENE, $tf_file_name) || die "Couldn't open file $tf_file_name\n";
	
	my %tf;

	# Remove TFBS that overlap with blacklisted regions (100 bp window on each side of midpoint)
	my $blist_filename = "../ENCODE.tf.bound.union_known_sorted_BLIST.bed";
	open ( BLIST, $blist_filename ) || die "Couldn't open Blacklist file $blist_filename\n";
	my %blacklist;
	while ( my $line = <BLIST> )
	{
		chomp $line;
		$blacklist{$line} = 1;
	}
	
	while( my $line = <GENE> )
	{
		chomp($line);
		if ( $line eq "" )
		{	next;	}
		# skip blacklisted lines
		if ( exists $blacklist{$line} && $blacklist{$line} == 1 )
		{	
			print STDERR "Blacklisted TFBS: $line\n";
			next;
		}

		my @fields = split /\t/, $line;
		push @{$tf{$fields[0]}}, [$fields[1], $fields[2], $fields[3], $fields[5]];
	}

	close ( GENE );
	$self->{'tf'} = \%tf;
	$self->{'tf_file_name'} = $tf_file_name;
	return $self;	
}

sub get_tf_boundaries
{
	my ($self) = @_;

	return %{$self->{'tf'}};

}

sub get_tf_filename
{
        my ($self) = @_;

        return $self->{'tf_file_name'};
}

1;	
