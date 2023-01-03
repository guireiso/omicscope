package Circos::RGB;

=pod

=head1 NAME

Circos::RGB - RGB color routines

=head1 SYNOPSIS

This module is not meant to be used directly.

=head1 DESCRIPTION

Circos is an application for the generation of publication-quality,
circularly composited renditions of genomic data and related
annotations.

Circos is particularly suited for visualizing alignments, conservation
and intra and inter-chromosomal relationships. However, Circos can be
used to plot any kind of 2D data in a circular layout - its use is not
limited to genomics. Circos' use of lines to relate position pairs
(ribbons add a thickness parameter to each end) is effective to
display relationships between objects or positions on one or more
scales.

All documentation is in the form of tutorials at L<http://www.circos.ca>.

=cut

# -------------------------------------------------------------------

use strict;
use warnings;

use base 'Exporter';
our @EXPORT = qw(
									rgb_to_rgb255
							 );

use Carp qw( carp confess croak );
use Digest::MD5 qw(md5_hex);
use Math::Round;
use Math::VecStat qw(min max);
use Math::Trig;

use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/lib";

use POSIX qw(pow);

use Circos::Constants;
use Circos::Debug;
use Circos::Error;
use Circos::Utils;

################################################################
# Applies only to first three coordinates
sub rgb_to_rgb255 {
	my $rgb   = shift;
	die "RGB coordinates malformed ".join(",",@$rgb) unless @$rgb == 3 || @$rgb == 4;
	my $alpha = pop @$rgb if @$rgb == 4;
	# make sure values are in [0,1]
	my @rgb    = map { put_between($_,0,1) } @$rgb;
	my @rgb255 = map { round( 255 * $_) } @rgb;
	push @rgb255, $alpha if defined $alpha;
	return @rgb255;
}

1;
