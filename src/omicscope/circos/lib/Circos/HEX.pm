package Circos::HEX;

=pod

=head1 NAME

Circos::HSV - HSV to RGB conversion

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
									hex_to_rgb
							 );

use Carp qw( carp confess croak );
use Digest::MD5 qw(md5_hex);
use Math::Round;
use Math::VecStat qw(min max);
use Math::Trig;

use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/lib";

use POSIX qw(floor);

use Circos::Constants;
use Circos::Debug;
use Circos::Error;
use Circos::Utils;

sub hex_to_rgb {
	my $str = shift;
	my @rgb = unpack 'C*', pack 'H*', $str;
	return @rgb;
}

1;
