package Circos::ColorTransform;

=pod

=head1 NAME

Circos::ColorTransform - Transform a color

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
									transform_color
							 );

use Carp qw( carp confess croak );
use FindBin;
#use File::Basename;
#use File::Spec::Functions;
#use File::Temp qw(tempdir);
use GD;
use Memoize;
use Math::Round;
#use Math::VecStat qw(min max);
use Params::Validate qw(:all);
#use Regexp::Common;
#use Storable;
#use Math::Trig;
#use Sys::Hostname;

#use Time::HiRes qw(gettimeofday tv_interval);
#use List::Util qw( max min );

use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/lib";

use Circos::Configuration;
use Circos::Constants;
use Circos::Debug;
use Circos::Error;
use Circos::Utils;
use Circos::LCH;

our $transform_table;

################################################################
# Transform the RGB coordinates of a color.
sub transform_color {
	my ($rgb,$name,$method) = @_;
	if(defined $transform_table->{$name}{$method}) {
		printinfo("random","lookup",$name);
		return $transform_table->{$name}{$method};
	}
	if($method =~ /rand/) {
		printinfo("random","gen",$name);
		printinfo(@$rgb);
		$rgb->[0] = int(rand(256));
		$rgb->[1] = int(rand(256));
		$rgb->[2] = int(rand(256));
		printinfo(@$rgb);
		printinfo();
	} elsif ($method =~ /wmmn/) {
    if($rgb->[0] == 255 && $rgb->[1] == 255 && $rgb->[2] == 255) {
      $rgb = [0,0,0];
    } else {
      my $d1  = deltae([255,0,0],$rgb);
      my $d2  = deltae([0,255,0],$rgb);
      my @lch = rgb_to_lch(@$rgb);
      if($d1 < $d2) {
        $lch[0] = remap($d1,0,200,0,150);
        $lch[1] = remap($d1,0,200,150,0);
        $lch[2] = 40;
      } else {
        $lch[0] = remap($d2,0,200,0,150);
        $lch[1] = remap($d2,0,200,150,0);
        $lch[2] = 150;
      }
      $rgb = [lch_to_rgb(@lch)];
    }
	}
	return $transform_table->{$name}{$method} = $rgb;
}

1;
