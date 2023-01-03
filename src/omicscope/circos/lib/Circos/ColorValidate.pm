package Circos::ColorValidate;

=pod

=head1 NAME

Circos::ColorValidate - Validate color syntax

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
									validate_color
									validate_rgb
									validate_hsv
									validate_lch
									validate_hex
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
use Regexp::Common;
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
use Circos::HEX;
use Circos::HSV;
use Circos::RGB;
use Circos::ColorTransform;

################################################################
# Examine a color definition and return RGB values for it,
# if the definition is well-formed.
#
# If asked for, also transform the color.
#
sub validate_color {
	my ($color_definition,$color_name,$strict) = @_;
	my @rgb;
	if (my @hsv  = validate_hsv($color_definition,$strict)) {
		my @rgb255 = hsv_to_rgb(@hsv);
		@rgb       = rgb_to_rgb255(\@rgb255);
		printdebug_group("color","parsing_color hsv",$color_definition,"rgb",@rgb);
	} elsif (my @lch = validate_lch($color_definition,$strict)) {
		@rgb = lch_to_rgb(@lch);
	} elsif (my ($hex,$alpha) = validate_hex($color_definition,$strict)) {
		@rgb = hex_to_rgb($hex);
		push @rgb, $alpha if defined $alpha;
		printdebug_group("color","parsing_color hex",$color_definition,"rgb",@rgb);
	} elsif (@rgb = validate_rgb($color_definition,$strict)) {
		printdebug_group("color","parsing_color rgb",$color_definition);
	}
	if(@rgb && validate_rgb_components(\@rgb,$strict,$color_definition)) {
		# we're good, transform below if asked for
	} else {
		fatal_error("color","malformed_rgb",$color_definition) if $strict;
		return;
	}

	my $rgb_t;
	if(defined fetch_conf("randomcolor")) {
		if(! grep($color_name eq $_,split($COMMA,fetch_conf("randomcolor")))) {
			$rgb_t = transform_color(\@rgb,$color_definition,"random");
		}
	} elsif (defined fetch_conf("wmmn")) {
		if(! grep($color_name eq $_,split($COMMA,fetch_conf("wmmn")))) {
			$rgb_t = transform_color(\@rgb,$color_definition,"wmmn");
		}
	}
	return $rgb_t ? @$rgb_t : @rgb;

}

################################################################
# Validate HEX color
# HEX or hex(HEX)
# e.g. FFAA06 ffaa06 hex(FFAA06) hex(ffaa06)
#
# with optional transprency as second argument
#
# FFAA06,0.25
# hex(FFAA06,0.25)
sub validate_hex {
	my ($definition,$strict) = @_;
	$strict = 1 if ! defined $strict;
	my ($str,$alpha);
	my $rx = qr/([0-9A-F]{6})(?:\s*,\s*($RE{num}{real}))?/i;
	if( $definition =~ /hex?\s*\(\s*$rx\s*\)/i ) {
		($str,$alpha) = ($1,$2);
		fatal_error("color","malformed_hex",$definition) if $str !~ /^$rx$/i;
	} elsif ($definition =~ /^\s*$rx\s*$/) {
		($str,$alpha) = ($1,$2);
	}
	return ($str,$alpha) if defined $str;
	fatal_error("color","malformed_rgb",$definition) if $strict;
	return;
}

################################################################
# Validate rgb 3 or 4 value list.
sub validate_rgb {
	my ($definition,$strict) = @_;
	$strict = 1 if ! defined $strict;
	my @rgb;
	if ( ref $definition eq "ARRAY") {
		@rgb = @$definition;
	} elsif ( $definition =~ /rgba?\s*\(\s*([\d.,]+)\s*\)/i ) {
		@rgb = split(/\s*,\s*/,$1);
	} elsif ( $definition =~ /,/ ) {
		@rgb = split(/\s*,\s*/,$definition);
	}
	# if @rgb is defined then we have found what looks like a possible RGB color
	# set $strict so that if the string is malformed, error is produced
	$strict = 1 if @rgb;
	if (@rgb && validate_rgb_components(\@rgb,$strict,$definition)) {
		return @rgb;
	}
	fatal_error("color","malformed_rgb",$definition) if $strict;
	return;
}

################################################################
# Verify that a list is allowable RGB components
#
# r,g,b[,a] in rgb:[0,255] and a:[0,1]
sub validate_rgb_components {
	my ($rgb,$strict,$str) = @_;
  my $n = grep(defined $_,@$rgb);
	my ($r,$g,$b,$a) = @$rgb;
	return unless $n == 3 || $n == 4;
	for ($r,$g,$b) {
		if(! is_number($_,"real",$strict,0,255)) {
			fatal_error("color","bad_color_component","rgb",$str,$_,0,255);
		}
	}
	if (defined $a && ! is_number($a,"real",$strict,0,1)) {
		fatal_error("color","bad_alpha",$str,$a);
	}
	return 1;
}

################################################################
# Verify that a string is a valid HSV color definition
sub validate_hsv {
	my ($definition,$strict) = @_;
	$strict = 1 if ! defined $strict;
	my @hsv;
	if ( ref $definition eq "ARRAY" ) {
		@hsv = @$definition;
	} elsif ( $definition =~ /hsv\s*\(\s*([\d.,]+)\s*\)/i ) {
		@hsv = split(/\s*,\s*/,$1);
	}
	$strict = 1 if @hsv;
	if (@hsv && validate_hsv_components(\@hsv,$strict,$definition)) {
		return @hsv;
	}
	fatal_error("color","malformed_hsv",$definition) if $strict;
	return;
}

################################################################
# Verify that a list is an allowable RGB or RGBA list.
#
# A = alpha (0..1), 0 = transparent, 1 = opaque
sub validate_hsv_components {
  my ($hsv,$strict,$str) = @_;
  my $n    = grep(defined $_,@$hsv);
	my ($h,$s,$v,$a) = @$hsv;
	return unless $n == 3 || $n == 4;
	for ([$h,0,360],[$s,0,1],[$v,0,1]) {
		if(! is_number($_->[0],"real",$strict,$_->[1],$_->[2])) {
			fatal_error("color","bad_color_component","hsv",$str,@$_);
		}
	}
	if (defined $a && ! is_number($a, "real", $strict, 0,1)) {
		fatal_error("color","bad_alpha",$a);
	}
	return 1;
}

################################################################
# Validate LCH

sub validate_lch {
	my ($definition,$strict) = @_;
	$strict = 1 if ! defined $strict;
	my @lch;
	if ( ref $definition eq "ARRAY") {
		@lch = @$definition;
	} elsif ( $definition =~ /lch\s*\(\s*([-\d.,]+)\s*\)/i ) {
		@lch = split(/\s*,\s*/,$1);
	}
	$strict = 1 if @lch;
	if (@lch && validate_lch_components(\@lch,$strict,$definition)) {
		#printinfo(@lch);
		return @lch;
	}
	fatal_error("color","malformed_lch",$definition) if $strict;
	return;
}

################################################################
# Verify that a list is an allowable LCH or LCHA list.
#
# a = alpha (0..1), 0 = transparent, 1 = opaque
sub validate_lch_components {
	my ($lch,$strict,$str) = @_;
	my $n = grep(defined $_,@$lch);
	my ($l,$c,$h,$a) = @$lch;
	return unless $n == 3 || $n == 4;
	for ([$l,0,150],[$c,0,150],[$h,-180,360]) {
		if(! is_number($_->[0],"real",$strict,$_->[1],$_->[2])) {
			fatal_error("color","bad_color_component","lch",$str,@$_);
		}
	}
	if (defined $a && ! is_number($a, "real", $strict, 0,1)) {
		fatal_error("color","bad_alpha",$str,$a);
	}
	return 1;
}

1;
