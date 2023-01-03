package Circos::Colors;

=pod

=head1 NAME

Circos::Colors - Color handling for Circos

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
									allocate_colors
									allocate_color
									color_to_list
									find_transparent
									rgb_color
									rgb_color_opacity
									rgb_color_transparency
									rgb_to_color
									fetch_color
									aa_color
							 );

use Carp qw( carp confess croak );
use Digest::MD5 qw(md5_hex);
use FindBin;
use File::Basename;
use File::Spec::Functions;
use File::Temp qw(tempdir);
use List::MoreUtils qw( uniq );
use GD;
use Memoize;
use Math::Round;
use Math::VecStat qw(min max);
use Params::Validate qw(:all);
use Regexp::Common;
use Storable;
use Math::Trig;
use Sys::Hostname;

#use Time::HiRes qw(gettimeofday tv_interval);
#use List::Util qw( max min );

use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/lib";

use POSIX qw(pow);

use Circos::Configuration;
use Circos::Constants;
use Circos::Debug;
use Circos::Error;
use Circos::Image;
use Circos::Utils;
use Circos::HSV;
use Circos::LCH;
use Circos::RGB;
use Circos::ColorValidate;

for my $f ( qw (rgb_color_opacity rgb_color_transparency ) ) {
  memoize($f);
}

# -------------------------------------------------------------------
sub allocate_colors {

	my $image            = shift;

	my $allocated_colors = 0;
	my $colors           = {};

	# scan the <colors> block and first allocate all colors
	# specified as
	#   r,g,b   or rgb(r,g,b)
  #   r,g,b,a or rgb(r,g,b,a) or rgba(r,g,b,a)
	#
	#   lch(l,c,h)
  #   lch(l,c,h,a)
  #   hsv(h,s,v)
  #   hsv(h,s,v,a)
  #   hex(59311c)
	#
	# resolution of name lookups or lists is deferred until later
	
	start_timer("colordefinitions");

	# First, allocate color that are defined based on color space coordinates. 

	for my $color_name ( sort keys %{ $CONF{colors} } ) {
		if (ref $CONF{colors}{$color_name} eq "ARRAY") {
	    my @unique_definitions = uniq @{$CONF{colors}{$color_name}};
	    if (@unique_definitions == 1) {
				my $txt = join($SPACE,@unique_definitions);
				printwarning("The color [$color_name] has multiple identical definitions: $txt");
				$CONF{colors}{$color_name} = $unique_definitions[0];
	    } else {
				fatal_error("color","multiple_defn",$color_name,
										join($NEW_LINE, map { " $_" } @unique_definitions));
	    }
		} elsif ( my $cref = ref $CONF{colors}{$color_name}) {
	    fatal_error("color","malformed_structure",$color_name,$cref);
		}
		my $color_definition = $CONF{colors}{$color_name};
		next if $color_definition =~ /\|/;
		if(my @rgb = validate_color($color_definition,$color_name,0)) {
	    allocate_color($color_name,\@rgb,$colors,$image);
		}
	}

	stop_timer("colordefinitions");

	# now resolve name lookups
	start_timer("colorlookups");
	for my $color_name ( sort keys %{ $CONF{colors} } ) {
		my $color_definition = lc $CONF{colors}{$color_name};
		# if this color has already been allocated, skip it
		next if exists $colors->{$color_name};
		my %lookup_seen;
		while ( exists $CONF{colors}{$color_definition} ) {
	    printdebug_group("color","colorlookup",$color_definition);
	    if ($lookup_seen{$color_definition}++) {
				fatal_error("color","circular_defn",$color_definition,$CONF{color}{$color_definition});
	    }
	    $colors->{$color_name} = $colors->{$color_definition};
	    printdebug_group("color","colorlookupassign",
											 $color_name,$color_definition,
											 $CONF{colors}{$color_definition});
	    $color_definition = $CONF{colors}{$color_definition};
		}
	}
	stop_timer("colorlookups");

	# automatic transparent colors
	start_timer("colortransparency");
	create_transparent_colors($colors,$image);
	stop_timer("colortransparency");

	# now resolve lists - employ caching since this can be slow (2-5 seconds);

	goto SKIPLISTS if defined_and_zero(fetch_conf("color_lists_use"));

	start_timer("colorlists");
	my $hostname = hostname;
	my $user     = $ENV{USERNAME} ? $ENV{USERNAME} . $PERIOD : $EMPTY_STRING;
	my $cache_file;
	my $cache_file_root = sprintf("%s.%s.%sdat",
																Circos::Configuration::fetch_configuration("color_cache_file") || "circos.colorlist",
																$hostname,
																$user);				  
	if (my $cache_file_dir = Circos::Configuration::fetch_configuration("color_cache_dir")) {
		$cache_file = catfile($cache_file_dir,$cache_file_root);
	} else {
		# use File::Temp to temporarily create a directory and use this
		# to figure out the system's temporary directory root (e.g. /tmp)
		my $cache_dir = tempdir();
		rmdir($cache_dir);
		if (! $cache_dir) {
	    fatal_error("io","temp_dir_not_created","color_cache_dir");
		}
		my $cache_dir_root = dirname($cache_dir);
		printdebug_group("cache","temporary file dir",$cache_dir_root);
		$cache_file = catfile($cache_dir_root,$cache_file_root);
	}
	my $allocated_color_list = [keys %$colors];
	my $list_cache;
	my $cache_ok;
	my $is_cache_static = Circos::Configuration::fetch_configuration("color_cache_static");
	my $rebuild_cache   = Circos::Configuration::fetch_configuration("color_cache_rebuild");
	if ($rebuild_cache) {
		printdebug_group("cache","colorlist cache rebuild forced");
	} elsif (-e $cache_file) {
		start_timer("colorcache");
		printdebug_group("cache","colorlist cache",$cache_file,"found");
		if ($is_cache_static || -M $cache_file < -M $CONF{configfile}) {
			printdebug_group("cache","colorlist cache",$cache_file,"useable - static or more recent than configfile");
			# cache file younger than config file, read cache
			eval {
		    $list_cache = retrieve($cache_file);
			};
			if ($@) {
		    printwarning("Problem reading color cache file $cache_file");
		    $cache_ok = 0;
			} else {
		    printdebug_group("cache","colorlist cache",$cache_file,"read in");
		    my $target_hash = Digest::MD5::md5_hex(join("", sort keys %{$CONF{colors}}));
		    if ($list_cache->{colorhash} eq $target_hash) {
					printdebug_group("cache","color list hash",$target_hash,"matches that of cache file - using cache file");
					$cache_ok = 1;
		    } elsif ($is_cache_static) {
					printdebug_group("cache","color list hash",$target_hash,"doesn't match that of cache file - using cache anyway because it is static");
					$cache_ok = 1;
		    } else {
					printdebug_group("cache","color list hash",$target_hash,"does not match - colors changed? - recomputing file");
		    }
			}
		} else {
			printdebug_group("cache","colorlist cache",$cache_file,"older than configfile - recreating cache");
		}
		stop_timer("colorcache");
	} else {
		printdebug_group("cache","colorlist cache",$cache_file,"not found");
	}
	if (! $cache_ok) {
		# create cache
		$list_cache->{colorhash} = Digest::MD5::md5_hex(join("", sort keys %{$CONF{colors}}));
		printdebug_group("cache","creating colorlist cache, hash",$list_cache->{colorhash});
		for my $color_name ( sort keys %{ $CONF{colors} } ) {
	    # skip if this color has already been allocated
	    next if exists $colors->{$color_name};
	    my @color_definitions = str_to_list($CONF{colors}{$color_name});
	    my @match_set;
	    for my $color_definition (@color_definitions) {
				# do a very quick match to narrow down the colors with fast grep()
				my $rx = $color_definition;
				if ($rx =~ /rev\((.+)\)/) {
					$rx  = $1;
				}
				my @early_matches = grep($_ =~ /$rx/, @$allocated_color_list);
				my @matches;
				# now do a full match, including sorting results
				if (@early_matches) {
					@matches = sample_list($color_definition,\@early_matches); #$allocated_color_list);
				}
				if (! @matches) {
					fatal_error("color","bad_name_in_list",$color_name,$color_definition);
				}
				push @match_set, @matches;
	    }
	    $list_cache->{list2color}{$color_name} = \@match_set;
	    printdebug_group("color","colorlist",$color_name,@match_set);
		}
		# store cache
		my $create_cache_file = Circos::Configuration::fetch_configuration("color_cache_create");
		if ( ! defined $create_cache_file || $create_cache_file ) {
	    eval { 
				printdebug_group("cache","writing to colorlist cache file [$cache_file]");
				store($list_cache,$cache_file);
	    };
		} else {
	    printdebug_group("cache","skipping creating cache file [$cache_file]");
		}
		if ($@) {
	    printwarning("Could not write to color list cache file $cache_file - store() gave error");
	    printinfo($@);
		} elsif ($create_cache_file) {
	    if (-e $cache_file) {
				printdebug_group("cache","wrote to colorlist cache file [$cache_file]");
	    } else {
				printwarning("Could not find the cache file we supposedly just created $cache_file");
	    }
		}
	}
	for my $color (keys %{$list_cache->{list2color}}) {
		$colors->{$color} = $list_cache->{list2color}{$color};
		push @$allocated_color_list, $color;
	}
	stop_timer("colorlists");
 SKIPLISTS:
	return $colors;
}

################################################################
# Returns the opacity of a color, based on its name.
#
# Colors with a trailing _aNN have a transparency level in the range
# 0..auto_alpha_steps.
#
# Otherwise, color space coordinates are checked for the fourth coordinate
#
# 255,0,0,0.5
# rgb(255,0,0,0.5)
# lch(50,50,20,0.5)
# hsv(0,1,1,0.5)

sub rgb_color_opacity {
	my $color = shift;
	my $alpha = 1; # default
	return $alpha if ! defined $color;
	if ( $color =~ /(.+)_a(\d+)$/ ) {
		unless ( $CONF{image}{auto_alpha_colors} && $CONF{image}{auto_alpha_steps} ) {
	    die "you are trying to process a transparent color ($color) ",
				"but do not have auto_alpha_colors or auto_alpha_steps defined";
		}
		my $color_root = $1;
		$alpha         = $2 / (1+$CONF{image}{auto_alpha_steps});
	} elsif ($color =~ /\(?\s*(?:$RE{num}{real}\s*,\s*){3}($RE{num}{real})\s*\)?/) {
		$alpha = $1;
		if(! is_number($alpha,"real",0,0,1)) {
			fatal_error("color","bad_alpha",$color,$alpha);
		}
	}
	return $alpha;
}

################################################################
sub rgb_color_transparency {
  my $color = shift;
  $color    = lc $color;
  return 1 - rgb_color_opacity($color);
}

################################################################
sub allocate_color {
	my ($name,$definition,$colors,$image) = @_;
	my @rgb = ref $definition eq "ARRAY" ? @$definition : split(",",$definition);
	my $idx;
	printdebug_group("color","allocate_color","rgb",@rgb);
	if ( @rgb == 3 ) {
		if ($name =~ /.+_a\d+$/) {
			fatal_error("color","reserved_name_a",$name,$definition);
		}
		eval {
			my $color_index = $image->colorExact(@rgb);
			if ( $color_index == -1 ) {
		    $colors->{$name} = $image->colorAllocate(@rgb);
			} else {
		    $colors->{$name} = $color_index;
			}
		};
		printdebug_group("color","allocate_color",@rgb,$image->colorExact(@rgb));
		if ($@) {
			fatal_error("color","cannot_allocate",$name,$definition,$@);
		}
	} elsif ( @rgb == 4 ) {
		my $alpha = $rgb[3];
		if (! is_number($alpha,"real",0,0,1)) {
			fatal_error("color","bad_alpha",sprintf("%s [%s]",join($COMMA,@rgb),$name),$alpha);
		}
		# gd's colorAllocateAlpha(r,g,b,a) uses a=0 opaque and a=127 transparent
		$rgb[3] = (1-$alpha)*127;
		#$rgb[3] = round($rgb[3]);
		eval {
			printdebug_group("color","allocate_color","gdrgb",@rgb);
			$colors->{$name} = $image->colorAllocateAlpha(@rgb);
			#printinfo($name,$colors->{$name},@rgb);
		};
		if ($@) {
			fatal_error("color","cannot_allocate",$name,$definition,$@);
		}
	}
	printdebug_group("color","allocate_color","idx",$colors->{$name},$name,@rgb,"now have",int(keys %$colors),"colors");
}

################################################################
# Automatically allocate colors with alpha values, if asked for,
# and only for colors which don't already have transparency.
#
# The number of steps is determined by auto_alpha_steps in the
# <image> block
#
# Colors with alpha values have names COLOR_aN for N=1..num_steps
# The alpha value (0,1) 0=transparent 1=opaque for step i
#
# 1-i/(num_steps+1)
#
# For example, if the number of steps is 5, then for the color
# chr19=153,0,204, the follow additional 5 colors will be
# allocated (see full list in lines with 'auto_alpha_color' with -debug).
#
# Now add automatic transparency levels to all the defined colors
# using _aN suffix
#
sub create_transparent_colors {
	my ($colors,$image) = @_;
	return unless fetch_conf("image","auto_alpha_colors");
	my $nsteps = fetch_conf("image","auto_alpha_steps");
	for my $color_name (keys %$colors) {
		# if this color is already transparent, skip it
		my $alpha = rgb_color_opacity($color_name);
		# if color already has transparency, do not generate auto-alpha colors
		next if $alpha < 1;
		my @rgb = $image->rgb($colors->{$color_name});
		# provide _a0 synonym
		$colors->{ sprintf("%s_a0",$color_name) } = $colors->{ $color_name };
		for my $i (1..$nsteps) {
			my $alpha            = 1-$i/($nsteps+1);
	    my $color_name_alpha = sprintf("%s_a%d",$color_name,$i);
	    printdebug_group("color","allocate","auto_alpha_color",$color_name_alpha,@rgb,$alpha);
	    allocate_color($color_name_alpha,[@rgb,$alpha],$colors,$image);
		}
	}
}

################################################################
# Return the RGB coordinates for a color based on name.
sub rgb_color {
  my $color = shift;
  return undef if ! defined $color;
  $color = lc $color;
  if ( $color =~ /(.+)_a(\d+)/ ) {
		my $color_root = $1;
		return rgb_color($color_root);
  } else {
		return undef unless defined $color;
		my @rgb;
		if ( defined $COLORS->{$color} ) {
			@rgb = $IM->rgb( $COLORS->{$color} );
		} else {
			my $cnew = fetch_color( $color, $COLORS, $IM );
			@rgb = $IM->rgb( $cnew );
		}
		return @rgb;
		my $colordef  = $COLORS->{$color};
		if ($COLORS->{$colordef}) {
			return rgb_color($colordef);
		}
		@rgb = split( $COMMA, $colordef );
		return @rgb;
  }
}

################################################################
# Given an RGB value, return the color name
sub rgb_to_color_name {
	my @args = @_;
	my ($r,$g,$b,$a,$no_error);
	if(@args == 3) {
		($r,$g,$b) = @args;
	} elsif (@args == 4) {
		($r,$g,$b,$no_error) = @args;
	} elsif (@args == 5) {
		($r,$g,$b,$a,$no_error) = @args;
	} else {
		die "Bad number of arguments.";
	}
	# If alpha is found, then assume the color has not been defined
	return;
	for my $color (keys %$COLORS) {
		next if $color =~ /_a\d+$/;
		my @crgb = $IM->rgb( fetch_color($color) );
		if ($r == $crgb[0] &&	$g == $crgb[1] &&	$b == $crgb[2]) {
	    return $color;
		}
	}
	return if $no_error;
	fatal_error("color","bad_rgb_lookup",$r,$g,$b);
}

################################################################
# Fetch the color. 

sub fetch_color {
	my ($color_name,$color_table,$im) = @_;
	$color_table ||= $COLORS;
	$im          ||= $IM;
	start_timer("colorfetch");
	if (exists $color_table->{$color_name}) {
		stop_timer("colorfetch");
		printdebug_group("color","fetch",$color_name,$color_table->{$color_name});
		return $color_table->{$color_name};
	} elsif ($color_table->{lc $color_name}) {
		my $lc_color = lc $color_name;
		printwarning("Circos colors should be lowercase. You have asked for color [$color_name] and it was interpreted as [$lc_color]");
		stop_timer("colorfetch");
		return $color_table->{lc $color_name};
	} elsif (my @rgb = validate_color($color_name,$color_name,0)) {
		allocate_color($color_name,\@rgb,$color_table,$im);
		return $color_table->{$color_name};
	} elsif (my $default_color = fetch_conf("default_color")) {
		return $color_table->{$default_color};
	} else {
		fatal_error("color","undefined",$color_name);
	}
}

################################################################
# Return a color object, depending on whether antialiasing is set or not.
#
# Antialiasing in GD works only for opaque colors.
sub aa_color {
	my ($color_name,$im,$imcolors) = @_;
	my $color = fetch_color($color_name,$imcolors);
	if (not_defined_or_one(fetch_conf("anti_aliasing")) && rgb_color_opacity($color_name) == 1) {
		$im->setAntiAliased($color);
		return gdAntiAliased;
	} else {
		$color;
	}
}

################################################################
# Find the first unallocated RGB color that we can use
# as the transparent color. GD is stupid this way -- you actually
# need to define a color first and then make it transparent :/ 
sub find_transparent {
	my @rgb = (0,0,0);
	my $idx = 0;
	my $color;
	do {
		$rgb[ $idx % 3]++;
		$idx++;
		eval { $color = rgb_to_color_name(@rgb,1) };
	} while ($color);
	return @rgb;
}

################################################################
# 
sub color_to_list {
	my $color = shift;
	return if ! defined $color;
	my @color_names = split(/[\s+,]+/,$color);
	my @colors;
	for my $color_name (@color_names) {
		if (ref $COLORS->{$color_name} eq "ARRAY") {
	    push @colors, @{$COLORS->{$color_name}};
		} else {
	    push @colors, $color_name;
		}
	}
	return @colors;
}

1;
