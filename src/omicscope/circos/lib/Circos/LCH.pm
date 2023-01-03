package Circos::LCH;

=pod

=head1 NAME

Circos::LCH - LCH (ab) to RGB conversion

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
									deltae
									rgb_to_lch
									lch_to_rgb
							 );

use Carp qw( carp confess croak );
use Digest::MD5 qw(md5_hex);
use Math::Round;
use Math::VecStat qw(sum min max);
use Math::Trig;

use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/lib";

use POSIX qw(pow);

use Circos::Constants;
use Circos::Debug;
use Circos::Error;
use Circos::Utils;
use Circos::RGB;

my $white = { 'D65' => [ 0.312713, 0.329016 ] }; # Daylight 6504K
my $srgb  = {
						 white_point => 'D65',
						 gamma => 'sRGB', # 2.4,
						 m     => [ [  0.4124237575757575,  0.2126560000000000,  0.0193323636363636 ], [  0.3575789999999999,  0.7151579999999998,  0.1191930000000000 ], [  0.1804650000000000,  0.0721860000000000,  0.9504490000000001 ] ],
						 mstar => [ [  3.2407109439941704, -0.9692581090654827,  0.0556349466243886 ], [ -1.5372603195869781,  1.8759955135292130, -0.2039948042894247 ], [ -0.4985709144606416,  0.0415556779089489,  1.0570639858633826 ] ],
					 };

sub deltae {
	my ($rgb1,$rgb2) = @_;
	my $lab1 = RGB_to_Lab([ map { $_/255 } @$rgb1 ]);
	my $lab2 = RGB_to_Lab([ map { $_/255 } @$rgb2 ]);
	return sqrt( sum ( map { ($lab1->[$_]-$lab2->[$_])**2 } (0..2) ) );
}

################################################################
# LCH to RGB
sub rgb_to_lch {
	my ($r,$g,$b,$alpha) = @_;
	my $lab = RGB_to_Lab([$r/255,$g/255,$b/255]);
	my $lch = Lab_to_LCHab($lab);
	push @$lch,$alpha if defined $alpha;
	return @$lch;
}

sub lch_to_rgb {
	my ($L,$C,$H,$alpha) = @_;
	my $lab = LCHab_to_Lab([$L,$C,$H]);
	my $rgb = Lab_to_RGB($lab);
	push @$rgb, $alpha if defined $alpha;
	return rgb_to_rgb255($rgb);
}

sub xyY_to_XYZ {
	my ($xyy) = @_;
	my ($x, $y, $Y) = @{$xyy};
	my ($X, $Z);
	if (! ($y == 0)) {
		$X = $x * $Y / $y;
		$Z = (1 - $x - $y) * $Y / $y;
	} else {
		$X = 0; $Y = 0; $Z = 0;
	}
	return [ $X, $Y, $Z ];
}

################################################################
# LCH to Lab
sub LCHab_to_Lab {
	my $lch = shift;
	my ($L, $C, $H) = @{$lch};
	my ($a, $b);

	$H *= $TWOPI/360;
	my $th = tan($H);
	$a = $C / sqrt( $th * $th + 1 );
	$b = sqrt($C*$C - $a*$a);

	#$H = $H - 2*pi*int($H / 2*pi); # convert H to 0..2*pi - this seems to be wrong
	if ($H < 0) { $H = $H + 2*$PI; }
	if ($H > $PI_HALF && $H < 3*$PI_HALF) { $a = - $a; }
	if ($H > $PI) { $b = - $b; }

	return [ $L, $a, $b ];
}

################################################################
# Lab to RGB
sub Lab_to_RGB {
	my $lab = shift;
	my $xyz_white = &RGB_to_XYZ([ 1.0, 1.0, 1.0 ], "sRGB");
	my $xyz = &Lab_to_XYZ($lab, $xyz_white);
	return &XYZ_to_RGB($xyz, "sRGB");
}

sub RGB_to_XYZ {
	my $rgb = shift;
	my $rgb_lin = &RGB_to_linear_RGB($rgb,$srgb);
	my $xyz = &_mult_v3_m33($rgb_lin, $srgb->{m});
	return ($xyz);
}

sub XYZ_to_RGB {
	my ($xyz, $space) = @_;
	my $rgb_lin = &_mult_v3_m33($xyz, $srgb->{mstar});
	my $rgb = &linear_RGB_to_RGB($rgb_lin, $space);
	return ($rgb);
}

sub linear_RGB_to_RGB {
	my $rgb = shift;
	my ($R, $G, $B) = @{$rgb};

	my $s = $srgb;
	if ($s->{gamma} eq 'sRGB') {
		# handle special sRGB gamma curve
		if ( abs($R) <= 0.0031308 ) { $R = 12.92 * $R; }
		else { $R = 1.055 * &_apow($R, 1/2.4) - 0.055; };
		
		if ( abs($G) <= 0.0031308 ) { $G = 12.92 * $G; }
		else { $G = 1.055 * &_apow($G, 1/2.4) - 0.055; }
		
		if ( abs($B) <= 0.0031308 ) { $B = 12.92 * $B; }
		else { $B = 1.055 * &_apow($B, 1/2.4) - 0.055; }
	} else {
		$R = &_apow($R, 1/$s->{gamma});
		$G = &_apow($G, 1/$s->{gamma});
		$B = &_apow($B, 1/$s->{gamma});
	}
	return [ $R, $G, $B ];
}

sub RGB_to_linear_RGB {
	my $rgb = shift;
	my ($R, $G, $B) = @{$rgb};

	my $s = $srgb;
	if ($s->{gamma} eq 'sRGB') # handle special sRGB gamma curve
	{
		if ( abs($R) <= 0.04045 ) { $R = $R / 12.92; }
		else { $R = &_apow( ( $R + 0.055 ) / 1.055 , 2.4 ); }

		if ( abs($G) <= 0.04045 ) { $G = $G / 12.92; }
		else { $G = &_apow( ( $G + 0.055 ) / 1.055 , 2.4 ); }

		if ( abs($B) <= 0.04045 ) { $B = $B / 12.92; }
		else { $B = &_apow( ( $B + 0.055 ) / 1.055 , 2.4 ); }
	} else {
		$R = &_apow($R, $s->{gamma});
		$G = &_apow($G, $s->{gamma});
		$B = &_apow($B, $s->{gamma});
	}
	return [ $R, $G, $B ];
}

sub RGB_to_Lab
{
	my $rgb = shift;
	my $xyz_white = &RGB_to_XYZ([ 1.0, 1.0, 1.0 ], $srgb);
	my $xyz = &RGB_to_XYZ($rgb, $srgb);
	return &XYZ_to_Lab($xyz, $xyz_white);
}

sub Lab_to_LCHab
{
	my $lab = shift;
	my ($L, $a, $b) = @{$lab};
	my ($C, $H);
	$C = sqrt( $a*$a + $b*$b );
	$H = atan2( $b, $a );
	$H = rad2deg($H);
	return [ $L, $C, $H ];
}

sub XYZ_to_Lab
{
	my ($xyz, $xyz_white) = @_;
	my ($X, $Y, $Z) = @{$xyz};
	my ($Xw, $Yw, $Zw) = @{$xyz_white};
	my ($L, $a, $b);

	my $epsilon =  0.008856;
	my $kappa = 903.3;

	my ($fx, $fy, $fz);
	my ($xr, $yr, $zr) = ( $X /  $Xw, 
						   $Y /  $Yw, 
						   $Z /  $Zw );

	if ($xr > $epsilon) { $fx = pow($xr, 1/3); } else { $fx = ($kappa*$xr + 16)/116; }
	if ($yr > $epsilon) { $fy = pow($yr, 1/3); } else { $fy = ($kappa*$yr + 16)/116; }
	if ($zr > $epsilon) { $fz = pow($zr, 1/3); } else { $fz = ($kappa*$zr + 16)/116; }

	$L = 116 * $fy - 16;
	$a = 500 * ($fx - $fy);
	$b = 200 * ($fy - $fz);

	return [ $L, $a, $b ];
}


sub Lab_to_XYZ
{
	my ($lab, $xyz_white) = @_;
	my ($L, $a, $b)    = @{$lab};
	my ($Xw, $Yw, $Zw) = @{$xyz_white};
	my ($X, $Y, $Z);

	my $epsilon =  0.008856;
	my $kappa   = 903.3;

	my ($fx, $fy, $fz);
	my ($xr, $yr, $zr);

	if ($L > $kappa*$epsilon) { $yr = pow( ($L + 16)/116, 3 ); } else { $yr = $L / $kappa; }
	if ( $yr > $epsilon )     { $fy =  ($L + 16)/116; } else { $fy  =  ($kappa*$yr + 16)/116; }

	$fx = ($a / 500) + $fy;
	$fz = $fy - ($b / 200);

	if (pow($fx, 3) > $epsilon) { $xr = pow($fx, 3); } else { $xr = (116 * $fx - 16)/$kappa; }
	if (pow($fz, 3) > $epsilon) { $zr = pow($fz, 3); } else { $zr = (116 * $fz - 16)/$kappa; }
	if ($L > $kappa*$epsilon) { $yr = pow(($L + 16)/116, 3);  } else { $yr = $L/$kappa; }

	$X = $xr * $Xw;
	$Y = $yr * $Yw;
	$Z = $zr * $Zw;
				
	return [ $X, $Y, $Z ];
}

sub _mult_v3_m33
{
	my ($v, $m) = @_;
	my $vout = [
				 ( $v->[0] * $m->[0]->[0] + $v->[1] * $m->[1]->[0] + $v->[2] * $m->[2]->[0] ), 
				 ( $v->[0] * $m->[0]->[1] + $v->[1] * $m->[1]->[1] + $v->[2] * $m->[2]->[1] ), 
				 ( $v->[0] * $m->[0]->[2] + $v->[1] * $m->[1]->[2] + $v->[2] * $m->[2]->[2] )
				 ];
	return $vout;
}

sub _apow {
	my ($v, $p) = @_;
	return ($v >= 0 ?
					pow($v, $p) : 
					-pow(-$v, $p));
}


1;
