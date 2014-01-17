# $Id$

package Sanger::Graphics::Glyph::Triangle;

### Usage:
###   new Sanger::Graphics::Glyph::Triangle({
###     width     => 5,       # Length of baseline
###     height    => 10,      # Distance from point to baseline
###     direction => 'up',    # Can be up, down, left, right - refers to the direction of the triangle
###     mid_point => [ 0, 1 ] # x, y coordinates of the point of the triangle
###   });

use strict;

use Sanger::Graphics::Glyph::Rect;

use base qw(Sanger::Graphics::Glyph::Poly);

sub new {
  my $class     = shift;
  my $args      = shift;
  my $width     = $args->{'width'} / 2;
  my $height    = $args->{'height'};
  my $direction = $args->{'direction'};
  my ($x, $y)   = @{$args->{'mid_point'}};
  my $rectangle;
  
  my %dir = (
    up    => sub { return ([ $x - $width,  $y + $height, $x + $width,  $y + $height ], { x => $x - $width,  y => $y           }); },
    down  => sub { return ([ $x - $width,  $y - $height, $x + $width,  $y - $height ], { x => $x - $width,  y => $y - $height }); },
    left  => sub { return ([ $x + $height, $y - $width,  $x + $height, $y + $width  ], { x => $x,           y => $y - $width  }); },
    right => sub { return ([ $x - $height, $y - $width,  $x - $height, $y + $width  ], { x => $x - $height, y => $y - $width  }); }
  );
  
  my ($points, $rect) = exists $dir{$direction} ? &{$dir{$direction}} : ();
  
  if (!$args->{'no_rectangle'}) {
    # Height refers to the distance between the point of the triangle and the baseline, which is actually the width of the rectangle if the triangle points left or right
    ($args->{'height'}, $args->{'width'}) = ($args->{'width'}, $args->{'height'}) if $direction =~ /^(left|right)$/;
    
    # Make an invisible rectangle to go on top of the triangle to make clicking for z-menus easier
    $rectangle = new Sanger::Graphics::Glyph::Rect({ %$args, %$rect, bordercolour => undef, colour => undef });
  }
  
  $args->{'points'} = [ $x, $y, @$points ];
  
  delete $args->{$_} for qw(width height mid-point direction no_rectangle);
  
  my $triangle = $class->SUPER::new($args);
  
  return $rectangle ? ($triangle, $rectangle) : $triangle;
}

1;