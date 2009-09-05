#!/usr/bin/perl
#
# Example: Modify all the pixels in an image (like -fx)
#
# Anthony Thyssen   5 October 2007
#
use strict;
use Image::Magick;

# read original image
my $orig = Image::Magick->new();
my $w = $orig->Read('rose:');
warn("$w")  if $w;
exit  if $w =~ /^Exception/;


# make a clone of the image for modifications
my $dest = $orig->Clone();

# You could enlarge destination image here if you like.

# Iterate over destination image...
my ($width, $height) = $dest->Get('width', 'height');

for( my $j = 0; $j < $height; $j++ ) {
  for( my $i = 0; $i < $width; $i++ ) {

    # read original image color
    my @pixel = $orig->GetPixel( x=>$i, y=>$j );

    # modify the pixel values (as normalized floats)
    $pixel[0] = $pixel[0]/2;      # darken red

    # write pixel to destination
    # (quantization and clipping happens here)
    $dest->SetPixel(x=>$i,y=>$j,color=>\@pixel);
  }
}

# display the result (or you could save it)
$dest->Write('win:');
$dest->Write('pixel_fx.gif');

