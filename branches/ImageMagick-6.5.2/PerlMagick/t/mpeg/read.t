#!/usr/bin/perl
#
# Test reading MPEG files
#
# Whenever a new test is added/removed, be sure to update the
# 1..n ouput.
#
BEGIN { $| = 1; $test=1; print "1..2\n"; }
END {print "not ok $test\n" unless $loaded;}
use Image::Magick;
$loaded=1;

require 't/subroutines.pl';

chdir 't/mpeg' || die 'Cd failed';

#
# Motion Picture Experts Group file interchange format (version 2)
#
testRead( 'input.m2v',
  'b3c01688247e4722e5e978668fa374c74ee5e97fb05af12d620791434b380814' );

#
# Motion Picture Experts Group file interchange format
#
++$test;
testRead( 'input.mpg',
  '9a8b39aa0896e5f2c14bffd2fc71c6179dd17fab87c33a9fafd6b4061a793e2c' );

1;
