#!/usr/bin/perl

## use strict;

#####################################################################################################
## Converts the post processed circular RNA data from STAR + filterCirc.awk to bed format:
##
## - chromosome
## - start (first nucleotide of circle, 0 based)
## - end (last nucleotide of circle, 0 based)
## - name (unique nr + "s" if the circle matches a splice junction)
## - score (nr of supporting reads)
## - strand

## ex: prog/starCirclesToBed.pl MHcircCollapsed.txt


#############
## Arguments
my $starPostFile = $ARGV[0];


################
## 'Main'

open IN, $starPostFile or die "Cannot open $starPostFile\n";
my $i = 0;

my %reads;

while(<IN>){
	chomp;
	my @fields = split /\s+/;
	my $Identifier = $fields[0].$fields[1].$fields[2].$fields[3];
	$reads{$Identifier} .= $fields[7].',';	
}


seek IN, 0, 0;

while(<IN>){
  chomp;
  my @fields = split /\s+/;

  my $Identifier = $fields[0].$fields[1].$fields[2].$fields[3];
  my $score = scalar(split(/,/,$reads{$Identifier}));
  my $chr = $fields[0];
  my $start = $fields[1];
  my $end = $fields[2] - 1;
  my $strand = $fields[3];

  print $chr."\t".$start."\t".$end."\t".$score."\t".$strand."\t".$reads{$Identifier}."\n";
}
