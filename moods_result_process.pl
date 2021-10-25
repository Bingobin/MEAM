#!/usr/bin/perl

use warnings;
use strict;

my $moods = shift or die $!;
my $tf = shift or die $!;
my $merge = shift or die $!;

open IN, "$moods" or die $!;
open OUT, ">./tmp_moods.$tf.bed" or die $!;
while(<IN>){
	chomp;
	my @tmp = split /[:|\-|,|(|)]/;
	my $start = $tmp[1] + $tmp[6];
	my $end = $start + length($tmp[-1]) - 1;
	print OUT "$tmp[0]\t$start\t$end\n";
}
close IN;
close OUT;

system("bedtools sort  -i tmp_moods.$tf.bed  > tmp_moods.$tf.sort.bed") == 0 or die $!;
system("bedtools merge -d -$merge -i tmp_moods.$tf.sort.bed -c 1,2,3 -o count,collapse,collapse  > $tf.position.bed") == 0 or die $!;
`rm tmp_moods.$tf.bed tmp_moods.$tf.sort.bed`;

#if($merge eq "TRUE"){
#	system("bedtools sort  -i tmp_moods.$tf.bed  > tmp_moods.$tf.sort.bed") == 0 or die $!;
#	system("bedtools merge -d 1 -i tmp_moods.$tf.sort.bed -c 1,2,3 -o count,collapse,collapse  > $tf.position.bed") == 0 or die $!;
#	`rm tmp_moods.$tf.bed tmp_moods.$tf.sort.bed`;
#}elsif($merge eq "FALSE"){
#	system("bedtools sort  -i tmp_moods.$tf.bed  > $tf.position.bed") == 0 or die $!;
#	`rm tmp_moods.$tf.bed`;
#}else{
#	die "Please specify whether the bed file needs to be merged: $!";
#}
