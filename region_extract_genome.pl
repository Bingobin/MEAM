#!/usr/bin/perl

use strict;
use warnings;

die "perl $0  <chr_bed_file>  <basd_num> :$!" if (@ARGV != 2);

my $file = shift;
my $base = shift;

my $ref = "/lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/bwaindex2/Homo_sapiens_assembly38.fasta";

my @regin;
open IN, "$file" or die $!;
#<IN>;
while(<IN>){
	chomp;
	push @regin, $_;
}

$/ = "\n>";
open IN, "$ref" or die $!;
while(<IN>){
	chomp;
	s/^>//;
	my ($id,$seq) = split(/\n/, $_, 2);
	$id =~ /(\S+?)\s+/;
	$id = $1;
	$seq =~ s/\n//g;
	for my $i (@regin){
		my @tmp = split(/\t/,$i);
		if($id eq $tmp[0]){
			my $region = substr($seq, $tmp[1]-1, ($tmp[2]-$tmp[1]+1));
			if($base > 0){
				my $r_front = substr($seq, $tmp[1]-$base-1, $base);
				my $r_back = substr($seq, $tmp[1], $base);
				print ">$i\n$r_front\t$region\t$r_back\n";
			}else{
#				print ">". join("_", @tmp) .  "\n$region\n";
				print ">$tmp[0]:$tmp[1]-$tmp[2]($tmp[3])\n$region\n";
			}
		}
	}
}
close IN;
