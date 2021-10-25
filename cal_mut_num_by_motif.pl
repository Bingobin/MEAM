#!/usr/bin/perl

use warnings;
use strict;

my $motif_pos = shift or die $!;
my $motif_len = shift or die $!;
my $maf = shift or die $!;

my $flank = 500;
my $bin = 10;

my %hash;
my $motif_num = 0;
open IN, "$motif_pos" or die $!;
while(<IN>){
	chomp;
#	my @tmp = split /[:|\-|,|(|)]/;
#	print join("\t", @tmp) . "\n";
	my @tmp = split /\t/;
#	my $start = $tmp[1] + $tmp[6];
#	my $end = $start + $motif_len - 1;
#	print "$tmp[0]:$start-$end\n";
	my $start = $tmp[1];
	my $end = $tmp[2];
	$hash{$tmp[0]}{$start} = $end;
	$motif_num += 1;
}
close IN;
#my $weight = 5000 / $motif_num;
my $weight = 1;

my %mutation;
my @tvti = qw( C>A C>G C>T T>A T>C T>G INS DEL);
$mutation{"0"}{"num"} = 0;
for my $k (@tvti){ 
	$mutation{"0"}{$k} =0;
}
for my $n (1 .. $flank/$bin){
	$mutation{-$n*$bin}{"num"} = 0;
	$mutation{$n*$bin}{"num"} = 0;
	for my $k (@tvti){
		$mutation{-$n*$bin}{$k} = 0;
		$mutation{$n*$bin}{$k} = 0;
	}
}


open IN, "$maf" or die $!;
<IN>;
while(<IN>){
	chomp;
	my @tmp = split /\t/;
	my $flag = 0;
	for my $i (sort {$a<=>$b} keys %{$hash{$tmp[4]}}){
		my $start = $i;
		my $end = $hash{$tmp[4]}{$i};
#		if($flag == 0){
			if($tmp[5] >= $start && $tmp[5] <= $end){
				$mutation{"0"}{"maf"} = $_;
				$mutation{"0"}{"num"} += 1;
				cal_tvti_indel($tmp[9], $tmp[10], $tmp[12], \%mutation, 0);
				$flag +=1;
				print STDERR "$tmp[4]\t$start\t$end\t$_\n";
			}
#		}else{
#			next;
#		}
		for my $n (1 .. $flank/$bin){
#			if( ( $tmp[5] >= $start - $n*$bin && $tmp[5] < $start - ($n-1)*$bin ) || ($tmp[5] <= $end + $n*$bin && $tmp[5] > $end + ($n-1)*$bin) ){
#				$mutation{$n*$bin}{$_} = 1;
#			}
#			if($flag == 0 ){
				if( $tmp[5] >= $start - $n*$bin && $tmp[5] < $start - ($n-1)*$bin ){
					$mutation{-$n*$bin}{"maf"} = $_;
					$mutation{-$n*$bin}{"num"} += 1;
					cal_tvti_indel($tmp[9], $tmp[10], $tmp[12], \%mutation, -$n*$bin);
					$flag += 1;
				}
#			}else{
#				next;
#			}
#			if($flag == 0 ){
				if( $tmp[5] <= $end + $n*$bin && $tmp[5] > $end + ($n-1)*$bin ){
					$mutation{$n*$bin}{"maf"} = $_;
					$mutation{$n*$bin}{"num"} += 1;
					cal_tvti_indel($tmp[9], $tmp[10], $tmp[12], \%mutation, $n*$bin);
					$flag += 1;
				}
#			}else{
#				next;
#			}
		}
	}
#	print STDERR "$flag\n";
}
close IN;

my $format = 1;
if($format == 0){
	print "POS\t" . join("\t", @tvti) . "\tNUM\n";
	for my $i (sort {$a <=> $b} keys %mutation){
		print "$i";
		for my $k (@tvti){
			print "\t$mutation{$i}{$k}*$weight";
		}
		print "\t$mutation{$i}{'num'}*$weight\n";
	}
}else{
	print "POS\tNUM\tTVTI\n";
	for my $i (sort {$a <=> $b} keys %mutation){
		for my $k (@tvti){
			print "$i\t" . $mutation{$i}{$k}*$weight . "\t$k\n";
		}
	}
}

sub cal_tvti_indel{
	my $type = shift;
	my $ref = shift;
	my $var = shift;
	my $hash = shift;
	my $key = shift;
	if($type eq "SNP"){
		if( ($ref eq "C" && $var eq "A") || ($ref eq "G" && $var eq "T") ){ $hash->{$key}->{"C>A"} += 1 }
		if( ($ref eq "C" && $var eq "G") || ($ref eq "G" && $var eq "C") ){ $hash->{$key}->{"C>G"} += 1 }
		if( ($ref eq "C" && $var eq "T") || ($ref eq "G" && $var eq "A") ){ $hash->{$key}->{"C>T"} += 1 }
		if( ($ref eq "T" && $var eq "A") || ($ref eq "A" && $var eq "T") ){ $hash->{$key}->{"T>A"} += 1 }
		if( ($ref eq "T" && $var eq "C") || ($ref eq "A" && $var eq "G") ){ $hash->{$key}->{"T>C"} += 1 }
		if( ($ref eq "T" && $var eq "G") || ($ref eq "A" && $var eq "C") ){ $hash->{$key}->{"T>G"} += 1 }
	}elsif($type eq "INS"){
		$hash->{$key}->{"INS"} += 1;
	}elsif($type eq "DEL"){
		$hash->{$key}->{"DEL"} += 1;	
	}
}
