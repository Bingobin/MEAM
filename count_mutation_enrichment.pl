#!/usr/bin/perl

use warnings;
use strict;

my $stat = shift or die $!;
my $motif_len = shift or die $!;
my $bg_mut = shift or die $!;

my @bg;
open IN, "$bg_mut" or die $!;
while(<IN>){
	chomp;
	my @tmp = split /\t/;
	push @bg, $tmp[0];
}
close IN;
my $bg_mean = mean(\@bg);
my $bg_sd = standard_deviation(\@bg, $bg_mean);
#my $bg_sd = 1;

#print "$bg_mean\t$bg_sd\n";

sub mean{
	my $arr = shift;
	my $sum = 0;
	my $num = scalar(@$arr);
	for my $i (@$arr){
		$sum += $i;
	}
	return sprintf("%.2f", $sum/$num);
}

sub standard_deviation{
	my $arr = shift;
	my $mean = shift;
	my $num = scalar(@$arr);
	my $var = 0;
	for my $i (@$arr){
		$var += ($i - $mean) ** 2
	}
	$var = $var / $num;
	return sqrt($var);
}


my $flank = 500;
my $bin = 10;
#my $motif_len = 8;

my %hash;
open IN, "$stat" or die $!;
<IN>;
while(<IN>){
	chomp;
	my @tmp = split /\t/;
	if(defined $hash{$tmp[0]}){
		$hash{$tmp[0]} += $tmp[1];
	}else{
		$hash{$tmp[0]} = 0;
		$hash{$tmp[0]} += $tmp[1];
	}
}
close IN;

#for my $i ( sort {$a<=>$b} keys %hash ){
#	print "$i\t$hash{$i}\n";
#}

my %result;

print "POS\tNUM\tZS\n";
my $val_0 = $hash{0}/$motif_len*10;
$result{0} = $val_0;
#print "0\t$val_0\t". ($val_0 - $bg_mean)/$bg_sd ."\n";
#for my $i (1..$flank/$bin){
#	my $key = $i*$bin;
#	my $val = ($hash{$key} + $hash{-$key}) / 2;
#	print "$key\t$val\n";
#}
my $val_accumulate = $val_0;
for my $i (1..5){
	my $key = $i*$bin;
	my $val = ($hash{$key} + $hash{-$key}) / 2;
#	print "$key\t$val\n";
	$val_accumulate += $val;
	$result{$i*$bin} = $val_accumulate/($i+1);
#	print "$key\t" . $val_accumulate/($i+1) . "\n";
}

my $val_100 = 0;
for my $i (6..10){
	my $key = $i*$bin;
	my $val = ($hash{$key} + $hash{-$key}) / 2;
	$val_100 +=  $val;
}
$val_accumulate += $val_100;
$result{100} = $val_accumulate/(100/$bin+1);
#print "100\t$val_acc_100\t" . $val_accumulate/(100/$bin+1) . "\n";

my ($val_200, $val_300, $val_400, $val_500) = (0, 0, 0, 0);

for my $i(11..$flank/$bin){
	my $key = $i*$bin;
	my $val = ($hash{$key} + $hash{-$key}) / 2;
	if($i <= 20 && $i > 10 ){
		$val_200 += $val;		
	}elsif( $i <= 30 && $i > 20 ){
		$val_300 += $val;
	}elsif( $i <= 40 && $i > 30 ){
		$val_400 += $val; 
	}elsif( $i <= 50 && $i > 40 ){
		$val_500 += $val; 
	}
}
$val_accumulate += $val_200; $result{200} = $val_accumulate/(200/$bin+1);
#print "200\t" . $val_accumulate/(200/$bin+1) . "\n";
$val_accumulate += $val_300; $result{300} = $val_accumulate/(300/$bin+1);
#print "300\t" . $val_accumulate/(300/$bin+1) . "\n";
$val_accumulate += $val_400; $result{400} = $val_accumulate/(400/$bin+1);
#print "400\t" . $val_accumulate/(400/$bin+1) . "\n";
$val_accumulate += $val_500; $result{500} = $val_accumulate/(500/$bin+1);
#print "500\t" . $val_accumulate/(500/$bin+1) . "\n";
my @pos = qw( 0 10 20 30 40 50 100 200 300 400 500 );
for my $i (@pos){
	print "$i\t$result{$i}\t" . ($result{$i} - $bg_mean)/$bg_sd . "\n";
}
