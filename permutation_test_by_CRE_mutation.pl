#!/usr/bin/perl

use warnings;
use strict;

my $cre = shift or die $!;
my $maf = shift or die $!;
my $region = shift or die $!;
my $tf = shift or die $!;

my $iter = 5000;

my @cre_pos;
open IN, "$cre" or die $!;
while(<IN>){
	chomp;
	my @tmp = split /\t/;
	for my $i ($tmp[1]..$tmp[2]){
		push @cre_pos, "$tmp[0]\t$i";
	}	
}
close IN;

my $all = scalar(@cre_pos);
my @mut_num;

my @perm_num;
my $start_num;
if( -e "$tf.bg.txt"){
	open IN, "$tf.bg.txt" or die $!;
	while(<IN>){
		chomp;
		my @tmp = split /\t/;
		push @perm_num, $tmp[0] if ($tmp[2] == $region);		
	}
	close IN;
	`mv $tf.bg.txt $tf.bg.txt.bak`;
	open OUT, ">$tf.bg.txt" or die $!;
	for (my $m=0; $m < @perm_num; $m++){
		print OUT "$perm_num[$m]\t" . ($m+1) . "\t$region\n";
	}
	close OUT;
	$start_num = scalar(@perm_num) + 1;
}else{
	$start_num = 1;
}
open OUT, ">>$tf.bg.txt" or die $!;
###########################
for (my $k=$start_num; $k<=$iter; $k++){
	my %hash;
	my $num;
	for (my $n=1; $n<=$region; $n++){
		do{
			$num = int(rand($all));
		}while(defined $hash{$num});
		$hash{$num} = 1;
	}

	my %position;
	for my $i (sort {$a<=>$b} keys %hash){
		my @tmp = split("\t", $cre_pos[$i-1]);
		$position{$tmp[0]}{$tmp[1]} = $tmp[1] + 99;
#		$position{$tmp[0]}{$tmp[1]} = $tmp[1] + 9;
	}
	open OUT1, ">./tmp_rand_bg.$tf.bed" or die $!;
	for my $i (sort{$a cmp $b} %position){
		for my $j (sort {$a<=>$b} keys %{$position{$i}}){
			print OUT1 "$i\t$j\t$position{$i}{$j}\n";
		}
	}
	close OUT1;
	open IN, "$maf" or die $!;
	open OUT2, ">./tmp_mutation.$tf.bed" or die $!;
	<IN>;
	while(<IN>){
		chomp;
		my @tmp = split /\t/;
		print OUT2 "$tmp[4]\t$tmp[5]\t$tmp[6]\n";
	}
	close IN;
	close OUT2;
	`bedtools intersect -a tmp_rand_bg.$tf.bed  -b tmp_mutation.$tf.bed -wa -wb > tmp_intersect.$tf.bed`;
	chomp(my $mutation_num = `cut -f 4-6 tmp_intersect.$tf.bed | sort | uniq | wc -l | awk '{print \$1}'`);
	$mutation_num =  $mutation_num / 10;
#	chomp(my $mutation_num = `wc -l tmp_intersect.bed | awk '{print \$1}'`);
	print OUT "$mutation_num\t$k\t$region\n";
#	print STDERR "Iteration $k:Mutation num = $mutation_num;Region num = $region\n";
#	push @mut_num, $mutation_num;
}
close OUT;
`rm tmp_rand_bg.$tf.bed tmp_mutation.$tf.bed tmp_intersect.$tf.bed`;
#print "@mut_num\n";
