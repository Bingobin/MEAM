#!/usr/bin/perl

use warnings;
use strict;

my $tf_list = shift or die $!;
my $cre_sort = "/lustre/home/acct-medkkw/medlyb/project/00.process/combind_somatic/03.MEMOS_24s/APL_CRE.sort.bed";
chomp(my $cre_num = `wc -l $cre_sort | awk '{print \$1}'`);

print "TF_name\tType\tNumber\tPeak_num\tSqrt(jaccard)\tTF_CRE_Num\tTotal_CRE_Num\n";
open IN, "$tf_list" or die $!;
<IN>;
while(<IN>){
	chomp;
	my @tmp = split /\t/;
	chomp(my $peak_num = `wc -l $tmp[1] | awk '{print \$1}'`);
	`bedtools jaccard -a $cre_sort -b $tmp[1] > tmp_jaccard.bed`;
	my $jaccard;
	my $intersect;
	open INI, "./tmp_jaccard.bed" or die $!;
	<INI>;
	while(<INI>){
		chomp;
		my @tmp2 = split /\t/;
		$jaccard  = sprintf("%.4f",sqrt($tmp2[2]));
		$intersect = $tmp2[3];
	}
	close INI;
	`rm tmp_jaccard.bed`;
	my $cre_inter = "/lustre/home/acct-medkkw/medlyb/project/00.process/combind_somatic/03.MEMOS_24s/$tmp[0]_MEMOS/$tmp[0]_APL_CRE.bed";
	chomp(my $incre_num = `wc -l $cre_inter | awk '{print \$1}'`);
	print "$tmp[0]\tinCRE\t$intersect\t$peak_num\t$jaccard\t$incre_num\t$cre_num\n";
	print "$tmp[0]\tnoCRE\t". ($peak_num - $intersect) . "\t$peak_num\t$jaccard\t$incre_num\t$cre_num\n";
}
close IN;
