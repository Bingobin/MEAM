#!/usr/bin/perl

use warnings;
use strict;

my $tf_list = shift or die $!;
#my $cre_sort = "/lustre/home/acct-medkkw/medlyb/project/00.process/combind_somatic/03.MEMOS_24s/APL_CRE.sort.bed";
my $mut_cre = "/lustre/home/acct-medkkw/medlyb/project/00.process/combind_somatic/03.MEMOS_24s/APL_CRE.mutation.sort.bed";
my $wt_cre = "/lustre/home/acct-medkkw/medlyb/project/00.process/combind_somatic/03.MEMOS_24s/APL_CRE.nonmut.sort.bed";
chomp(my $mut_cre_num = `wc -l $mut_cre | awk '{print \$1}'`);
chomp(my $wt_cre_num = `wc -l $wt_cre | awk '{print \$1}'`);

print "TF_name\tType\tTF_CRE_Num\tOR_CRE_Num\tPercent\n";
open IN, "$tf_list" or die $!;
<IN>;
while(<IN>){
	chomp;
	my @tmp = split /\t/;
	my $tf_mut_cre = `bedtools intersect -a $mut_cre -b $tmp[1] -wa | sort | uniq | wc -l | awk '{printf \$1}'`;
	my $tf_wt_cre = `bedtools intersect -a $wt_cre -b $tmp[1] -wa | sort | uniq | wc -l | awk '{printf \$1}'`;
	my $or_mut_cre = $mut_cre_num - $tf_mut_cre;
	my $or_wt_cre = $wt_cre_num - $tf_wt_cre;
	my $mut_pt = sprintf("%.4f",$tf_mut_cre / $mut_cre_num);
	my $wt_pt = sprintf("%.4f",$tf_wt_cre / $wt_cre_num);
	print "$tmp[0]\tmutCRE\t$tf_mut_cre\t$or_mut_cre\t$mut_pt\n";
	print "$tmp[0]\twtCRE\t$tf_wt_cre\t$or_wt_cre\t$wt_pt\n";
}
close IN;
