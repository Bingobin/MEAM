#!/usr/bin/perl

use warnings;
use strict;
use Cwd;
use Getopt::Long;
use FindBin qw($Bin);

=head1 MEAM_pipeline.pl
Usage: 
		perl MEAM_pipeline.pl  -t MYB -l 10  <APL_CRE.mutation.bed>  <WL1109_APL_D0_wgs_mutation.mutect2.maf>
options:
		-tf_name|t    <str>    transcript name
		-motif_dir|d  <str>    the pfm database dir of motif
		-mood_pv|p    <flo>    the pvalue of MOOD scan the motif position
		-merge|m      <int>    Maximum distance between features allowed for features to be merged, default 10
		-step|s       <int>    Step,default 0, total 6 steps.
#		-out|o        <str>    the output dirtory
		-help
=cut

GetOptions(
		'tf_name|t=s' => \my $tf_name,
		'motif_dir|d=s' => \my $motif_dir,
		'mood_pv|p=f' => \my $mood_pv,
		'merge|m=i' => \my $merge,
		'step|s=i' => \my $step,
#		'out|o=s' => \my $out,
		'help|h' => \my $Help		
);


die `pod2text $0` if($Help);

die "Must specify <interest profile of genome region bed file> and <the WGS mutation MAF file>, use -h to view help!\n" if @ARGV < 2;

$motif_dir ||= "/lustre/home/acct-medkkw/medlyb/wl_proj/WL234_Lib/database/JASPAR/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar_MOOD";
$mood_pv ||= 0.0001;
$merge ||= 10;
$step ||= 0;
#$out ||= "./OUT_MEAM";
#`mkdir $out` unless (-e $out);

my $region_file = shift;
my $mutation_maf = shift;

##############################Step01
my $date = localtime(); 
print STDERR "[$date]:<Step01> extrcat the sequence of interest region.\n";
my $cmd1 = "perl $Bin/region_extract_genome.pl $region_file  0 > tmp_interest_region.$tf_name.fa";
print STDERR "[$date]:<CMD> $cmd1\n";
if($step <=1){
	system($cmd1) == 0 or die $!;
}else{
	print STDERR "[$date]:<Step01> skiped!\n";
}
#############################Step02
$date = localtime(); 
print STDERR "[$date]:<Step02> find the motif position of spcified TF by MOODs.\n"; 
my $cmd2 = "moods-dna.py  -m $motif_dir/*_$tf_name*pfm -s tmp_interest_region.$tf_name.fa  -p $mood_pv  > $tf_name.position.txt";
print STDERR "[$date]:<CMD> $cmd2\n";
if($step <= 2){
	system($cmd2) == 0 or die $!;
}else{
	print STDERR "[$date]:<Step02> skiped!\n";
}
#############################Step03
$date = localtime();
print STDERR "[$date]:<Step03> transform MOODS result to BED formats.\n";
my $cmd3 = "perl $Bin/moods_result_process.pl $tf_name.position.txt $tf_name $merge";
print STDERR "[$date]:<CMD> $cmd3\n";
if($step <= 3){
	system($cmd3) == 0 or die $!;
}else{
	print STDERR "[$date]:<Step03> skiped!\n";
}
my $size = `wc -l $tf_name.position.bed | awk '{printf \$1}'`;
my $motif_len = `awk '{i+=\$3-\$2+1; j+=1}END{k=i/j;printf k}' $tf_name.position.bed`;
print STDERR "[$date]:Total num of $tf_name mitof binding: $size\n";
print STDERR "[$date]:The average of $tf_name mitof bingding length: $motif_len\n";
#############################Step04
$date = localtime();
print STDERR "[$date]:<Step04> calculated the mutation num at motif flank region.\n";
my $cmd4 = "perl $Bin/cal_mut_num_by_motif.pl $tf_name.position.bed $motif_len $mutation_maf 1>$tf_name.mutation.stat 2>$tf_name.mutation.info ";
print STDERR "[$date]:<Cmd> $cmd4\n";
if($step <=4){
	system($cmd4) == 0 or die $!;
}else{
	print STDERR "[$date]:<Step04> skiped!\n";	
}
#############################Step05
$date = localtime(); 
print STDERR "[$date]:<Step05> compute the background mutation distribution by permutation test.\n";
my $cmd5 = "perl $Bin/permutation_test_by_CRE_mutation_v2.pl $region_file  $mutation_maf  $size $tf_name";
print STDERR "[$date]:<CMD> $cmd4\n";
if($step <=5){
	system($cmd5) == 0 or die $!;
}else{
	print STDERR "[$date]:<Step05> skiped!\n";
}
############################Step06
$date = localtime(); 
print STDERR "[$date]:<Step06> calculated the cumulative distribution of mutation in motif flank region(10,20,30,40,50,100,200,300,400,500).\n";
my $cmd6 = "perl $Bin/count_mutation_enrichment.pl $tf_name.mutation.stat $motif_len $tf_name.bg.txt > $tf_name.mutation.count";
print STDERR "[$date]:<Cmd> $cmd6\n";
if($step <=6){
	system($cmd6) == 0 or die $!;
	`rm tmp_interest_region.$tf_name.fa`;
}else{
	print STDERR "[$date]:<Step06> skiped!\n";
}
############################
$date = localtime(); print STDERR "[$date]:All Done!\n";
