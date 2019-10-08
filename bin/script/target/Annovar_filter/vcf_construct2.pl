#!usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"txt:s","vcf:s");
my $usage=<<"USAGE";
	program: $0
	-txt	Input  .txt file
	-vcf			Output result
USAGE
die $usage unless($opts{txt} && $opts{vcf});


open(IN,$opts{txt}) or die "Can't find txt result";
open(OUT,">$opts{vcf}") or die "Can't write to vcf file";
print OUT "## my VCF\n";
print OUT "##INFO=<ID=Func.refGene,Number=.,Type=String,Description=Func.refGene annotation provided by ANNOVAR>\n";
print OUT "##INFO=<ID=Gene.refGene,Number=.,Type=String,Description=Gene.refGene annotation provided by ANNOVAR>\n";
print OUT "##INFO=<ID=GeneDetail.refGene,Number=.,Type=String,Description=GeneDetail.refGene annotation provided by ANNOVAR>\n";
print OUT "##INFO=<ID=ExonicFunc.refGene,Number=.,Type=String,Description=ExonicFunc.refGene annotation provided by ANNOVAR>\n";
print OUT "##INFO=<ID=AAChange.refGene,Number=.,Type=String,Description=AAChange.refGene annotation provided by ANNOVAR>\n";
print OUT "##INFO=<ID=cytoBand,Number=.,Type=String,Description=cytoBand annotation provided by ANNOVAR>\n";
print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
while(<IN>){
	chomp;
	my @temp=split("\t",$_);
	print OUT "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[10]\tFunc.refGene=$temp[5];Gene.refGene=$temp[6];GeneDetail.refGene=$temp[7];ExonicFunc.refGene=$temp[8];AAChange.refGene=$temp[9]\n";
}
close IN;
close OUT;



