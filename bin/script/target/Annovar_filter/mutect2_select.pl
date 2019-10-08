#!usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"annovar_txt:s","result:s","SNP:s","indel:s");
my $usage=<<"USAGE";
	program: $0
	-annovar_txt	Input annovar output .txt file
	-result			Output  total result
	-SNP 			Output  SNP only
	-indel			Output 	indel Only
USAGE
die $usage unless($opts{annovar_txt} && $opts{result} && $opts{SNP} && $opts{indel});
open(IN_annovar,$opts{annovar_txt}) or die "Can't find annovar result";
open(OUT,">$opts{result}") or die "Can't write to result file";
open(OUT_SNP,">$opts{SNP}") or die "Can't write SNP to result file";
open(OUT_indel,">$opts{indel}") or die "Can't write SNP to result file";
my $annovar_title=<IN_annovar>;
while(<IN_annovar>){
	chomp;
	my @temp_annovar=split("\t",$_);
	if($temp_annovar[150] eq "PASS"){
		print OUT "$temp_annovar[0]\t$temp_annovar[1]\t$temp_annovar[2]\t$temp_annovar[3]\t$temp_annovar[4]\t$temp_annovar[5]\t$temp_annovar[6]\t$temp_annovar[7]\t$temp_annovar[8]\t$temp_annovar[9]\t$temp_annovar[30]\n";
		if($temp_annovar[3] eq "-" || $temp_annovar[4] eq "-"){
			print OUT_indel "$temp_annovar[0]\t$temp_annovar[1]\t$temp_annovar[2]\t$temp_annovar[3]\t$temp_annovar[4]\t$temp_annovar[5]\t$temp_annovar[6]\t$temp_annovar[7]\t$temp_annovar[8]\t$temp_annovar[9]\t$temp_annovar[30]\n";
		}
		if($temp_annovar[3] ne "-" && $temp_annovar[4] ne "-"){
			print OUT_SNP "$temp_annovar[0]\t$temp_annovar[1]\t$temp_annovar[2]\t$temp_annovar[3]\t$temp_annovar[4]\t$temp_annovar[5]\t$temp_annovar[6]\t$temp_annovar[7]\t$temp_annovar[8]\t$temp_annovar[9]\t$temp_annovar[30]\n";
		}
	}	
}
close IN_annovar;
close OUT;
close OUT_SNP;
close OUT_indel;




