#!usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"annovar_txt:s","result:s");
my $usage=<<"USAGE";
	program: $0
	-annovar_txt	Input annovar output .txt file
	-result			Output result
USAGE
die $usage unless($opts{annovar_txt} && $opts{result});


open(IN_annovar,$opts{annovar_txt}) or die "Can't find annovar result";
open(OUT,">$opts{result}") or die "Can't write to result file";
my $annovar_title=<IN_annovar>;
while(<IN_annovar>){
	chomp;
	my @temp_annovar=split("\t",$_);
	my @somatic=split(";",$temp_annovar[151]);
	if($somatic[1] eq "SOMATIC"){
		my @somatic_temp=split("=",$somatic[5]);
		if($somatic_temp[1] < 0.01){
			print OUT "$temp_annovar[0]\t$temp_annovar[1]\t$temp_annovar[2]\t$temp_annovar[3]\t$temp_annovar[4]\t$temp_annovar[5]\t$temp_annovar[6]\t$temp_annovar[7]\t$temp_annovar[8]\t$temp_annovar[9]\t$temp_annovar[30]\n";
		}
	}
}
close IN_annovar;
close OUT;



