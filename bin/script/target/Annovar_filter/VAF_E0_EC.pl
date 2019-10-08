#!usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"SNP:s","bamreadcount_EC:s","bamreadcount_E0:s","output_file:s");
my $usage=<<"USAGE";
	program: $0
	-SNP	Input SNP file
	-bamreadcount_EC
	-bamreadcount_E0
	-output_file	Output smmary file
USAGE
die $usage unless($opts{SNP} && $opts{bamreadcount_EC} && $opts{bamreadcount_E0}&& $opts{output_file});

open(IN_SNP,$opts{SNP}) or die "Can't find SNP file";
open(IN_bamreadcount_EC,$opts{bamreadcount_EC}) or die "Can't find EC bamreadcount";
open(IN_bamreadcount_E0,$opts{bamreadcount_E0}) or die "Can't find E0 bamreadcount";
open(OUT,">$opts{output_file}") or die "Can't write result to output file";


#bamreadcount2hash
##E0
my %bamreadcount2hash_E0;
while(<IN_bamreadcount_E0>){
	chomp;
	my @bamreadcount=split("\t",$_);
	my $labe_E0=join(":",$bamreadcount[0],$bamreadcount[1]);
	$bamreadcount2hash_E0{$labe_E0}=$_;
}
close IN_bamreadcount_E0;

##EC
my %bamreadcount2hash_EC;
while(<IN_bamreadcount_EC>){
	chomp;
	my @bamreadcount=split("\t",$_);
	my $labe_EC=join(":",$bamreadcount[0],$bamreadcount[1]);
	$bamreadcount2hash_EC{$labe_EC}=$_;
}
close IN_bamreadcount_EC;

print OUT "CHROM\tSTART\tEND\tREF\tALT\tFunc_refGene\tGene_refGene\tgenomicSuperDups\tExonicFunc_refGene\tGeneDetail_refGene\trs_number\tDepth_E0\tDepth_alt_E0\tVAF_E0\tDepth_EC\tDepth_alt_EC\tVAF_EC\n";

while(<IN_SNP>){
	chomp;
	my @temp_snp=split("\t",$_);
	my $labe_final=join(":",$temp_snp[0],$temp_snp[1]);
	my $depth_E0_final=0;;my $depth_EC_final=0;
	my $alt_depth_E0_final=0;my $alt_depth_EC_final=0;
	my $VAF_E0_final=0;my $VAF_EC_final=0;

	if(exists $bamreadcount2hash_E0{$labe_final} ) {		
		my $information_E0= $bamreadcount2hash_E0{$labe_final};
		my @info_E0=split("\t",$information_E0);
		
		if ($#info_E0 >0) { 
			for (my $m=5;$m<=$#info_E0;$m++) {
				my @Indel=split(":",$info_E0[$m]);
				my $sequence=$Indel[0];	
				$depth_E0_final=$info_E0[3];
				#print $temp_snp[4],"\n";
				#my $sequence=substr($Indel[0],1);  ##substr the indel sequence in bamreadcount
				if ($temp_snp[3] eq "-" or $temp_snp[4] eq "-") {
					if ($Indel[0]=~ /\W+/) {
						my $sequence=substr($Indel[0],1);
						if (uc($sequence) eq $temp_snp[3] or uc($sequence) eq $temp_snp[4] ) {
							$alt_depth_E0_final=$Indel[1];
							if ($depth_E0_final >0) {     ##caculate VAF
								$VAF_E0_final=$alt_depth_E0_final/$depth_E0_final;
							}
							else {
								$VAF_E0_final=0;
							}
						}
					}
				}
				elsif (uc($sequence) eq $temp_snp[4] ) {
					$alt_depth_E0_final=$Indel[1];
					if ($depth_E0_final >0) {     ##caculate VAF
						$VAF_E0_final=$alt_depth_E0_final/$depth_E0_final;
					}
					else {
						$VAF_E0_final=0;
					}
				}
			}
		}

	}
	if(exists $bamreadcount2hash_EC{$labe_final} ) {		
		my $information_EC= $bamreadcount2hash_EC{$labe_final};
		my @info_EC=split("\t",$information_EC);
		
		if ($#info_EC >0) { 
			for (my $m=5;$m<=$#info_EC;$m++) {
				my @Indel=split(":",$info_EC[$m]);
				my $sequence=$Indel[0];	
				$depth_EC_final=$info_EC[3];
				#print $temp_snp[4],"\n";
				#my $sequence=substr($Indel[0],1);  ##substr the indel sequence in bamreadcount
				if ($temp_snp[3] eq "-" or $temp_snp[4] eq "-") {
					if ($Indel[0]=~ /\W+/) {
						my $sequence=substr($Indel[0],1);
						if (uc($sequence) eq $temp_snp[3] or lc($sequence) eq $temp_snp[3] or uc($sequence) eq $temp_snp[4] or lc($sequence) eq $temp_snp[4]) {
							$alt_depth_EC_final=$Indel[1];
							if ($depth_EC_final >0) {     ##caculate VAF
								$VAF_EC_final=$alt_depth_EC_final/$depth_EC_final;
							}
							else {
								$VAF_EC_final=0;
							}
						}
					}
				}
				elsif (uc($sequence) eq $temp_snp[4] or lc($sequence) eq $temp_snp[4] ) {
					$alt_depth_EC_final=$Indel[1];
					if ($depth_EC_final >0) {     ##caculate VAF
						$VAF_EC_final=$alt_depth_EC_final/$depth_EC_final;
					}
					else {
						$VAF_EC_final=0;
					}
				}
			}
		}

	}
	print OUT "$_\t$depth_E0_final\t$alt_depth_E0_final\t$VAF_E0_final\t$depth_EC_final\t$alt_depth_EC_final\t$VAF_EC_final\n";	
}	
	
close IN_SNP;
close OUT;



		
