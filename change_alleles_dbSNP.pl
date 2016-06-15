#!/usr/bin/perl
use warnings;

open(FILE,"mart_export.txt");
readline(FILE);
open(OUT,">snps_to_add_from_dbSNP.txt");
open(OUT2,">multiallelic_SNPs.txt");
while(<FILE>){
	chomp;
	@split=split(/,/,$_);
	@alleles=split(/\//,$split[5]);
	if(scalar(@alleles)>2){
		print OUT2 "$split[0]\t@alleles\n";
	}
	###Alleles from dbSNP should appear in the first columns
	else {
		$a1=$alleles[0];
		$a2=$alleles[1];
		###change the alleles
		$alleles[0]=~tr/ATGC/TACG/;
		$alleles[1]=~tr/ATGC/TACG/;
		print OUT "$split[0]\t$alleles[0]\t$alleles[1]\t$a1\t$a2\t$split[3]\n";
	}
}
close(FILE);
close(OUT);
close(OUT2);

	

