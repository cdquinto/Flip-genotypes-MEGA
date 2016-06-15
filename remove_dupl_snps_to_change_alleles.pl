#!/usr/bin/perl
use warnings;

open(FILE,"list_snps_to_change_alleles.txt2");
open(OUT,">list_snps_to_change_alleles.txt");

%ids=();
while(<FILE>){
	chomp;
	@split=split(/\t/,$_);
	
	if(!exists($ids{$split[0]})){
		$ids{$split[0]} = 1; 
		print OUT "$_\n";
	}	
}
close(FILE);
close(OUT); 