#!/usr/bin/perl -w
use warnings; 

%hash=();

open(OUT,">list_snps_to_change_alleles.txt");

open(FILE,"maavp1v1.raw_FinalReport.txt");
while(<FILE>){
	chomp; 
	@split=split(/\t/,$_);
	if (!exists ($hash{$split[0]})){
		if($split[2] ne $split[6] && $split[3] ne $split[7] && $split[2] ne $split[3]){ 
			$hash{$split[0]}=["$split[2]\t","$split[3]\t","$split[6]\t","$split[7]\t","$split[9]\t"];
			}
	}
	else{ 
		next;
	}
}
close(FILE);

for $key (keys %hash){
	print OUT "$key\t@{$hash{$key}}\n";
}
close(OUT);
