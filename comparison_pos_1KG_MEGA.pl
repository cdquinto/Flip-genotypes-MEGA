#!/usr/bin/perl 

use warnings;

%hash=();

open(FILE,"../list_pos_1KG.txt");
while(<FILE>){
	chomp;
	@split=split(/\t/,$_);
	$hash{$split[1]}=$split[0];
}
close(FILE);

open(FILE2,"list_pos_MEGA.txt");
while(<FILE2>){
	chomp;
	@split2=split(/\t/,$_);
	foreach $key (%hash){
		if ($key eq $split2[1]){
			print "$hash{$key}\t$split2[0]\t$key\n";
		}
	}
}
close(FILE2);
