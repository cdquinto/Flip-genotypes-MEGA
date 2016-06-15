#!/usr/bin/perl
use warnings;
use File::Slurp;

@array=read_file("test.tped");


for($j=0;$j<127680;$j++){
	chomp $array[$j];
	@split=split(/ /,$array[$j]);
	$string1=$split[4].$split[5];
	$string2=$split[6].$split[7];

	$matches = ($string1 ^ $string2) =~ tr/\0//;
	if ($matches != 2 && $matches !=1){
		print "$array[$j]\n";
	}
}

