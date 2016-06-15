#!/usr/local/bin/perl

use warnings; 

$ped=$ARGV[0];
@name=split(/\.ped/,$ped);

print "\n\n\n***INICIANDO PIPELINE***\n\n\n";

##Order the plink files by chromosome number and transpose the data
system("plink --file $name[0] --allow-no-sex --recode --not-chr 0 --out $name[0]");
system("plink --file $name[0] --allow-no-sex --recode transpose --not-chr 0 --out $name[0]");

print "\n\n\n***BUSCANDO DUPLICADOS***\n\n\n";

###Find duplicates using plink
system("plink --file $name[0] --list-duplicate-vars --out duplicates_plink");
system("sed -i 1d duplicates_plink.dupvar");
system("sed -i 's/ /\t/g' duplicates_plink.dupvar");

print "\n\n\n***IMPRIMIENDO INFO DE DUPLICADOS***\n\n\n";
print "\n\n\n***IGNORAR ERRORES DE PERL***\n\n\n";

##Get summary of the duplicates
%ids=();

open(IDS,">duplicate_snps_${name[0]}.txt");
open(POS,">duplicate_positions_${name[0]}.txt");

open(DUP,"duplicates_plink.dupvar");
while(<DUP>){
	chomp; 
	@split=split(/\t/,$_);
	
	if($split[5] ne ""){
		if($split[3] eq $split[4] && $split[3] eq $split[5]){
			if(!exists($ids{$split[3]})){ 
				$ids{$split[3]} = 3;
			} 	
		}
		elsif($split[3] ne $split[4] && $split[3] ne $split[5] && $split[4] eq $split[5]){
			print POS "$split[0] $split[1] $split[2] $split[3] $split[4]\n"; 
		}
		elsif($split[3] ne $split[4] && $split[3] ne $split[5] && $split[4] ne $split[5]){
			print POS "$split[0] $split[1] $split[2] $split[3] $split[4] \$split[5]\n"; 
		}		
	}
	
	else{
		if($split[3] eq $split[4]){
			$ids{$split[3]} = 2;
		}
	
		if ($split[3] ne $split[4]){
			if ($split[3]=~/^rs/){
				print POS "$split[0] $split[1] $split[2] $split[4] $split[3]\n";
			}
			else{
				print POS "$split[0] $split[1] $split[2] $split[3] $split[4]\n";
			}
		}	
	}
}

foreach $key ( sort keys %ids ){
	if( $ids{$key} != 1){
		print IDS "$key\t$ids{$key}\n"; 
	}     
}
close(IDS);
close(POS);
close(DUP);

print "\n\n\n***ELIMINANDO SNPS DUPLICADOS***\n\n\n";
####This part removes the SNPs that appears more than one time (SNPs that have the same id)
open(FILE1,"${name[0]}.tped");
open(OUT,">temp_${name[0]}.tped");

%ids=();

while (<FILE1>){
	chomp;
	@split=split(/ /,$_);
	
	if(!exists($ids{$split[1]})){
		$ids{$split[1]} = 1; 
		print OUT "$_\n";
	}
	else{
		next;
	}	
}
close(FILE1);
close(OUT); 

system("cp ${name[0]}.tfam temp_${name[0]}.tfam");

###Now we have to remove one of the SNPs that have the same physical position

open(FILE1,"duplicate_positions_${name[0]}.txt");
open(OUT,">to_remove.txt");

while(<FILE1>){
	chomp;
	@split=split(/ /,$_);
	if(scalar(@split)==5){
		print OUT "$split[4]\n";
	}
	if(scalar(@split)>5){
		print OUT "$split[3]\n";
		print OUT "$split[4]\n";
	}
}
close(FILE1);
close(OUT);		

#system("cut -f4 -d' ' duplicate_positions_${name[0]}.txt > to_remove.txt");
system("plink --tfile temp_${name[0]} --exclude to_remove.txt --allow-no-sex --recode transpose --out temp2_${name[0]}");

print "\n\n\n***CAMBIADO CODIGO DE ALELOS***\n\n\n";

###Use the list of SNPs that comes from get_snps_to_change_alleles.pl that prints the list of the SNPs whose codes have to be changed

###With this function, we change the name of the alleles
###This final file does not have SNPs with the same id, or with duplicates with the same position
system("plink --tfile temp2_${name[0]} --update-alleles list_snps_to_change_alleles.txt --update-map new_positions_snps.txt --allow-no-sex --recode transpose --out temp3_${name[0]}");
system("plink --tfile temp3_${name[0]} --recode --out temp3_${name[0]}");

print "\n\n\n***RENOMBRANDO SNPS***\n\n\n";

###Now we have to rename the SNPs IDs 
system (qq(awk 'NR==FNR{if(\$2 != "." && \$2 !~ /,/){a[\$1]=\$2};next}{if(\$2 in a)sub(\$2,a[\$2]);print}' OFS="\t" /data/users/coquinto/test_langebio/MEGA_Consortium_v2_15070954_A1_b138_rsids.txt temp3_${name[0]}.map > temp3_${name[0]}_renamed.map));

system("mv temp3_${name[0]}_renamed.map temp3_${name[0]}.map");
system("plink --file temp3_${name[0]} --recode transpose --out temp3_${name[0]}");

###Have to remove duplicates after the renaming of the sites
open(FILE2,"temp3_${name[0]}.tped");
open(OUT,">temp4_${name[0]}.tped");

%ids=();

while (<FILE2>){
	chomp;
	@split=split(/ /,$_);
	
	if(!exists($ids{$split[1]}) && $split[4] ne "0"){
		$ids{$split[1]} = 1; 
		print OUT "$_\n";
	}	
}
close(FILE2);
close(OUT); 

system("cp temp3_${name[0]}.tfam temp4_${name[0]}.tfam");
system("plink --tfile temp4_${name[0]}  --update-alleles list_snps_to_change_alleles.txt --allow-no-sex --recode transpose --out temp5_${name[0]}");

system("plink --tfile temp5_${name[0]} --recode --out $name[0].flip");
#system("mv temp5_${name[0]}.tped ${name[0]}.flip.tped");
#system("mv temp5_${name[0]}.tfam ${name[0]}.flip.tfam");

system("rm temp*_*");
system("rm *.log");
system("rm *.nosex");
#system("rm to_remove.txt");

print "\n\n\n***LISTO!!!***\n\n\n";
