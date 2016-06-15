#!/usr/local/bin/perl

use warnings; 

$ped=$ARGV[0];
@name=split(/\.ped/,$ped);

print "\n\n\n***INICIANDO PIPELINE***\n\n\n";


###Order the plink files by chromosome number and transpose the data
system("plink2 --file $name[0] --allow-no-sex --recode --not-chr 0 --out $name[0]");
system("plink2 --file $name[0] --allow-no-sex --recode transpose --not-chr 0 --out $name[0]");

print "\n\n\n***BUSCANDO DUPLICADOS***\n\n\n";

###Find duplicates using plink
system("plink2 --file $name[0] --list-duplicate-vars --out duplicates_plink");
system("sed -i 1d duplicates_plink.dupvar");

print "\n\n\n***IMPRIMIENDO INFO DE DUPLICADOS***\n\n\n";
print "\n\n\n***IGNORAR ERRORES DE PERL***\n\n\n";

###Get summary of the duplicates
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

####This part removes the SNPs that appears more than one time (SNPs that have the same id)
print "\n\n\n***ELIMINANDO SNPS DUPLICADOS***\n\n\n";
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

system("plink2 --tfile temp_${name[0]} --exclude to_remove.txt --allow-no-sex --recode transpose --out temp2_${name[0]}");

print "\n\n\n***CAMBIADO CODIGO DE ALELOS***\n\n\n";

###Use the list of SNPs that comes from get_snps_to_change_alleles.pl that prints the list of the SNPs whose codes have to be changed
###With this function, we change the name of the alleles
###This final file does not have SNPs with the same id, or with duplicates with the same position
system("plink2 --tfile temp2_${name[0]} --update-alleles list_snps_to_change_alleles.txt --update-map new_positions_snps.txt --allow-no-sex --recode transpose --out temp3_${name[0]}");
system("plink2 --tfile temp3_${name[0]} --recode --out temp3_${name[0]}");

print "\n\n\n***RENOMBRANDO SNPS***\n\n\n";
###Now we have to rename the SNPs IDs 
system (qq(awk 'NR==FNR{if(\$2 != "." && \$2 !~ /,/){a[\$1]=\$2};next}{if(\$2 in a)sub(\$2,a[\$2]);print}' OFS="\t" /data/users/coquinto/test_langebio/MEGA_Consortium_v2_15070954_A1_b138_rsids.txt temp3_${name[0]}.map > temp3_${name[0]}_renamed.map));

system("mv temp3_${name[0]}_renamed.map temp3_${name[0]}.map");
system("plink2 --file temp3_${name[0]} --recode transpose --out temp3_${name[0]}");

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

###Update the new files
system("cp temp3_${name[0]}.tfam temp4_${name[0]}.tfam");
system("plink2 --tfile temp4_${name[0]}  --update-alleles list_snps_to_change_alleles.txt --allow-no-sex --recode transpose --out temp5_${name[0]}");

system("plink2 --tfile temp5_${name[0]} --recode --out $name[0].flip");
 
system("rm temp*_*");
system("rm *.log");
system("rm *.nosex");

###Calculate concordance with the controls
print "\n\n\n***CHECAR CONTROLES***\n\n\n";
@controls_1KG=(NA12878,HG01938,HG01941,NA19088,HG01935);
@controls_HP=(NA21402,NA12878,NA21405,NA21685,NA19088);

open(FILE,"$name[0].flip.ped"); ##or $name[0].flip.fam
while(<FILE>){
	chomp;
	foreach $control (@controls_1KG){
		$_=~/$control/;
		$flag="1KG";
	}
	foreach $control (@controls_HP){
		$_=~/$control/;
		$flag="HP";
	}		
}
close(FILE);

###Check the dataset
if ($flag eq "1KG"){
	print "$flag\n";
	
	#find common snps with 1KG_Affy6.0_training_samples
	%RS_NUMBER=();
	
	@files=("$name[0].flip.map","1KG_Affy6.0_training_samples.map");
	
	foreach $file (@files){
   		open (INFILE,"<$file" ) || die ( "Cannot open the file  $file  for reading: $!" );
		while(<INFILE>){
        		@split=split(/\t/,$_);
        		$rsNumber = $split[1];
        	
			if(!exists($RS_NUMBER{$rsNumber})){ $RS_NUMBER{$rsNumber}=1;}
			else { $RS_NUMBER{$rsNumber}++; }
   		}
   		close(INFILE);
	}
	
	open(OUT,">common_snps_1KG.txt");
	foreach $rsNumber (sort keys %RS_NUMBER){
		if( $RS_NUMBER{$rsNumber} == scalar(@files) ){ print OUT "$rsNumber\n"; }
	}
	
	##make the merge with plink 
	system("plink2 --file 1KG_Affy6.0_training_samples --extract common_snps_1KG.txt --make-bed --out 	1KG_Affy6.0_training_samples_overlap");
	system("plink2 --file $name[0].flip --extract common_snps_1KG.txt --make-bed --out $name[0].flip_overlap");
	
	system("plink2 --bfile 1KG_Affy6.0_training_samples_overlap --bmerge $name[0].flip_overlap.bed $name[0].flip_overlap.bim $name[0].flip_overlap.fam --recode --out merged_data");
		           
}

if ($flag eq "HP"){
	print "$flag\n";
	#hapmap3_b37_autosomes_training_samples.map
	
	%RS_NUMBER=();
	
	@files=("$name[0].flip.map","hapmap3_b37_autosomes_training_samples.map");
	
	foreach $file (@files){
   		open (INFILE,"<$file" ) || die ( "Cannot open the file  $file  for reading: $!" );
		while(<INFILE>){
        		chomp;
        		@split=split(/\t/,$_);
        		$rsNumber = $split[1];

			if(!exists($RS_NUMBER{$rsNumber})){ $RS_NUMBER{$rsNumber}=1;}
			else { $RS_NUMBER{$rsNumber}++; }
   		}
   		close(INFILE);
	}
	
	open(OUT,">common_snps_HP.txt");
	foreach $rsNumber (sort keys %RS_NUMBER){
		if( $RS_NUMBER{$rsNumber} == 2 ){ print OUT "$rsNumber\n"; }
	}
	
	##make the merge with plink 
	system("plink2 --file hapmap3_b37_autosomes_training_samples --extract common_snps_HP.txt --make-bed --out 	hapmap3_b37_autosomes_training_samples_overlap");
	system("plink2 --file $name[0].flip --extract common_snps_HP.txt --make-bed --out $name[0].flip_overlap");
	
	system("plink2 --bfile hapmap3_b37_autosomes_training_samples_overlap --bmerge $name[0].flip_overlap.bed $name[0].flip_overlap.bim $name[0].flip_overlap.fam --recode --out merged_data");	
	
	system("plink2 --bfile hapmap3_b37_autosomes_training_samples_overlap --bmerge $name[0].flip_overlap.bed $name[0].flip_overlap.bim $name[0].flip_overlap.fam --exclude-snps merged_data.missnp --recode --out merged_data");	

	
	#plink2 --bfile entrenamiento_16sample_nomissing_noduplicates2_overlap --flip merged_data.missnp --make-bed --out source2

###Merge again 
#plink2 --bfile 1KG_Affy6.0_training_samples_overlap --bmerge source2.bed source2.bim source2.fam --out merged_data2

	
}

##########
print "\n\n\n***LISTO!!!***\n\n\n";
