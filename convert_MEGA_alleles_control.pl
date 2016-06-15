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
system("sed -i 's/ /\t/g' duplicates_plink.dupvar");

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
	
	#print "$split[5]\n";
	
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
system("plink2 --tfile temp5_${name[0]} --make-bed --out $name[0].flip");

###Calculate concordance with the controls
print "\n\n\n***REVISANDO CONTROLES***\n\n\n";

###Remove the sites that are monomorphic, and multiallelic
open(TXT,">monomorphic.txt");

open(BIM,"$name[0].flip.bim");
while(<BIM>){
	chomp;
	@split=split(/\t/,$_);
	if ($split[4] eq "0"){
		print TXT "$split[1]\n";		
	}
}
close(BIM);
close(TXT);

######

if (-s "monomorphic.txt"){
	print "archivo no vacio";
	system("plink2 --bfile $name[0].flip --exclude monomorphic.txt --make-bed --out $name[0].flip_mono");
	system("plink2 --bfile $name[0].flip_mono --exclude multiallelic_SNPs.txt --make-bed --out $name[0].flip_mono_nomulti");

	####find common snps with 1KG_Affy6.0_training_samples
	%RS_NUMBER=();
	
	@files=("$name[0].flip_mono_nomulti.bim","1KG_Omni.bim");
	
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
	
	###Make the merge with plink 
	system("plink2 --bfile 1KG_Omni --extract common_snps_1KG.txt --make-bed --out 1KG_Omni_overlap");
	system("plink2 --bfile $name[0].flip_mono_nomulti --extract common_snps_1KG.txt --make-bed --out $name[0].flip_mono_nomulti_overlap");
	
	system("plink2 --bfile 1KG_Omni_overlap --bmerge $name[0].flip_mono_nomulti_overlap.bed $name[0].flip_mono_nomulti_overlap.bim $name[0].flip_mono_nomulti_overlap.fam --recode --out merged_data");

	print "\n\nSi no hay un archivo merged_data.bed; entonces hay un problema!!!\n\n";

	####Find controls in merged_data
	@controls_1KG=(HG01946,HG01935,NA21440,NA21317);

	%ids_uniq=();

	open(PED,"merged_data.fam");
	while(<PED>){
		chomp;
		@split2=split(/\t/,$_);
		for $id (@controls_1KG){
			if ($split2[1] =~ /$id/){
				if(!exists($ids_uniq{$id})){ $ids_uniq{$id}=1;}
				else { $ids_uniq{$id}++; }
			}	
		}
	}
	close(PED);

	$control="";
	foreach $element (sort keys %ids_uniq){
		print "$element\t$ids_uniq{$element}\n";
		if( $ids_uniq{$element} >= 2 ){
			print "$element\n";
			$control=$element;
		}
	}
	print "CONTROL $control\n";
	
	###Calculate concordance
	system("sed -nr '/${control}/p' merged_data.ped > merged_data_2.ped");
	system("cp merged_data.map merged_data_2.map");
	system("plink2 --file merged_data_2 --distance allele-ct square ");

	$difference=`cut -f1 plink.dist | tail -1`;
	chomp($difference);
	$total=`wc -l merged_data_2.map | cut -d' ' -f1`;
	chomp($total);
	$concordance=1-($difference/$total);

	print "\n\n\n***CONCORDANCIA CON $control***\n";
	print "Concordance $concordance\n";
	print "Si la concordancia es menor de 0.90, hay un problema!!!\n";

}	

else{ 
	
	system("plink2 --bfile $name[0].flip --exclude multiallelic_SNPs.txt --make-bed --out $name[0].flip_nomulti");

	####find common snps with 1KG_Affy6.0_training_samples
	%RS_NUMBER=();
	
	@files=("$name[0].flip_nomulti.bim","1KG_Omni.bim");
	
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
	
	###Make the merge with plink 
	system("plink2 --bfile 1KG_Omni --extract common_snps_1KG.txt --make-bed --out 1KG_Omni_overlap");
	system("plink2 --bfile $name[0].flip_nomulti --extract common_snps_1KG.txt --make-bed --out $name[0].flip_nomulti_overlap");
	
	system("plink2 --bfile 1KG_Omni_overlap --bmerge $name[0].flip_nomulti_overlap.bed $name[0].flip_nomulti_overlap.bim $name[0].flip_nomulti_overlap.fam --recode --out merged_data");

	print "\n\nSi no hay un archivo merged_data.bed; entonces hay un problema!!!\n\n";

	####Find controls in merged_data
	@controls_1KG=(HG01946,HG01935,NA21440,NA21317);

	%ids_uniq=();

	open(PED,"merged_data.fam");
	while(<PED>){
		chomp;
		@split2=split(/\t/,$_);
		for $id (@controls_1KG){
			if ($split2[1] =~ /$id/){
				if(!exists($ids_uniq{$id})){ $ids_uniq{$id}=1;}
				else { $ids_uniq{$id}++; }
			}	
		}
	}
	close(PED);

	$control="";
	foreach $element (sort keys %ids_uniq){
		#print "$element\t$ids_uniq{$element}\n";
		if( $ids_uniq{$element} >= 2 ){
			#print "$element\n";
			$control=$element;
		}
	}
	print "CONTROL $control\n";
	
	###Calculate concordance
	system("sed -nr '/${control}/p' merged_data.ped > merged_data_2.ped");
	system("cp merged_data.map merged_data_2.map");
	system("plink2 --file merged_data_2 --distance allele-ct square ");

	$difference=`cut -f1 plink.dist | tail -1`;
	chomp($difference);
	$total=`wc -l merged_data_2.map | cut -d' ' -f1`;
	chomp($total);
	$concordance=1-($difference/$total);

	print "\n\n\n***CONCORDANCIA CON $control***\n";
	print "Concordance $concordance\n";
	print "Si la concordancia es menor de 0.90, hay un problema!!!\n";
}	
		           
##########
print "\n\n\n***LISTO!!!***\n\n\n";

##########
system("rm temp*_*");
system("rm *.log");
system("rm *.nosex");
system("rm to_remove.txt");
system("rm merged_data*");
system("rm plink.*");
system("rm *overlap.*");
system("rm common_snps_1KG.txt");
system("rm duplicates_plink.dupvar");
system("rm monomorphic.txt");
system("rm *_mono*");
system("rm $name[0].tped");
system("rm $name[0].tfam");
system("rm $name[0].flip.bed");
system("rm $name[0].flip.bim");
system("rm $name[0].flip.fam");

