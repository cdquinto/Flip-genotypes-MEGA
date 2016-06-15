#!/usr/bin/perl
use warnings;
use File::Slurp;

$file=$ARGV[0];
#$file="/home/cdquinto/test_langebio/analysis_entrenamiento_nocluster/concordance_hapmap3_MEGA/merged_data_hapmap_3_ok.ped"; 

@array=read_file($file);

open(OUT,">results.txt");
open(OUT2,">resultsR.txt");

#@matrix = ([0,6],[1,4],[2,5],[3,7] ); 
#@matrix=([0,7],[1,5],[2,8],[3,6],[4,9]);
#@matrix=([0,2]);
@matrix=([0,1]);

#for($i=0;$i<4;$i++){
#for($i=0;$i<5;$i++){
for($i=0;$i<1;$i++){
	chomp $array[$matrix[$i][0]];
	chomp $array[$matrix[$i][1]];
	#print $array[$matrix[$i][1]],"\n";
	@split1=split(/\-9/,$array[$matrix[$i][0]]);
	$split1[1] =~ s/ //g;
	@string1=split(//,$split1[1]);
	print scalar(@string1),"\n";

	@split2=split(/\-9/,$array[$matrix[$i][1]]);
	@name=split(/ /,$split2[0]);
	$split2[1] =~ s/ //g;
	@string2=split(//,$split2[1]);

	$cont=0;
	
	for($j=0;$j<scalar(@string1);$j++){
  		if($string1[$j] eq $string2[$j]){
  			$cont++;
  		}
  		elsif ($string1[$j] eq 'G' && $string2[$j] eq 'C'){
  			$cont++;
  		}
  		elsif ($string1[$j] eq 'C' && $string2[$j] eq 'G'){
  			$cont++;
  		}
  		elsif ($string1[$j] eq 'T' && $string2[$j] eq 'A'){
  			$cont++;
  		}
  		elsif ($string1[$j] eq 'A' && $string2[$j] eq 'T'){
  			$cont++;
  		}
	}
	
	$long=length($split1[1]);
	$ratio=$cont/$long;
	$minus_ratio=1-$ratio;
	#print OUT "$name[0]\t$matches\t$long\t$ratio\t$minus_ratio\n";
	#print OUT2 "$name[0]\t$ratio\t$minus_ratio\n";
	print OUT "$name[1]\t$cont\t$long\t$ratio\t$minus_ratio\n";
	print OUT2 "$name[1]\t$ratio\t$minus_ratio\n";
}

die;

open(RCMD,">commandsR.txt");
print RCMD "data <- read.delim(\"resultsR.txt\", header=FALSE, row.names=1)\n";
#print RCMD "pdf(\"Concordance_1KG_MEGA_clusterGLOBAL.pdf\")\n";
print RCMD "pdf(\"Concordance_HAPMAP_MEGA_clusterGLOBAL.pdf\")\n";
print RCMD "at_tick <- seq(0.1,0.9,0.2)\n";
#print RCMD "barplot(t(as.matrix(data)), ylim=c(0,1), col=c(\"chocolate1\",\"cadetblue1\"), main=\"Concordance between 1KG and MEGA\")\n";
print RCMD "barplot(t(as.matrix(data)), ylim=c(0,1), col=c(\"chocolate1\",\"cadetblue1\"), main=\"Concordance between HAPMAP3 and MEGA\")\n";
print RCMD "abline(h = mean(data\$V2), col = \"black\", lwd=5)\n";
print RCMD "axis(2, at = at_tick,tcl = -0.25)\n";
print RCMD "dev.off()\n";
close(RCMD);

system ( "R --slave < commandsR.txt" );
