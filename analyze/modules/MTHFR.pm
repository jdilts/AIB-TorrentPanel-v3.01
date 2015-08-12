##############################################################
### MTHFR.pm                                               ###
### April 14, 2014                                         ###
### AIBiotech LLC                                          ###
### William Budd; PhD                                      ###
### Version: 2.00                                          ###
### Improvements Needed:                                   ###
##############################################################
package MTHFR;
use strict;
use warnings;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(callMTHFR callMTHFR677 callMTHFR1298);

#-----PROTOTYPE-----
sub callMTHFR($$$);
sub callMTHFR677($$$);
sub callMTHFR1298($$$);
sub buildVariations();
sub calculate_ratios($$);
sub checkPositions($$);
sub callGeno($$);

#-------------------

##### This method accepts three inputs, reference to sample Hash, number of mapped reads from alignment barcode summary file, and the minimum SNP coverage.
##### Checks for missing data after building a hash with the variants of interest for the gene. If reads are less than 5500, any allele is not covered or covered
##### with less than 20X coverage, the sample errors out and does not proceed genotyping. For each position of interest, the method calls other sub-routines to
##### examine the genetic makeup at that position and returns a reference to an array with the ratios of each nucleotide. 

sub callMTHFR($$$){
        my($temp_sample,$reads,$snpcov) = @_;	
        my %temp_sample = %$temp_sample;
		#-----VARIABLE DECLARATION-----
        my ($result, $returnRef, $variantRef, $SNPerror, $message, $variant, $call,$finalCall);
		my (@ordered, @genoCalls);
		my (%variations);      	
		@ordered = ("MTHFRC677T,80","MTHFRA1298C,61");	
		
		########### Main method of callMTHFR 
		$variantRef= buildVariants();		
		%variations = %$variantRef;
		
		### Checks for errors in read coverage and ensures that each position has reads	
		$message="";
		if($reads <5500){
			$message = "Error: Reads less than 5500";
		}
		$SNPerror = checkPositions($variantRef, $temp_sample);
		if($SNPerror ne""){
			$message = "$SNPerror";
			return "MTHFR|$message||||||||\n";
		}
        
		### Determines if each position is a Het, Mut or WT. This also determines if each position exceeds 20X coverage at each allele. 
		### If not at least 20X coverage, issues an error message. 
		foreach my $key (@ordered){	   		
			$variant = $variations{$key};			
			($returnRef) = calculate_ratios($temp_sample{$key},$snpcov);
			$call = callGeno($variant, $returnRef);			
			if($call =~/Error/){
				$message = $call;
				return "MTHFR|$message $key||||||||\n";
			}
			else{
				push(@genoCalls, $call);			
			}
		}# end call WT, HET, and MUT
		
		
		$finalCall = $genoCalls[0]."|".$genoCalls[1];
		
		
		###### Return statements to Genotyper for printing results ! 
		if($message eq ""){
			return "MTHFR|$message|$finalCall|||||||\n";
		}
		else{
			return "MTHFR|$message||||||||\n";
		}
}# end callMTHFR           


sub callMTHFR677($$$){
        my($temp_sample,$reads,$snpcov) = @_;	
        my %temp_sample = %$temp_sample;
		#-----VARIABLE DECLARATION-----
        my ($result, $returnRef, $variantRef, $SNPerror, $message, $variant, $call,$finalCall);
		my (@ordered, @genoCalls);
		my (%variations);      	
		@ordered = ("MTHFRC677T,80");	
		
		########### Main method of callMTHFR 
		$variantRef= buildVariants();		
		%variations = %$variantRef;
		
		### Checks for errors in read coverage and ensures that each position has reads	
		$message="";
		if($reads <5500){
			$message = "Error: Reads less than 5500";
		}
		$SNPerror = checkPositions($variantRef, $temp_sample);
		if($SNPerror ne""){
			$message = "$SNPerror";
			return "MTHFRC677T|$message||||||||\n";
		}
        
		### Determines if each position is a Het, Mut or WT. This also determines if each position exceeds 20X coverage at each allele. 
		### If not at least 20X coverage, issues an error message. 
		   		
		$variant = $variations{$key};			
		($returnRef) = calculate_ratios($temp_sample{$key},$snpcov);
		$call = callGeno($variant, $returnRef);			
		if($call =~/Error/){
			$message = $call;
			return "MTHFRC677T|$message $key||||||||\n";
		}
		else{
			push(@genoCalls, $call);			
		}		
		$finalCall = $genoCalls[0];		
		
		###### Return statements to Genotyper for printing results ! 
		if($message eq ""){
			return "MTHFRC677T|$message|$finalCall|||||||\n";
		}
		else{
			return "MTHFRC677T|$message||||||||\n";
		}
}# end callMTHFR677

sub callMTHFR1298($$$){
        my($temp_sample,$reads,$snpcov) = @_;	
        my %temp_sample = %$temp_sample;
		#-----VARIABLE DECLARATION-----
        my ($result, $returnRef, $variantRef, $SNPerror, $message, $variant, $call,$finalCall);
		my (@ordered, @genoCalls);
		my (%variations);      	
		@ordered = ("MTHFRA1298C,61");	
		
		########### Main method of callMTHFR 
		$variantRef= buildVariants();		
		%variations = %$variantRef;
		
		### Checks for errors in read coverage and ensures that each position has reads	
		$message="";
		if($reads <5500){
			$message = "Error: Reads less than 5500";
		}
		$SNPerror = checkPositions($variantRef, $temp_sample);
		if($SNPerror ne""){
			$message = "$SNPerror";
			return "MTHFRA1298C|$message||||||||\n";
		}
        
		### Determines if each position is a Het, Mut or WT. This also determines if each position exceeds 20X coverage at each allele. 
		### If not at least 20X coverage, issues an error message. 
		   		
		$variant = $variations{$key};			
		($returnRef) = calculate_ratios($temp_sample{$key},$snpcov);
		$call = callGeno($variant, $returnRef);			
		if($call =~/Error/){
			$message = $call;
			return "MTHFRA1298C|$message $key||||||||\n";
		}
		else{
			push(@genoCalls, $call);			
		}		
		$finalCall = $genoCalls[0];		
		
		###### Return statements to Genotyper for printing results ! 
		if($message eq ""){
			return "MTHFRA1298C|$message|$finalCall|||||||\n";
		}
		else{
			return "MTHFRA1298C|$message||||||||\n";
		}
}# end callMTHFR1298


##### Calculate rations accepts a reference fo an array that contains nucleotide counts for each position. This method transforms the counts to a proportion/ ratio
##### The method also ensures that the coverage exceeds 20 fold for each position of interest. It returns a reference to an array with a flag and all ration
##### flag, A, C, G, T in order of nucleotide freqeuncy table 
sub calculate_ratios($$){
my ($ref,$snpcov) = @_;
my ($ratioA,$ratioC,$ratioG,$ratioT,$ratioDEL,$reads, $reference);
my (@return);
my $VarCovFLAG = 0; #variant coverage Flag
my @VAR = @$ref;
        #checks coverage	
        if ($VAR[1] >= $snpcov){
                $VarCovFLAG = 1;
                $ratioA = ($VAR[2]/$VAR[1]);
                $ratioC = ($VAR[3]/$VAR[1]);
                $ratioG = ($VAR[4]/$VAR[1]);
                $ratioT = ($VAR[5]/$VAR[1]);
                $ratioDEL = ($VAR[6]/$VAR[1]);
        }
		@return = ($VarCovFLAG, $ratioA, $ratioC, $ratioG, $ratioT, $ratioDEL);
        return (\@return);
}# end calculate_ratios 

##### sub buildVariants builds a hash of the alternative alleles for each position of interest. A reference to the hash is created and returned.
sub buildVariants(){
my(%variations);
%variations = (
					"MTHFRC677T,80" => "T",	
					"MTHFRA1298C,61" => "C",
				);
	return(\%variations);
} # end buildVariants


##### checkPositions accepts two hash references. One containing all positions of interest and the other containing the positions sequenced. 
##### This method ensures that all positions are included in the nucelotide_frequency table. It is a QA/ QC check. 
sub checkPositions($$){
my ($variantRef, $sampleRef) = @_;
##variable declarations ##
my($key, $error);
my (%variants, %samples);
	%variants = %$variantRef;
	%samples= %$sampleRef;
	$error="";
	foreach $key (keys %variants){
		if(! exists $samples{$key}){
			$error = "No call for $key";
		}		
	}
return $error;
}# end sub checkPositions

##### callGeno method accepts the nucleotide variant and a reference to an array containing the ratio/ proportions of the calls at each position.
##### Method is set up such that Position 0 is the less than 20X flag, 1-> A proportion, 2-> C proportion, 3-> G proportion, 4-> T proportion, 5-> deletion proportion
##### Loops over the array of ration and focuses on the variant(s) of interest. If variant proportion is >= 0.8 , allele called mutant
##### Allele proportions >-0.2 and <0.8 are heterozygous. All others are wild type. Call is returned to main method. 
sub callGeno($$){
my ($variant, $ratioRef) = @_;
## variable declarations 
my ($genotype, $A, $T, $C, $G, $del, $flag, $call);
my (@ratios);
my (%conversion);	
	$genotype="";
	@ratios = @$ratioRef;
	if($ratios[0]==0){return "Error: less than 20 reads";}
	%conversion=(	"1" => "A",
					"2" => "C",
					"3" => "G",
					"4" => "T",
					"5" => "del",
				);

	for(my $i=1; $i<scalar @ratios; $i++){
		if($conversion{$i} eq $variant){
			if($ratios[$i]>=0.8){
				$genotype = "MT";
			}
			elsif($ratios[$i] >=0.20){
				$genotype = "HET";
			}
			else{
				$genotype = "WT";
			}
		}
	}
	return $genotype;
}# end callGeno




1; #return statement
