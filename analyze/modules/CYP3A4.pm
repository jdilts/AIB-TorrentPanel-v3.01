##############################################################
### CYP3A4.pm                                              ###
### April 14, 2014                                         ###
### AIBiotech LLC                                          ###
### William Budd; PhD                                      ###
### Version: 2.00                                          ###
### Improvements Needed:                                   ###
##############################################################
package CYP3A4;
use strict;
use warnings;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(callCYP3A4);

#-----PROTOTYPE-----
sub callCYP3A4($$$);
sub buildVariations();
sub calculate_ratios($$);
sub checkPositions($$);
sub callGeno($$);
sub prioritizeMutations($);
sub Format($$);
#-------------------

##### This method accepts three inputs, reference to sample Hash, number of mapped reads from alignment barcode summary file, and the minimum SNP coverage.
##### Checks for missing data after building a hash with the variants of interest for the gene. If reads are less than 5500, any allele is not covered or covered
##### with less than 20X coverage, the sample errors out and does not proceed genotyping. For each position of interest, the method calls other sub-routines to
##### examine the genetic makeup at that position and returns a reference to an array with the ratios of each nucleotide. 
###### Case 1 = 0 mut and >0 het 
###### Case 2 = >0 mut and 0 het
###### Case 3 = >0 mut and >0 het
###### Case 4 = 0 mut and 0 het ->  WT 
sub callCYP3A4($$$){
        my($temp_sample,$reads,$snpcov) = @_;	
        my %temp_sample = %$temp_sample;
		#-----VARIABLE DECLARATION-----
        my ($result, $returnRef, $variantRef, $SNPerror, $message, $variant, $star2, $call, $starRef, $star, $value, $position,
			$finalCall, $SNP, $priorityRef, $call1, $call2, $allele1, $allele2);
		my (@ordered, @genoCalls, @het, @mut, @mutations, @heteros, @stack, @alleles);
		my (%variations, %starCalls, %priorityHash);      	
		@ordered = ("CYP3A4s12,60","CYP3A4s17,84","CYP3A4s1B,44","CYP3A4s2,105","CYP3A4s3,80");	
		
		########### Main method of callCYP2C9 
		$variantRef= buildVariants();
		$starRef = buildStarCalls();
		%variations = %$variantRef;
		%starCalls = %$starRef;	
		### Checks for errors in read coverage and ensures that each position has reads	
		$message="";
		if($reads <5500){
			$message = "Error: Reads less than 5500";
		}
		$SNPerror = checkPositions($variantRef, $temp_sample);
		if($SNPerror ne""){
			$message = "$SNPerror";
			return "CYP3A4|$message||||||||\n";
		}
        
		### Determines if each position is a Het, Mut or WT. This also determines if each position exceeds 20X coverage at each allele. 
		### If not at least 20X coverage, issues an error message. 
		foreach my $key (@ordered){	   			
			$variant = $variations{$key};			
			($returnRef) = calculate_ratios($temp_sample{$key},$snpcov);
			$call = callGeno($variant, $returnRef);
			$star = $starCalls{$key};
			$value = $star." ".$call;
			if($call =~/Error/){
				$message = $call;
				return "CYP3A4|$message $key||||||||\n";
			}			
			push(@genoCalls, $value);			
		}# end call WT, HET, and MUT
				
		foreach $position(@genoCalls){			
			if($position =~/Het/){
				push(@het, $position);
			}
			elsif($position =~/Mut/){
				push (@mut, $position);
			}			
		}# end build Call arrays
		
		###### Case 1 = 0 mut and > 0 het
		#Call heterogeneous genotype 
		if (scalar @het >0 && scalar @mut==0){		
			foreach my $element (@het){
				($SNP) = $element=~/(\d+\w*)\s*Het/;
				push(@heteros, $SNP);
			}
			$priorityRef = prioritizeMutations(\@heteros);
			%priorityHash = %$priorityRef;
			if(exists $priorityHash{2}){
				$finalCall = Format($priorityHash{1}, $priorityHash{2});			
			}
			elsif( not exists $priorityHash{2}){
				$finalCall = "*1|*$priorityHash{1}"; 
			}
		}
		
		###### Case 2 = >0 mut and 0 het
		## Call mutant genotype
		if(scalar@mut >0 && scalar @het==0){
			foreach my $element (@mut){		
				($SNP) = $element=~/(\d+\w*)\s*Mut/;
				push(@mutations, $SNP);
			}	
			$priorityRef = prioritizeMutations(\@mutations);
			%priorityHash = %$priorityRef;			
			$finalCall = "*$priorityHash{1}|*$priorityHash{1}"; 
		}
		
		#### Case 3 -> Distributed and prioritized if combination of het and mutant calls. 
		elsif(scalar @het>0 && scalar @mut>0){		
			foreach my $element (@het){			
				push(@stack,$element);
			} 			
			foreach my $element2 (@mut){
				push(@stack, $element2);
				push(@stack, $element2);
			}		
			foreach my $element(@stack){
				($SNP) = $element=~/(\d+\w*)\s*\w+/;
				push(@alleles, $SNP);
			}	
			$priorityRef = prioritizeMutations(\@alleles);
			%priorityHash = %$priorityRef;	
			$call1 = $priorityHash{1};
			$call2 = $priorityHash{2};			
			$finalCall = Format($call1, $call2);
		}# end Case 3
	
		###### Case 4 = 0 mut and 0 het ->  WT 	
		elsif (scalar @het ==0 && scalar @mut==0){
			$finalCall = "*1|*1" 
		}
		
		###### Return statements to Genotyper for printing results ! 
		if($message eq ""){
			return "CYP3A4|$message|$finalCall|||||||\n";
		}
		else{
			return "CYP3A4|$message||||||||\n";
		}
}# end callCYP3A4             

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
%variations = (		"CYP3A4s12,60" => "T",
					"CYP3A4s17,84"=> "C",
					"CYP3A4s1B,44"=> "G",
					"CYP3A4s2,105"=> "G",
					"CYP3A4s3,80"=> "C",			
				);
	return(\%variations);
} # end buildVariants

##### buildStarCalls builds a hash that transforms each position into the appropriate genetic designation associated with each position. Example "*7, *9 etc."
sub buildStarCalls(){
my(%variations);
%variations = (		"CYP3A4s12,60" => "12",
					"CYP3A4s17,84"=> "17",
					"CYP3A4s1B,44"=> "1B",
					"CYP3A4s2,105"=> "2",
					"CYP3A4s3,80"=> "3",				
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
					"5" => "-",
				);	
	for(my $i=1; $i<scalar @ratios; $i++){
		if($conversion{$i} eq $variant){
			if($ratios[$i]>=0.8){
				$genotype = "Mut";
			}
			elsif($ratios[$i] >=0.20){
				$genotype = "Het";
			}
			else{
				$genotype = "WT";
			}
		}
	}
	return $genotype;
}# end callGeno

###### prioritizeMutations accepts a reference to an array of converted star calls. It establishes a priority hash then loops over the array of converted calls
###### The priorities are identified and pushed into an array. After all calls are prioritized, the array of priorities is sorted in numerical order. 
###### After sorting the ordered array, a new hash is created from 1 through the number of calls and the value is the converted call. 
###### Example input: 17, 3, 2, 9
###### Example output: 9, 3, 2, 17  
sub prioritizeMutations($){
my ($ref) = @_;
my @array = @$ref;
## Variable declaration
my($query, $element, $counter, $pri);
my (@ordered, @sorted);
my (%hash);
my %priorityHash =( "12"  => "1",
					"17" => "2",
					"1B"  => "3",
					"2"  => "4",
					"3"  => "5",													
					);
					
my %reversedHash = ( "1"  => "12",
					"2" => "17",
					"3"  => "1B",
					"4"  => "2",
					"5"  => "3",													
					);					
$counter=1;
foreach $element(@array){
	$pri = $priorityHash{$element};
	push(@ordered, $pri);	
}

@sorted = sort {$a <=> $b} @ordered;
for(my $i=1; $i<= scalar @sorted; $i++){
	$hash{$i}=$reversedHash{$sorted[$i-1]};
	}  	
return \%hash;

}

##### sub Format is a subroutine that orients the genotyped calls into the correct format. Example *2 *17 not *17 *2 The sub accepts two scalar values and sorts
##### them numerically and formats them into a genotype call and returns to the main block. 
sub Format($$){
my($call1, $call2) = @_;
my (@unordered, @ordered, $call);
push(@unordered, $call1);
push(@unordered, $call2);
@ordered = sort {$a <=> $b} @unordered;
$call = "*$ordered[0]|*$ordered[1]"; 
return $call;
}# end checkFormatting 

1; #return statement
