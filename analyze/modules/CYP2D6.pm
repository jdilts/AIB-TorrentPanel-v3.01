##############################################################
### CYP2D6.pm                                              ###
### April 8, 2014                                          ###
### AIBiotech LLC                                          ###
### William Budd; PhD                                      ###
### Version: 2.00                                          ###
### Improvements Needed:                                   ###
##############################################################
package CYP2D6;
use strict;
use warnings;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(callCYP2D6);

#-----PROTOTYPE-----
sub callCYP2D6($$$);
sub buildVariations();
sub calculate_ratios($$);
sub checkPositions($$);
sub callGeno($$$);
sub prioritizeMutations($);
sub Format($$);
sub checkLinkage($);
sub check2($$$);
#-------------------

##### This method accepts three inputs, reference to sample Hash, number of mapped reads from alignment barcode summary file, and the minimum SNP coverage.
##### Checks for missing data after building a hash with the variants of interest for the gene. If reads are less than 5500, any allele is not covered or covered
##### with less than 20X coverage, the sample errors out and does not proceed genotyping. For each position of interest, the method calls other sub-routines to
##### examine the genetic makeup at that position and returns a reference to an array with the ratios of each nucleotide. 
###### Case 1 = 0 mut and >0 het 
###### Case 2 = >0 mut and 0 het
###### Case 3 = >0 mut and >0 het
###### Case 4 = 0 mut and 0 het ->  WT 
sub callCYP2D6($$$){
        my($temp_sample,$reads,$snpcov) = @_;	
        my %temp_sample = %$temp_sample;
		#-----VARIABLE DECLARATION-----
        my ($result, $returnRef, $variantRef, $SNPerror, $message, $variant, $star2, $call, $starRef, $star, $value, $position,
			$finalCall, $SNP, $priorityRef, $call1, $call2, $allele1, $allele2, $temp, $pos, $interpretation, $star9, $valid,
			$checkRef, $checkRef2);
		my (@ordered, @genoCalls, @het, @mut, @mutations, @heteros, @stack, @alleles, @calls);
		my (%variations, %starCalls, %priorityHash);      	
		@ordered = ("CYP2D6s10s12,99","CYP2D6s10s12,123","CYP2D6s17,74","CYP2D6s2s7s41,15","CYP2D6s2s7s41,100","CYP2D6s2s7s41,153",
					"CYP2D6s8s14s6s29,33","CYP2D6s8s14s6s29,81","CYP2D6s8s14s6s29,132","CYP2D6s3s9,78","CYP2D6s3s9,142","CYP2D6s3s9,143",
					"CYP2D6s3s9,144","CYP2D6s4,90");			
	
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
			return "CYP2D6|$message||||||||\n";
		}
		if($message ne ""){
			return "CYP2D6|$message||||||||\n";
		}
        
		### Determines if each position is a Het, Mut or WT. This also determines if each position exceeds 20X coverage at each allele. 
		### If not at least 20X coverage, issues an error message. 
		foreach my $key (@ordered){				
			$variant = $variations{$key};			
			($returnRef) = calculate_ratios($temp_sample{$key},$snpcov);		
			$call = callGeno($variant, $returnRef, $key);
			if($call =~/Error/){
				$message = $call;
				return "CYP2D6|$message $key||||||||\n";
			}					
			elsif($call =~/\|/){         ####### Used for positions that have two possible variations! 
				(@calls) = split(/\|/,$call);
				foreach $temp (@calls){
					($pos,$interpretation) = split(/\s+/,$temp);
					$star = $starCalls{$pos};
					$value = $star." ".$interpretation;
					push(@genoCalls, $value);
				}
			}
			else{
				$star = $starCalls{$key};
				$value = $star." ".$call;
				push(@genoCalls, $value);
			}					
		}# end call WT, HET, and MUT
#		print "@genoCalls\n";
		if($genoCalls[11] eq $genoCalls[12] && $genoCalls[12] eq $genoCalls[13]){ ## Processed star 9 variants
			$valid = 1;
			splice(@genoCalls,11,2);
		}
		else{			
			$message = "Manual review: Possible star 9 failure";
		}
#		print "@genoCalls\n";
		foreach $position(@genoCalls){			
			if($position =~/Het/){
				push(@het, $position);
			}
			elsif($position =~/Mut/){
				push (@mut, $position);
			}			
		}# end build Call arrays
		
		############## Linkage must be checked before making calls; consolidate linked calls to a single call 
		if(scalar @het>1){
			$checkRef = checkLinkage(\@het);
			@het = @$checkRef;
		}	
		if(scalar @mut>1){
			$checkRef2 = checkLinkage(\@mut);
			@mut = @$checkRef2;
		}
		
		###### Case 1 = 0 mut and > 0 het
		#Call heterogeneous genotype 
		if (scalar @het >0 && scalar @mut==0){		
			foreach my $element (@het){
				($SNP) = $element=~/(\d+)\s*Het/;
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
				($SNP) = $element=~/(\d+)\s*Mut/;
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
				($SNP) = $element=~/(\d+)\s*\w+/;
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
			return "CYP2D6|$message|$finalCall|||||||\n";
		}
		else{
			return "CYP2D6|$message||||||||\n";
		}
}# end callCYP2C19               

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
					"CYP2D6s10s12,99" => "T",
					"CYP2D6s10s12,123" => "A",
					"CYP2D6s17,74" => "T",
					"CYP2D6s2s7s41,15" => "T",
					"CYP2D6s2s7s41,100" => "C",
					"CYP2D6s2s7s41,153" => "A",
					"CYP2D6s8s14s6s29,33" => "A",
					"CYP2D6s8s14s6s29,81" => "-",
					"CYP2D6s8s14s6s29,132"=> "T A",
					"CYP2D6s3s9,78" => "-",
					"CYP2D6s3s9,142" =>"-",
					"CYP2D6s3s9,143" =>"-",
					"CYP2D6s3s9,144" => "-",	
					"CYP2D6s4,90" => "A",
				);
	return(\%variations);
} # end buildVariants

##### buildStarCalls builds a hash that transforms each position into the appropriate genetic designation associated with each position. Example "*7, *9 etc."
sub buildStarCalls(){
my(%variations);
%variations = (		"CYP2D6s10s12,99" => "10",
					"CYP2D6s10s12,123" => "12",
					"CYP2D6s17,74" => "17",
					"CYP2D6s2s7s41,15" => "2",
					"CYP2D6s2s7s41,100" => "7",
					"CYP2D6s2s7s41,153" => "41",
					"CYP2D6s8s14s6s29,33" => "29",
					"CYP2D6s8s14s6s29,81" => "6",
					"CYP2D6s8s14s6s29,132_T"=> "8",
					"CYP2D6s8s14s6s29,132_A"=> "14",					
					"CYP2D6s3s9,78" => "3",	
					"CYP2D6s3s9,142" =>"9",
					"CYP2D6s3s9,143" =>"9",					
					"CYP2D6s3s9,144" => "9",	
					"CYP2D6s4,90" => "4",
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
sub callGeno($$$){
my ($variant, $ratioRef, $key) = @_;
## variable declarations 
my ($genotype, $j,$A, $T, $C, $G, $del, $flag, $call, $variant1, $variant2, $WTfreq, $hetFreq);
my (@ratios);
my (%conversion);
	$variant1="";
	$variant2="";
	$genotype="";
	if($key =~/CYP2D6s3s9,14/){
		$WTfreq =0.7;
		$hetFreq = 0.3;
	}
	else{
		$WTfreq =0.8;
		$hetFreq=0.2;
	}
	@ratios = @$ratioRef;
	if($ratios[0]==0){return "Error: less than 20 reads";}
	%conversion=(	"1" => "A",
					"2" => "C",
					"3" => "G",
					"4" => "T",
					"5" => "-",
				);	
	## Consider variations that have only alternate nucleotide			
	if($variant !~/\s+/){			
		for(my $i=1; $i<scalar @ratios; $i++){
			if($conversion{$i} eq $variant ){
				if($ratios[$i]>= $WTfreq){
					$genotype = "Mut";
				}
				elsif($ratios[$i] >= $hetFreq){
					$genotype = "Het";
				}
				else{
					$genotype = "WT";
				}
			}
		}
		return $genotype;
	}# end one alternative nucleotide
	
	## Consider variations that have more than one alternative nucleotide
	if($variant =~/\s+/){
		 $genotype = check2($variant, \@ratios, $key);		
		return $genotype;
	} #end if more than one alternative nucleotide	
}# end callGeno

###### check2 subroutine is used to evaluate two possible variants for nucleotides that have more than one identified pathological variant.
######  Example: Position 132 of the CYP2D6s8s14s29 can be either a T or A 
sub check2($$$){
my ($variant, $ratioREF, $key) = @_;
## Variable declarations
my($variant1, $variant2, $genotype, $j, $call1, $call2, $returnCall);
my (@ratios, @returns);
my(%conversion);
	@ratios = @$ratioREF;
	%conversion=(	"1" => "A",
					"2" => "C",
					"3" => "G",
					"4" => "T",
					"5" => "-",
				);	
	($variant1, $variant2) = $variant =~/(\w)\s(\w)/;	
	for ($j=1; $j< scalar @ratios; $j++){
		if($conversion{$j} eq $variant2 ||$conversion{$j} eq $variant1){
#			print "$conversion{$j}\t$ratios[$j]\n";
			if($ratios[$j]>=0.8){
				$genotype = $key."_"."$conversion{$j} Mut";
				push(@returns,$genotype);
			}
			elsif($ratios[$j] >=0.20){
				$genotype = $key."_"."$conversion{$j} Het";
				push(@returns,$genotype);
			}
			else{
				$genotype = $key."_"."$conversion{$j} WT";
				push(@returns,$genotype);
			}
		}		
	} #loop over all ratios in row

	($call1, $call2) = @returns;
	$returnCall = $call1."|".$call2;
return $returnCall;
}# end sub call132

##____________________________________________________________________________________________________________________________________________________________##
#### sub prioritizeMutations accepts a reference to an array that contains all of the possible SNPs. A priority Hash and a reversed priority hash are generated
#### according to the instructions laid out by John Woody. Can be modified as needed. The array is looped over and the priority of the call is determined
#### by looking at the priorityHash and the priority is pushed into an array. The array is sorted and used as a guide to build a hash with the priority as key
#### and the value is the SNP. A reference to the hash is returned to the main block. 
sub prioritizeMutations($){
my ($ref) = @_;
my @array = @$ref;
## Variable declaration
my($query, $element, $counter, $pri);
my (@ordered, @sorted);
my (%hash);
my %priorityHash =( "5"  => "1",
					"14" => "2",
					"12" => "3",
					"8"  => "4",
					"7"  => "5",
					"6"  => "6",
					"3"  => "7",
					"4"  => "8",
					"41" => "9",
					"29" => "10",
					"17" => "11",
					"10" => "12",
					"9"  => "13",
					"2"  => "14"
					);
					
my %reversedHash = ( "1"  => "5",
					"2" => "14",
					"3" => "12",					
					"4"  => "8",
					"5"  => "7",
					"6"  => "6",
					"7"  => "3",
					"8"  => "4",
					"9" => "41",
					"10" => "29",
					"11" => "17",
					"12" => "10",
					"13"  => "9",
					"14"  => "2"
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
}# end prioritize mutations sub-routine 

##____________________________________________________________________________________________________________________________________________________________##
#### sub checkLinkage accepts a reference to an array that contains all SNP calls and checks for possible Linkages. A hash is created that contains the linkage
#### as a key and the value is the element that will not be called. Two nested loops create all possible ordered combinations and checks the linkage hash. If it
#### exists, the value that needs removed from the array is pushed into a hash called comboHash. The next step is to loop over the original array and remove 
#### the elements that are linked together. 
sub checkLinkage($){
my ($ref) =@_;
my @data =@$ref;

my($combo, $reverseCombo, $element, $element2, $pair);
my (@linkedCall, @unordered, @ordered);

my (%comboHash);
##linked Hash has as a value the element that needs removed!!! 
my %linked =(	"2_17" =>"2",
				"2_29" => "2",
				"2_41" => "2",
				"4_10"=> "10"
			);
		
foreach my $call (@data){	
	($element) = $call =~/(\d+)/;		
	foreach my $call2(@data){		
		($element2) = $call2 =~/(\d+)/;
		@unordered =($element,$element2);
		@ordered = sort {$a <=> $b} @unordered;
		$pair = $ordered[0]."_".$ordered[1];
		if(exists $linked{$pair}){
			$comboHash{$linked{$pair}}="";
		}
	}
}	
for (my $i=0; $i<scalar @data; $i++){
	my ($comparison) = $data[$i] =~/(\d+)/;
	if(exists $comboHash{$comparison}){
		splice(@data, $i,1);
	}
}
return \@data;
}# end sub checkLinkage


##### sub Format is a subroutine that orients the genotyped calls into the correct format. Example *2 *17 not *17 *2 The sub accepts two scalar values and sorts
##### them numerically and formats them into a genotype call and returns to the main block. 
sub Format($$){
my($call1, $call2) = @_;
my (@unordered, @ordered, $call);
push(@unordered, $call1);
push(@unordered, $call2);
@ordered = sort {$a <=> $b} @unordered;
$call = "*$ordered[0]|*$ordered[1]"; #<--Version 1:02
return $call;
}# end checkFormatting 

1; #return statement
