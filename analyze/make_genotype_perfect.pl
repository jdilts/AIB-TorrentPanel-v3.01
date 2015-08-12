#####################################################
## Author: Snehal K. Talati                        ##
## Date: June 10th, 2014                           ##
## Company: American International Biotechnology   ##
## Department: Bioinformatics                      ##
## Version: Version 1                              ##
#####################################################


use strict;
use warnings;
use File::Slurp;

######################################################################################################################
my ($input_file,$tsp) = @ARGV;
my $PerfectFile = $tsp."/AIB_RESULTS/genotype.txt";
my $ErrorFile = $tsp."/AIB_RESULTS/genotype_errors.txt"; 
my %MainHash;
my @Info2;

#################################################################################################################################

#################################################################################################################################
sub buildHash ($);
##################################################### MAIN #######################################################################

my ($hashIndex,$header) = buildHash($input_file);
%MainHash = %$hashIndex;
open (OFH, ">$PerfectFile") or die "Cannot open output file$!";
open (OFH2, ">$ErrorFile") or die "Cannot open output file$!";

## print header
my @header = @$header;
foreach(@header){print OFH "$_";}
print OFH "\n";

foreach my $key (keys %MainHash) {
	my $substring = "Error";
	if (not ($key =~ /^POS/)) {
		if ($MainHash{$key} =~ /Error.*|No Call.*/ig) {
			@Info2 = split ('\t', $MainHash{$key});	
			print OFH2 "$key", "\n", join("\n", @Info2), "\n\n\n";	
				
		}
		else {
			@Info2 = split ('\t', $MainHash{$key});
			print OFH "$key", "\n", join("\n", @Info2), "\n\n\n";	
		}
	}
}

close OFH;
close OFH2;



##################################################### SUBROUTINES ################################################################
sub buildHash ($) {
my ($input) = @_;
my ($barcode, $flag1, $flag2, $line, $Access, $newline);
$flag1 = 0;
$flag2 = 0;
my $flag = 0;
my @temp_sample;
my @contents = read_file($input);
my %Hash;
my @header;

foreach (@contents) {

	if($_=~/^#/){
		push(@header,$_);
	}


	unless(($_=~/^#/) || ($_=~/^$/) && ($flag1==0)){ #used to pick up first sample

                $flag1=1;
                chomp $_;
				#print "$_\n";
				
				if ($_ =~/^$/){$flag2 = 1;} #finds blank lines
                else {$flag2 = 0;}
				
				push (@temp_sample, $_) unless ($flag2 == 1);
				
				  if ($flag2 == 1){

                        $barcode = $temp_sample[0];
                        shift @temp_sample;

                        
						$line = "";
                        foreach(@temp_sample){
                                
								my $substring = "Error";
								if (index($_, $substring) != -1) {
									if (($_ ne "")&&($flag==1)){
										$line .= "$_\t";
										#print "$_\n";
										}
									if (($_ ne "")&&($flag==0)){
										$line .= "$_\t";
										#print "$_\n";
										$flag=1;
										}
								}
								else {
									if (($_ ne "")&&($flag==1)){
										$line .= "$_\t";
										#print "$_\n";
										}
									if (($_ ne "")&&($flag==0)){
										$line .= "$_\t";
										#print "$_\n";
										$flag=1;
										}
								}
							
                        }
							
							$Hash{$barcode}= $line;

                        # clear
                        @temp_sample = (); #clear patient data
                        $flag2 = 0; # and set flag back to zero
                }
				
				
}
}
return (\%Hash,\@header);
}














