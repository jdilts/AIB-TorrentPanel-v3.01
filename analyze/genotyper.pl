#!/usr/bin/perl -w


## Author: Jeri Dilts
## Date: April 23rd, 2014
## Company: American International Biotechnology
## Department: Bioinformatics


my ($plugin_name,$plugin_filepath,$sample_filter,$snp_coverage,$ntfreq,$genotype,$sample_cov) = @ARGV;


use strict;
use lib ("/results/plugins/$ARGV[0]/analyze/modules");
use File::Slurp; #reads a file into an array
use APOE qw(callAPOE);
use CYP2C19 qw(callCYP2C19);
use CYP2C9 qw(callCYP2C9);
use CYP2D6 qw(callCYP2D6);
use VKORC1 qw(callVKORC1);
use FactorII qw(callFactorII);
use FactorV qw(callFactorV);
use MTHFR qw(callMTHFR callMTHFR677 callMTHFR1298);
use CYP3A4 qw(callCYP3A4);
use CYP3A5 qw(callCYP3A5);
use SAMPLE_ID qw(report_sampleID);
use printTorrent qw(printAccession printCalls);

## prototypes
sub grab_sample_cov($);
sub genotype_samples($$$$$$);

## main
my $coverages = grab_sample_cov($sample_cov);
open(my $OUT, ">>$genotype");
genotype_samples($plugin_filepath,$sample_filter,$snp_coverage,$ntfreq,$coverages,$OUT);
close($OUT);

## subroutines
sub grab_sample_cov($){

	my $sample_cov = shift;
	my @samplecov = read_file($sample_cov); #sample coverage data
	my %coverages;
	
	foreach(@samplecov){

		my @line = split(/,/,$_);
		$coverages{$line[0]} = $line[1];
	}
	return(\%coverages);
}

sub genotype_samples($$$$$$){

	## variables
	my ($plugin_filepath,$sample_filter,$snp_coverage,$ntfreq,$coverages,$OUT) = @_;
	my @ntfreq = read_file($ntfreq);
	my %coverages = %$coverages;
	my @temp_sample;
	my $flag = 0;
	my $flag2 = 0;
	
	foreach(@ntfreq){

	unless ((($_=~/^#/) || ($_=~/^$/)) && ($flag2 ==0)) {

		$flag2 = 1;
		chomp $_;

		if ($_ =~/^$/){$flag = 1;} #finds blank lines
		else {$flag = 0;}
	
		push (@temp_sample, $_) unless ($flag == 1);
		
		if ($flag == 1){
			
			my $barcode = $temp_sample[0];

			unless($barcode eq ""){

				printAccession($barcode, $OUT);
			
				shift @temp_sample;
				shift @temp_sample;
			
				# places the amplicon and each position as a key in hash %temp_sample
				my %SNP;
				foreach(@temp_sample){
			
					my $line = $_;
					my @splitLine = split(/\t/,$line);
					my $ampliconID = $splitLine[0];
					my $ampliconPosition = $splitLine[1];		
					my $key = $ampliconID.",".$ampliconPosition;
				
					my @nucCalls;
					my $value;
					for(my $j=2; $j< scalar @splitLine; $j++){
						my $temp = $splitLine[$j];
						push(@nucCalls, $temp);
						$value = \@nucCalls;
					}
				
					$SNP{$key} = $value;
				}

				my $reads = $coverages{$barcode};
				
				#////////////////////////////////////////
				# EXECUTE GENE MODULES
				#////////////////////////////////////////
				my $APOE = callAPOE(\%SNP,$reads,$snp_coverage);
				my $C19 = callCYP2C19(\%SNP,$reads,$snp_coverage);
				my $C9 = callCYP2C9(\%SNP,$reads,$snp_coverage);
				my $D6 = callCYP2D6(\%SNP,$reads,$snp_coverage);
				my $VKORC1 = callVKORC1(\%SNP,$reads,$snp_coverage);
				my $factorII = callFactorII(\%SNP,$reads,$snp_coverage);
				my $factorV = callFactorV(\%SNP,$reads,$snp_coverage);
				my $MTHFR = callMTHFR(\%SNP,$reads,$snp_coverage);
				my $MTHFR677 = callMTHFR677(\%SNP,$reads,$snp_coverage);
				my $MTHFR1298 = callMTHFR1298(\%SNP,$reads,$snp_coverage);
				my $CYP3A4 = callCYP3A4(\%SNP,$reads,$snp_coverage);
				my $CYP3A5 = callCYP3A5(\%SNP,$reads,$snp_coverage);
				my $sample_id = report_sampleID($sample_filter,\%SNP,$reads,$snp_coverage,$plugin_filepath,$barcode);
							
				my @calls = ($APOE,$C19,$C9,$D6,$VKORC1,$factorII,$factorV,$MTHFR,$MTHFR677,$MTHFR1298,$CYP3A4,$CYP3A5);
				printCalls(\@calls, $OUT);
				print $OUT "$sample_id";
				print $OUT "\n";

				@temp_sample = (); #clear patient data
				$flag = 0; # and set flag back to zero
			}
		}
	}
	}
}


