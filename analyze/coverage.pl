#!/usr/bin/perl -w


## Author: Jeri Dilts
## Date: 
## Company: American International Biotechnology
## Department: Bioinformatics
## 


use strict;
use File::Slurp;


## prototypes
sub calc_covPERsample($$$);

## variables
my ($non_dropout_barcodes,$sample_coverage,$tsp) = @ARGV;

## main
calc_covPERsample($non_dropout_barcodes,$sample_coverage,$tsp);

## subroutines
sub calc_covPERsample($$$){

	## variables
	my $non_dropout_barcodes = shift;
	my $sample_coverage = shift;
	my $tsp = shift;

	my @barcodes = read_file $non_dropout_barcodes; 
	
	my %SAMPLE_ID = (
	
		ATP13A4 => 'default',
		PALLD => 'default',
		ADAMTS2 => 'default',
		PTN => 'default',
		TRDMT1 => 'default',
		LEPREL2 => 'default',
		RAB31 => 'default',
		SYN3 => 'default',
		AmelX => 'default',
		AmelY => 'default'
	);

	open(OUT,">>$sample_coverage") || die $!;
	
	foreach (@barcodes){
	
	unless(($_=~/^#/)||($_=~/^$/)){
		my $sample = $_;
		chomp $sample;
		my $filtered_sam = "$tsp/$sample/$sample"."_rawlib.sam";
	
		open(IN, "<$filtered_sam");
		my $i = 0;
		while(<IN>){
	
			my $line = $_;
		
			unless($line=~/^@/ ){
		
				my @line = split(/\s+/, $line);
				chomp $line[2];
				if(exists $SAMPLE_ID{$line[2]}){$i++;}
			}
		}
		close(IN);
	
		

		
		
		print OUT "$sample,$i\n";
		
	}

	}	
	close(OUT);
}
