#!/usr/bin/perl -w


## Author: Jeri Dilts
## Date: 
## Company: American International Biotechnology
## Department: Bioinformatics


use strict;
use File::Slurp;

## prototypes
sub create_table($$$);

## variables
my ($non_barcode_dropouts,$nt_freq,$tsp) = @ARGV;

## main
create_table($non_barcode_dropouts,$nt_freq,$tsp);

## subroutines
sub create_table($$$){

	## variables
	my $non_barcode_dropouts = shift;
	my $nt_freq = shift;
	my $tsp = shift;

	my @barcodes = read_file $non_barcode_dropouts; 
	
	open(OUT,">>$nt_freq") || die $!;
	
	foreach (@barcodes){
	
		unless(($_=~/^#/) || ($_=~/^$/)){
			my $sample = $_;
			chomp $sample;

			unless ($sample =~/.*blank.*/){ #just in case data shows up in blank, which does happen time to time

				my $allele_counts = "$tsp/$sample/allele_counts.xls";
	
				open(IN, "<$allele_counts") || die $!;
				print OUT "$sample\n";
				print OUT "Chrom\tPosition\tRef\tCov\tA Reads\tC Reads\tG Reads\tT Reads\tDeletions\n";
		
				while(<IN>){
	
					my $line = $_;
					unless($line=~/^Chrom/){
						my @line = split(/\t/,$line);
						print OUT "$line[0]\t$line[1]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\t$line[12]";
					}
				}
			}
		
		unless($sample=~/.*blank.*/){print OUT "\n";} #just in case data shows up in the blank, which does happen time to time
		close(IN);
		}
	}	
	close(OUT);
}
