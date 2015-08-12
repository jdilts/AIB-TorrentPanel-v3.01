#!/usr/bin/perl -w

my ($run,$barcode_dropouts_path,$non_barcode_dropouts_path) = @ARGV;

use strict;
use File::Slurp;
use List::Compare;

#----- PROTOTYPES ---------------------------------
sub grabBams($);
sub grabBarcodeList($);
sub createOrderBarcodeHash($);
sub compareIDs($$$);
sub printToFLAT($$$);

#----- EXECUTE -----------------------------------
my ($bamID_r) = grabBams($run);
my ($barcodeID_r) = grabBarcodeList($run);
my ($order_r ) = createOrderBarcodeHash($barcodeID_r);
my ($intersection_r) = compareIDs($bamID_r,$barcodeID_r,$barcode_dropouts_path);
printToFLAT($intersection_r,$non_barcode_dropouts_path,$order_r);


#////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////
#---- SUBROUTINES-------------------------------------
#/////////////////////////////////////////////////////

sub grabBams($){

	my $run = shift;

	#grab present bams in run dir
	my @bam_id;
	my @bams = glob("/results/analysis/output/Home/$run/*bam");

	foreach my $file (@bams){
		if ($file=~/\/results\/analysis\/output\/Home\/$run\/(.*)_rawlib\.bam/){push(@bam_id, $1);}
	}

    return(\@bam_id);
}

sub grabBarcodeList($){

	my $run = shift;

	#grab barcodeList.txt
	my @barcode_id;
	my @barcodeList = read_file("/results/analysis/output/Home/$run/barcodeList.txt");
	chomp @barcodeList;
	splice(@barcodeList,0,3);

	foreach my $line (@barcodeList){

		my @p = split(',',$line);
		push(@barcode_id, $p[1]);
	}

	return(\@barcode_id);
}

sub createOrderBarcodeHash($){

	my $barcode_id = shift;
	my @barcode_id = @$barcode_id;
	my $i = 1;
	my %order;

	foreach(@barcode_id){
	
		$order{$i} = $_;
		$i++;
	}
	
	return(\%order);
}

sub compareIDs($$$){

	my ($bamID_r,$barcodeID_r,$barcode_dropouts_path) = @_;
	my @bam_id = @$bamID_r;
	my @barcode_id = @$barcodeID_r;

	#compare array contents
	my $lc = List::Compare->new(\@bam_id, \@barcode_id);
	my @diff = $lc->get_symmetric_difference;
	my @intersection = $lc->get_intersection;
	
	open (DIFF,">>",$barcode_dropouts_path) || die $!;
	foreach(@diff){chomp $_; print DIFF "$_\n";}
	close(DIFF);
		
	return(\@intersection);
}

sub printToFLAT($$$){

	my ($intersection,$non_barcode_dropouts_path,$order) = @_;
	my @intersection = @$intersection; #deref
	my %order = %$order; #deref
	
	open(DATA,">>",$non_barcode_dropouts_path) || die $!;

	my $i = 1;

	my $hash_size = keys(%order);
	
	until($i > $hash_size){

		foreach(@intersection){

			if ($order{$i} eq $_){
				print DATA "$order{$i}\n";
			}
		}
		$i++;
	}
	close(DATA);
}
