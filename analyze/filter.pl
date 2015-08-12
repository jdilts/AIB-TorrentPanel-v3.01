#!/usr/bin/perl -w


## Company: American International Biotechnology, LLC 
## Author: Jeri Dilts
## Date: April 2014


use strict;


#-----VARIABLES----------------------
my ($tsp,$sample) = @ARGV;

#-----PROTOTYPES----------------------
sub filter_bam($$$$);

#-----METRICS-------------------------
my %MappingQCmetrics = (

	APOEs1s2 => 'default',
	CYP2C19s17 => 'default',
	CYP2C19s2 => 'default',
	CYP2C19s3 => 'default',
	CYP2C19s4 => 'default',
	CYP2C19s7 => 'default',
	CYP2C19s6s8s9 => 'default',
	CYP2C9s2M => 'default',
	CYP2C9s6M => 'default',
	CYP2C9s5Ms3Ms11Ms4M	=> 'default',
	CYP2D6s10s12 => 'default',
	CYP2D6s17 => 'default',
	CYP2D6s2s7s41 => 'default',
	CYP2D6s8s14s6s29 => 'default',
	CYP2D6s3s9 => 'default',
	CYP2D6s4 => 'default', 
	VKORC1 => '18',
	CYP3A4s12 => 'default',
	CYP3A4s17 => 'default',
	CYP3A4s1B => 'default',
	CYP3A4s2 => 'default',
	CYP3A4s3 => 'default',
	CYP3A5s2 => 'default',
	CYP3A5s3 => 'default',
	CYP3A5s3Bs8 => 'default',
	CYP3A5s6 => 'default',
	CYP3A5s7 => 'default',
	CYP3A5s9 => 'default',
	FactorII => 'default',
	FVLeiden => 'default',
	MTHFRC677T => 'default',
	MTHFRA1298C => 'default',
	ATP13A4 => 'default',
	PALLD => 'default',
	ADAMTS2 => 'default',
	PTN => 'default',
	TRDMT1 => 'default',
	LEPREL2 => 'default',
	RAB31 => 'default',
	SYN3 => 'default',
	AmelX => 'default',
	AmelY => 'default',
    );
#-------------------------------------- *default=10? 


#-----EXECUTE-------------------------
my $run;
if($tsp=~/\/results\/analysis\/output\/Home\/(.*)\/plugin_out\//){$run = $1; chomp $run;}

my ($filtered_sam) = filter_bam($run,$sample,$tsp,\%MappingQCmetrics);

#-----SUBROUTINES----------------------
sub filter_bam($$$$){

	## variables
	my ($run,$sample,$MappingQCmetrics) = @_;
	my $bam = "/Home/$run/$sample"."_rawlib.bam";
	my $bai = $bam.".bai";
	my ($filtered_bam,$filtered_bai,$filtered_sam,$barcode);
		
	unless(($bam eq "/Home/$run/rawlib.bam")||($bam eq "/Home/$run/nomatch_rawlib.bam")){
		
		if($bam=~/(\/Home\/$run\/).*_rawlib.bam/){
				
			$filtered_bam = "$tsp/$sample/"."$sample"."_rawlib.bam";
			$filtered_bai = "$tsp/$sample/"."$sample"."_rawlib.bam.bai";
			$filtered_sam = "$tsp/$sample/"."$sample"."_rawlib.sam";
		}
		
		open(OUT, ">$filtered_sam");
		open(IN, "samtools view -h $bam |");
				
		while(<IN>){
				
			if( $_ =~ /^@/ ) {
				
				print OUT $_;
					
			} else {
				
				my @line = split( /\s+/ ); 
					
				chomp $line[2];
					
				if (exists $MappingQCmetrics{$line[2]}){
					
					if ($MappingQCmetrics{$line[2]} ne 'default'){
						
						my $filter_param = $MappingQCmetrics{$line[2]};			
							
						if ($line[4] >= $filter_param){print OUT $_;}
							
						else{print OUT $_;}
						
					} else {print OUT $_;}
				}
			}
		}
		
		close(OUT);
		close(IN);
		
		system("java -jar /opt/picard/picard-tools-current/SamFormatConverter.jar I=$filtered_sam O=$filtered_bam"); #convert to bam
		system("samtools index $filtered_bam $filtered_bai"); #index	
	}
}

























