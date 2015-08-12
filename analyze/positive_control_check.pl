#!/usr/bin/perl -w

## Author: Jeri Dilts 
## Date: May 8th, 2014 
## Company: American International 
## Biotechnology Department: Bioinformatics

use strict;
use File::Slurp;
use Switch; #core module

my ($genotype_file,$output) = @ARGV;
my @f = read_file($genotype_file); 
for(my $i = 1; $i < 11; $i++){shift(@f);} #shifts off header
my $pos_control_id = $f[0]; #capture positive control id
chomp $pos_control_id;
shift(@f);

## open out file
open OUT, ">>$output" || die $!;
print OUT "$pos_control_id\n";

for(my $i = 0; $i < 10; $i++){

	my $gene = "";
	my $error = "";
	my $hap1 =  "";
	my $hap2 = "";
	
	my $line = $f[$i];
	chomp $line;
	my @line = split(/\|/,$line);
	if(defined $line[0]){$gene = $line[0];}
	if(defined $line[1]){$error = $line[1];}
	if(defined $line[2]){$hap1 = $line[2];}
	if(defined $line[3]){$hap2 = $line[3];}

	#how to handle errors per gene

	switch($pos_control_id){

		case 'POS-NA16000' {
			switch($gene){
				case "APOE"{if (($hap1 eq "E3") && ($hap2 eq "E3")){print OUT "APOE looks good.\n";}else{print OUT "APOE is INCORRECT.\n";}}
				case "CYP2C19"{if (($hap1 eq "*17") && ($hap2 eq "*17")){print OUT "CYP2C19 looks good.\n";}else{print OUT "CYP2C19 is INCORRECT.\n";}}
				case "CYP2C9"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP2C9 looks good.\n";}else{print OUT "CYP2C9 is INCORRECT.\n";}}
				case "CYP2D6"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP2D6 looks good.\n";}else{print OUT "CYP2D6 is INCORRECT.\n";}}
				case "VKORC1"{if (($hap1 eq "*A") && ($hap2 eq "*B")){print OUT "VKORC1 looks good.\n";}else{print OUT "VKORC1 is INCORRECT.\n";}}
				case "FactorII"{if ($hap1 eq "Detected:Homozygotemutant"){print OUT "FactorII looks good.\n";}else{print OUT "FactorII is INCORRECT.\n";}}
				case "FVLeiden"{if ($hap1 eq "NotDetected"){print OUT "FVLeiden looks good.\n";}else{print OUT "FVLeiden is INCORRECT.\n";}}
				case "MTHFR"{if (($hap1 eq "HET") && ($hap2 eq "HET")){print OUT "MTHFR looks good.\n";}else{print OUT "MTHFR is INCORRECT.\n";}}
				case "CYP3A4"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP3A4 looks good.\n";}else{print OUT "CYP3A4 is INCORRECT.\n";}}
				case "CYP3A5"{if (($hap1 eq "*3") && ($hap2 eq "*3")){print OUT "CYP3A5 looks good.\n";}else{print OUT "CYP3A5 is INCORRECT.\n";}}
				#case "Sample ID"{if ($hap1 eq "E3"){next;}}
				#else{print OUT "$pos_control_id is INCORRECT.\n";}
			}
		}
		case 'POS-NA16689'{
			switch($gene){
				case "APOE"{if (($hap1 eq "E3") && ($hap2 eq "E3")){print OUT "APOE looks good.\n";}else{print OUT "APOE is INCORRECT.\n";}}
				case "CYP2C19"{if (($hap1 eq "*2") && ($hap2 eq "*2")){print OUT "CYP2C19 looks good.\n";}else{print OUT "CYP2C19 is INCORRECT.\n";}}
				case "CYP2C9"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP2C9 looks good.\n";}else{print OUT "CYP2C9 is INCORRECT.\n";}}
				case "CYP2D6"{if (($hap1 eq "*2") && ($hap2 eq "*10")){print OUT "CYP2D6 looks good.\n";}else{print OUT "CYP2D6 is INCORRECT.\n";}}
				case "VKORC1"{if (($hap1 eq "*A") && ($hap2 eq "*A")){print OUT "VKORC1 looks good.\n";}else{print OUT "VKORC1 is INCORRECT.\n";}}
				case "FactorII"{if ($hap1 eq "NotDetected"){print OUT "FactorII looks good.\n";}else{print OUT "FactorII is INCORRECT.\n";}}
				case "FVLeiden"{if ($hap1 eq "NotDetected"){print OUT "FVLeiden looks good.\n";}else{print OUT "FVLeiden is INCORRECT.\n";}}
				case "MTHFR"{if (($hap1 eq "WT") && ($hap2 eq "HET")){print OUT "MTHFR looks good.\n";}else{print OUT "MTHFR is INCORRECT.\n";}}
				case "CYP3A4"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP3A4 looks good.\n";}else{print OUT "CYP3A4 is INCORRECT.\n";}}
				case "CYP3A5"{if (($hap1 eq "*3") && ($hap2 eq "*3")){print OUT "CYP3A5 looks good.\n";}else{print OUT "CYP3A5 is INCORRECT.\n";}}
				#case "Sample ID"{if ($hap1 eq "E3"){next;}}
				#else{print OUT "$pos_control_id is INCORRECT.\n";}
			}
		}
		case 'POS-NA17052'{
			switch($gene){
				case "APOE"{if (($hap1 eq "E3") && ($hap2 eq "E3")){print OUT "APOE looks good.\n";}else{print OUT "APOE is INCORRECT.\n";}}
				case "CYP2C19"{if (($hap1 eq "*1") && ($hap2 eq "*3")){print OUT "CYP2C19 looks good.\n";}else{print OUT "CYP2C19 is INCORRECT.\n";}}
				case "CYP2C9"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP2C9 looks good.\n";}else{print OUT "CYP2C9 is INCORRECT.\n";}}
				case "CYP2D6"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP2D6 looks good.\n";}else{print OUT "CYP2D6 is INCORRECT.\n";}}
				case "VKORC1"{if (($hap1 eq "*A") && ($hap2 eq "*A")){print OUT "VKORC1 looks good.\n";}else{print OUT "VKORC1 is INCORRECT.\n";}}
				case "FactorII"{if ($hap1 eq "NotDetected"){print OUT "FactorII looks good.\n";}else{print OUT "FactorII is INCORRECT.\n";}}
				case "FVLeiden"{if ($hap1 eq "NotDetected"){print OUT "FVLeiden looks good.\n";}else{print OUT "FVLeiden is INCORRECT.\n";}}
				case "MTHFR"{if (($hap1 eq "HET") && ($hap2 eq "HET")){print OUT "MTHFR looks good.\n";}else{print OUT "MTHFR is INCORRECT.\n";}}
				case "CYP3A4"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP3A4 looks good.\n";}else{print OUT "CYP3A4 is INCORRECT.\n";}}
				case "CYP3A5"{if (($hap1 eq "*1") && ($hap2 eq "*3")){print OUT "CYP3A5 looks good.\n";}else{print OUT "CYP3A5 is INCORRECT.\n";}}
				#case "Sample ID"{if ($hap1 eq "E3"){next;}}
				#else{print OUT "$pos_control_id is INCORRECT.\n";}
			}
		}
		case 'POS-NA17058'{
			switch($gene){
				case "APOE"{if (($hap1 eq "E3") && ($hap2 eq "E3")){print OUT "APOE looks good.\n";}else{print OUT "APOE is INCORRECT.\n";}}
				case "CYP2C19"{if (($hap1 eq "*1") && ($hap2 eq "*2")){print OUT "CYP2C19 looks good.\n";}else{print OUT "CYP2C19 is INCORRECT.\n";}}
				case "CYP2C9"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP2C9 looks good.\n";}else{print OUT "CYP2C9 is INCORRECT.\n";}}
				case "CYP2D6"{if (($hap1 eq "*10") && ($hap2 eq "*10")){print OUT "CYP2D6 looks good.\n";}else{print OUT "CYP2D6 is INCORRECT.\n";}}
				case "VKORC1"{if (($hap1 eq "*A") && ($hap2 eq "*A")){print OUT "VKORC1 looks good.\n";}else{print OUT "VKORC1 is INCORRECT.\n";}}
				case "FactorII"{if ($hap1 eq "NotDetected"){print OUT "FactorII looks good.\n";}else{print OUT "FactorII is INCORRECT.\n";}}
				case "FVLeiden"{if ($hap1 eq "NotDetected"){print OUT "FVLeiden looks good.\n";}else{print OUT "FVLeiden is INCORRECT.\n";}}
				case "MTHFR"{if (($hap1 eq "HET") && ($hap2 eq "HET")){print OUT "MTHFR looks good.\n";}else{print OUT "MTHFR is INCORRECT.\n";}}
				case "CYP3A4"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP3A4 looks good.\n";}else{print OUT "CYP3A4 is INCORRECT.\n";}}
				case "CYP3A5"{if (($hap1 eq "*1") && ($hap2 eq "*3")){print OUT "CYP3A5 looks good.\n";}else{print OUT "CYP3A5 is INCORRECT.\n";}}
				#case "Sample ID"{if ($hap1 eq "E3"){next;}}
				#else{print OUT "$pos_control_id is INCORRECT.\n";}
			}
		}
		case 'POS-NA17114'{
			switch($gene){
				case "APOE"{if (($hap1 eq "E3") && ($hap2 eq "E4")){print OUT "APOE looks good.\n";}else{print OUT "APOE is INCORRECT.\n";}}
				case "CYP2C19"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP2C19 looks good.\n";}else{print OUT "CYP2C19 is INCORRECT.\n";}}
				case "CYP2C9"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP2C9 looks good.\n";}else{print OUT "CYP2C9 is INCORRECT.\n";}}
				case "CYP2D6"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP2D6 looks good.\n";}else{print OUT "CYP2D6 is INCORRECT.\n";}}
				case "VKORC1"{if (($hap1 eq "*B") && ($hap2 eq "*B")){print OUT "VKORC1 looks good.\n";}else{print OUT "VKORC1 is INCORRECT.\n";}}
				case "FactorII"{if ($hap1 eq "NotDetected"){print OUT "FactorII looks good.\n";}else{print OUT "FactorII is INCORRECT.\n";}}
				case "FVLeiden"{if ($hap1 eq "NotDetected"){print OUT "FVLeiden looks good.\n";}else{print OUT "FVLeiden is INCORRECT.\n";}}
				case "MTHFR"{if (($hap1 eq "WT") && ($hap2 eq "HET")){print OUT "MTHFR looks good.\n";}else{print OUT "MTHFR is INCORRECT.\n";}}
				case "CYP3A4"{if (($hap1 eq "*1") && ($hap2 eq "*1B")){print OUT "CYP3A4 looks good.\n";}else{print OUT "CYP3A4 is INCORRECT.\n";}}
				case "CYP3A5"{if (($hap1 eq "*1") && ($hap2 eq "*3")){print OUT "CYP3A5 looks good.\n";}else{print OUT "CYP3A5 is INCORRECT.\n";}}
				#case "Sample ID"{if ($hap1 eq "E3"){next;}}
				#else{print OUT "$pos_control_id is INCORRECT.\n";}
			}
		}
		case 'POS-NA17115'{
			switch($gene){
				case "APOE"{if (($hap1 eq "E3") && ($hap2 eq "E4")){print OUT "APOE looks good.\n";}else{print OUT "APOE is INCORRECT.\n";}}
				case "CYP2C19"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP2C19 looks good.\n";}else{print OUT "CYP2C19 is INCORRECT.\n";}}
				case "CYP2C9"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP2C9 looks good.\n";}else{print OUT "CYP2C9 is INCORRECT.\n";}}
				case "CYP2D6"{if (($hap1 eq "*1") && ($hap2 eq "*2")){print OUT "CYP2D6 looks good.\n";}else{print OUT "CYP2D6 is INCORRECT.\n";}}
				case "VKORC1"{if (($hap1 eq "*B") && ($hap2 eq "*B")){print OUT "VKORC1 looks good.\n";}else{print OUT "VKORC1 is INCORRECT.\n";}}
				case "FactorII"{if ($hap1 eq "NotDetected"){print OUT "FactorII looks good.\n";}else{print OUT "FactorII is INCORRECT.\n";}}
				case "FVLeiden"{if ($hap1 eq "NotDetected"){print OUT "FVLeiden looks good.\n";}else{print OUT "FVLeiden is INCORRECT.\n";}}
				case "MTHFR"{if (($hap1 eq "WT") && ($hap2 eq "WT")){print OUT "MTHFR looks good.\n";}else{print OUT "MTHFR is INCORRECT.\n";}}
				case "CYP3A4"{if (($hap1 eq "*1") && ($hap2 eq "*1B")){print OUT "CYP3A4 looks good.\n";}else{print OUT "CYP3A4 is INCORRECT.\n";}}
				case "CYP3A5"{if (($hap1 eq "*1") && ($hap2 eq "*3")){print OUT "CYP3A5 looks good.\n";}else{print OUT "CYP3A5 is INCORRECT.\n";}}
				#case "Sample ID"{if ($hap1 eq "E3"){next;}}
				#else{print OUT "$pos_control_id is INCORRECT.\n";}
			}
		}
		case 'POS-NA17129' {
			switch($gene){
				case "APOE"{if (($hap1 eq "E3") && ($hap2 eq "E3")){print OUT "APOE looks good.\n";}else{print OUT "APOE is INCORRECT.\n";}}
				case "CYP2C19"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP2C19 looks good.\n";}else{print OUT "CYP2C19 is INCORRECT.\n";}}
				case "CYP2C9"{if (($hap1 eq "*1") && ($hap2 eq "*2")){print OUT "CYP2C9 looks good.\n";}else{print OUT "CYP2C9 is INCORRECT.\n";}}
				case "CYP2D6"{if (($hap1 eq "*1") && ($hap2 eq "*4")){print OUT "CYP2D6 looks good.\n";}else{print OUT "CYP2D6 is INCORRECT.\n";}}
				case "VKORC1"{if (($hap1 eq "*B") && ($hap2 eq "*B")){print OUT "VKORC1 looks good.\n";}else{print OUT "VKORC1 is INCORRECT.\n";}}
				case "FactorII"{if ($hap1 eq "NotDetected"){print OUT "FactorII looks good.\n";}else{print OUT "FactorII is INCORRECT.\n";}}
				case "FVLeiden"{if ($hap1 eq "NotDetected"){print OUT "FVLeiden looks good.\n";}else{print OUT "FVLeiden is INCORRECT.\n";}}
				case "MTHFR"{if (($hap1 eq "WT") && ($hap2 eq "HET")){print OUT "MTHFR looks good.\n";}else{print OUT "MTHFR is INCORRECT.\n";}}
				case "CYP3A4"{if (($hap1 eq "*1") && ($hap2 eq "*1B")){print OUT "CYP3A4 looks good.\n";}else{print OUT "CYP3A4 is INCORRECT.\n";}}
				case "CYP3A5"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP3A5 looks good.\n";}else{print OUT "CYP3A5 is INCORRECT.\n";}}
				#case "Sample ID"{if ($hap1 eq "E3"){next;}}
				#else{print OUT "$pos_control_id is INCORRECT.\n";}
			}
		}
		case 'POS-NA17130'{
			switch($gene){
				case "APOE"{if (($hap1 eq "E2") && ($hap2 eq "E3")){print OUT "APOE looks good.\n";}else{print OUT "APOE is INCORRECT.\n";}}
				case "CYP2C19"{if (($hap1 eq "*1") && ($hap2 eq "*17")){print OUT "CYP2C19 looks good.\n";}else{print OUT "CYP2C19 is INCORRECT.\n";}}
				case "CYP2C9"{if (($hap1 eq "*1") && ($hap2 eq "*3")){print OUT "CYP2C9 looks good.\n";}else{print OUT "CYP2C9 is INCORRECT.\n";}}
				case "CYP2D6"{if (($hap1 eq "*1") && ($hap2 eq "*2")){print OUT "CYP2D6 looks good.\n";}else{print OUT "CYP2D6 is INCORRECT.\n";}}
				case "VKORC1"{if (($hap1 eq "*B") && ($hap2 eq "*B")){print OUT "VKORC1 looks good.\n";}else{print OUT "VKORC1 is INCORRECT.\n";}}
				case "FactorII"{if ($hap1 eq "NotDetected"){print OUT "FactorII looks good.\n";}else{print OUT "FactorII is INCORRECT.\n";}}
				case "FVLeiden"{if ($hap1 eq "NotDetected"){print OUT "FVLeiden looks good.\n";}else{print OUT "FVLeiden is INCORRECT.\n";}}
				case "MTHFR"{if (($hap1 eq "WT") && ($hap2 eq "WT")){print OUT "MTHFR looks good.\n";}else{print OUT "MTHFR is INCORRECT.\n";}}
				case "CYP3A4"{if (($hap1 eq "*1B") && ($hap2 eq "*1B")){print OUT "CYP3A4 looks good.\n";}else{print OUT "CYP3A4 is INCORRECT.\n";}}
				case "CYP3A5"{if (($hap1 eq "*1") && ($hap2 eq "*3")){print OUT "CYP3A5 looks good.\n";}else{print OUT "CYP3A5 is INCORRECT.\n";}}
				#case "Sample ID"{if ($hap1 eq "E3"){next;}}
				#else{print OUT "$pos_control_id is INCORRECT.\n";}
			}
		}
		case 'POS-NA17210'{
			switch($gene){
				case "APOE"{if (($hap1 eq "E3") && ($hap2 eq "E3")){print OUT "APOE looks good.\n";}else{print OUT "APOE is INCORRECT.\n";}}
				case "CYP2C19"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP2C19 looks good.\n";}else{print OUT "CYP2C19 is INCORRECT.\n";}}
				case "CYP2C9"{if (($hap1 eq "*1") && ($hap2 eq "*2")){print OUT "CYP2C9 looks good.\n";}else{print OUT "CYP2C9 is INCORRECT.\n";}}
				case "CYP2D6"{if (($hap1 eq "*1") && ($hap2 eq "*4")){print OUT "CYP2D6 looks good.\n";}else{print OUT "CYP2D6 is INCORRECT.\n";}}
				case "VKORC1"{if (($hap1 eq "*A") && ($hap2 eq "*A")){print OUT "VKORC1 looks good.\n";}else{print OUT "VKORC1 is INCORRECT.\n";}}
				case "FactorII"{if ($hap1 eq "NotDetected"){print OUT "FactorII looks good.\n";}else{print OUT "FactorII is INCORRECT.\n";}}
				case "FVLeiden"{if ($hap1 eq "NotDetected"){print OUT "FVLeiden looks good.\n";}else{print OUT "FVLeiden is INCORRECT.\n";}}
				case "MTHFR"{if (($hap1 eq "WT") && ($hap2 eq "Mut")){print OUT "MTHFR looks good.\n";}else{print OUT "MTHFR is INCORRECT.\n";}}
				case "CYP3A4"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP3A4 looks good.\n";}else{print OUT "CYP3A4 is INCORRECT.\n";}}
				case "CYP3A5"{if (($hap1 eq "*3") && ($hap2 eq "*3")){print OUT "CYP3A5 looks good.\n";}else{print OUT "CYP3A5 is INCORRECT.\n";}}
				#case "Sample ID"{if ($hap1 eq "E3"){next;}}
				#else{print OUT "$pos_control_id is INCORRECT.\n";}
			}
		}
		case 'POS-NA17221'{	
			switch($gene){
				case "APOE"{if (($hap1 eq "E3") && ($hap2 eq "E4")){print OUT "APOE looks good.\n";}else{print OUT "APOE is INCORRECT.\n";}}
				case "CYP2C19"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP2C19 looks good.\n";}else{print OUT "CYP2C19 is INCORRECT.\n";}}
				case "CYP2C9"{if (($hap1 eq "*2") && ($hap2 eq "*3")){print OUT "CYP2C9 looks good.\n";}else{print OUT "CYP2C9 is INCORRECT.\n";}}
				case "CYP2D6"{if (($hap1 eq "*1") && ($hap2 eq "*2")){print OUT "CYP2D6 looks good.\n";}else{print OUT "CYP2D6 is INCORRECT.\n";}}
				case "VKORC1"{if (($hap1 eq "*A") && ($hap2 eq "*B")){print OUT "VKORC1 looks good.\n";}else{print OUT "VKORC1 is INCORRECT.\n";}}
				case "FactorII"{if ($hap1 eq "NotDetected"){print OUT "FactorII looks good.\n";}else{print OUT "FactorII is INCORRECT.\n";}}
				case "FVLeiden"{if ($hap1 eq "NotDetected"){print OUT "FVLeiden looks good.\n";}else{print OUT "FVLeiden is INCORRECT.\n";}}
				case "MTHFR"{if (($hap1 eq "HET") && ($hap2 eq "WT")){print OUT "MTHFR looks good.\n";}else{print OUT "MTHFR is INCORRECT.\n";}}
				case "CYP3A4"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP3A4 looks good.\n";}else{print OUT "CYP3A4 is INCORRECT.\n";}}
				case "CYP3A5"{if (($hap1 eq "*3") && ($hap2 eq "*3")){print OUT "CYP3A5 looks good.\n";}else{print OUT "CYP3A5 is INCORRECT.\n";}}
				#case "Sample ID"{if ($hap1 eq "E3"){next;}}
				#else{print OUT "$pos_control_id is INCORRECT.\n";}
			}
		}
		case 'POS-NA17226'{
			switch($gene){
				case "APOE"{if (($hap1 eq "E2") && ($hap2 eq "E4")){print OUT "APOE looks good.\n";}else{print OUT "APOE is INCORRECT.\n";}}
				case "CYP2C19"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP2C19 looks good.\n";}else{print OUT "CYP2C19 is INCORRECT.\n";}}
				case "CYP2C9"{if (($hap1 eq "*1") && ($hap2 eq "*2")){print OUT "CYP2C9 looks good.\n";}else{print OUT "CYP2C9 is INCORRECT.\n";}}
				case "CYP2D6"{if (($hap1 eq "*4") && ($hap2 eq "*4")){print OUT "CYP2D6 looks good.\n";}else{print OUT "CYP2D6 is INCORRECT.\n";}}
				case "VKORC1"{if (($hap1 eq "*A") && ($hap2 eq "*B")){print OUT "VKORC1 looks good.\n";}else{print OUT "VKORC1 is INCORRECT.\n";}}
				case "FactorII"{if ($hap1 eq "NotDetected"){print OUT "FactorII looks good.\n";}else{print OUT "FactorII is INCORRECT.\n";}}
				case "FVLeiden"{if ($hap1 eq "NotDetected"){print OUT "FVLeiden looks good.\n";}else{print OUT "FVLeiden is INCORRECT.\n";}}
				case "MTHFR"{if (($hap1 eq "WT") && ($hap2 eq "WT")){print OUT "MTHFR looks good.\n";}else{print OUT "MTHFR is INCORRECT.\n";}}
				case "CYP3A4"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP3A4 looks good.\n";}else{print OUT "CYP3A4 is INCORRECT.\n";}}
				case "CYP3A5"{if (($hap1 eq "*3") && ($hap2 eq "*3")){print OUT "CYP3A5 looks good.\n";}else{print OUT "CYP3A5 is INCORRECT.\n";}}
				#case "Sample ID"{if ($hap1 eq "E3"){next;}}
				#else{print OUT "$pos_control_id is INCORRECT.\n";}
			}
		}
		case 'POS-NA17246'{
			switch($gene){
				case "APOE"{if (($hap1 eq "E3") && ($hap2 eq "E3")){print OUT "APOE looks good.\n";}else{print OUT "APOE is INCORRECT.\n";}}
				case "CYP2C19"{if (($hap1 eq "*8") && ($hap2 eq "*17")){print OUT "CYP2C19 looks good.\n";}else{print OUT "CYP2C19 is INCORRECT.\n";}}
				case "CYP2C9"{if (($hap1 eq "*1") && ($hap2 eq "*2")){print OUT "CYP2C9 looks good.\n";}else{print OUT "CYP2C9 is INCORRECT.\n";}}
				case "CYP2D6"{if (($hap1 eq "*2") && ($hap2 eq "*4")){print OUT "CYP2D6 looks good.\n";}else{print OUT "CYP2D6 is INCORRECT.\n";}}
				case "VKORC1"{if (($hap1 eq "*A") && ($hap2 eq "*B")){print OUT "VKORC1 looks good.\n";}else{print OUT "VKORC1 is INCORRECT.\n";}}
				case "FactorII"{if ($hap1 eq "NotDetected"){print OUT "FactorII looks good.\n";}else{print OUT "FactorII is INCORRECT.\n";}}
				case "FVLeiden"{if ($hap1 eq "NotDetected"){print OUT "FVLeiden looks good.\n";}else{print OUT "FVLeiden is INCORRECT.\n";}}
				case "MTHFR"{if (($hap1 eq "WT") && ($hap2 eq "WT")){print OUT "MTHFR looks good.\n";}else{print OUT "MTHFR is INCORRECT.\n";}}
				case "CYP3A4"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP3A4 looks good.\n";}else{print OUT "CYP3A4 is INCORRECT.\n";}}
				case "CYP3A5"{if (($hap1 eq "*3") && ($hap2 eq "*3")){print OUT "CYP3A5 looks good.\n";}else{print OUT "CYP3A5 is INCORRECT.\n";}}
				#case "Sample ID"{if ($hap1 eq "E3"){next;}}
				#else{print OUT "$pos_control_id is INCORRECT.\n";}
			}
		}
		case 'POS-NA17280'{
			switch($gene){
				case "APOE"{if (($hap1 eq "E3") && ($hap2 eq "E3")){print OUT "APOE looks good.\n";}else{print OUT "APOE is INCORRECT.\n";}}
				case "CYP2C19"{if (($hap1 eq "*1") && ($hap2 eq "*8")){print OUT "CYP2C19 looks good.\n";}else{print OUT "CYP2C19 is INCORRECT.\n";}}
				case "CYP2C9"{if (($hap1 eq "*1") && ($hap2 eq "*2")){print OUT "CYP2C9 looks good.\n";}else{print OUT "CYP2C9 is INCORRECT.\n";}}
				case "CYP2D6"{if (($hap1 eq "*2") && ($hap2 eq "*3")){print OUT "CYP2D6 looks good.\n";}else{print OUT "CYP2D6 is INCORRECT.\n";}}
				case "VKORC1"{if (($hap1 eq "*B") && ($hap2 eq "*B")){print OUT "VKORC1 looks good.\n";}else{print OUT "VKORC1 is INCORRECT.\n";}}
				case "FactorII"{if ($hap1 eq "NotDetected"){print OUT "FactorII looks good.\n";}else{print OUT "FactorII is INCORRECT.\n";}}
				case "FVLeiden"{if ($hap1 eq "NotDetected"){print OUT "FVLeiden looks good.\n";}else{print OUT "FVLeiden is INCORRECT.\n";}}
				case "MTHFR"{if (($hap1 eq "WT") && ($hap2 eq "Mut")){print OUT "MTHFR looks good.\n";}else{print OUT "MTHFR is INCORRECT.\n";}}
				case "CYP3A4"{if (($hap1 eq "*1") && ($hap2 eq "*1")){print OUT "CYP3A4 looks good.\n";}else{print OUT "CYP3A4 is INCORRECT.\n";}}
				case "CYP3A5"{if (($hap1 eq "*3") && ($hap2 eq "*3")){print OUT "CYP3A5 looks good.\n";}else{print OUT "CYP3A5 is INCORRECT.\n";}}
				#case "Sample ID"{if ($hap1 eq "E3"){next;}}
				#else{print OUT "$pos_control_id is INCORRECT.\n";}
			}
		}
		
	} # end of BIG switch

} # end of BIG for loop


close(OUT);			
