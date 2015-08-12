package SAMPLE_ID;
###############################################
# Company: American International Biotechnology 
# Author: Jeri Dilts 
# Date: March 31st, 2014
###############################################

#!/usr/bin/perl -w


use strict;
use CGI;
use Exporter;
our @ISA = qw(Exporter); 
our @EXPORT = qw(report_sampleID);


##### Prototypes ##############################
sub report_sampleID($$$$$$); 
sub evaluateID($$); 
sub evaluateID_RC($$); 
sub evaluateGENDER($$);



















###############################################
##### Subroutines #############################
sub report_sampleID($$$$$$){

	my($cov5500,$temp_sample,$reads,$snpcov,$plugin_filepath,$barcode) = @_;
	my %temp_sample = %{$temp_sample};
	
	#-----VARIABLES-----
	my ($ID_ATP13A4,$ID_PALLD,$ID_ADAMTS2,$ID_PTN,$ID_TRDMT1,$ID_LEPREL2,$ID_RAB31,$ID_SYN3,$GENDER);
	my $cov_ATP13A4 = 0;
	my $cov_PALLD = 0;
	my $cov_ADAMTS2 = 0;
	my $cov_PTN = 0;
	my $cov_TRDMT1 = 0;
	my $cov_LEPREL2 = 0;
	my $cov_RAB31 = 0;
	my $cov_SYN3 = 0;
	my $result;
	
	#-----VARIANT CHECK-----
	if(exists $temp_sample{"ATP13A4,500"}){($ID_ATP13A4,$cov_ATP13A4) = evaluateID($temp_sample{"ATP13A4,500"},$snpcov)};
	if(exists $temp_sample{"PALLD,500"}){($ID_PALLD,$cov_PALLD) = evaluateID($temp_sample{"PALLD,500"},$snpcov)};
	if(exists $temp_sample{"ADAMTS2,500"}){($ID_ADAMTS2,$cov_ADAMTS2) = evaluateID_RC($temp_sample{"ADAMTS2,500"},$snpcov)};
	if(exists $temp_sample{"PTN,500"}){($ID_PTN,$cov_PTN) = evaluateID($temp_sample{"PTN,500"},$snpcov)};
	if(exists $temp_sample{"TRDMT1,500"}){($ID_TRDMT1,$cov_TRDMT1) = evaluateID_RC($temp_sample{"TRDMT1,500"},$snpcov)};
	if(exists $temp_sample{"LEPREL2,500"}){($ID_LEPREL2,$cov_LEPREL2) = evaluateID($temp_sample{"LEPREL2,500"},$snpcov)};
	if(exists $temp_sample{"RAB31,500"}){($ID_RAB31,$cov_RAB31) = evaluateID($temp_sample{"RAB31,500"},$snpcov)};
	if(exists $temp_sample{"SYN3,500"}){($ID_SYN3,$cov_SYN3) = evaluateID($temp_sample{"SYN3,500"},$snpcov)};
	if(exists $temp_sample{"AmelY,487"}){($GENDER) = evaluateGENDER($temp_sample{"AmelY,487"},$snpcov);}else{$GENDER="F";}
	
	#print Dumper(\%temp_sample);
	
	open (FH, ">>$plugin_filepath/AIB_RESULTS/sampleID.txt") || die $!;
	chomp $reads;
	my @sampleID;
	if(($reads >= $cov5500)&&($cov_ATP13A4==1)&&($cov_PALLD==1)&&($cov_ADAMTS2==1)&&($cov_PTN==1)&&($cov_TRDMT1==1)&&($cov_LEPREL2==1)&&($cov_RAB31==1)&&($cov_SYN3==1)){
		my $gender_id = "$GENDER-".$ID_ATP13A4.$ID_PALLD.$ID_ADAMTS2.$ID_PTN.$ID_TRDMT1.$ID_LEPREL2.$ID_RAB31.$ID_SYN3;
		print FH "$barcode,$gender_id\n";
		my $id = "$ID_ATP13A4"."$ID_PALLD"."$ID_ADAMTS2"."$ID_PTN"."$ID_TRDMT1"."$ID_LEPREL2"."$ID_RAB31"."$ID_SYN3";
		$result = "Sample ID||$GENDER|$id|||||||\n";
		return ($result);
	}
	if($reads < $cov5500){
		my $id = "$ID_ATP13A4"."$ID_PALLD"."$ID_ADAMTS2"."$ID_PTN"."$ID_TRDMT1"."$ID_LEPREL2"."$ID_RAB31"."$ID_SYN3";
		my $error = "not enough sample coverage";
		print FH "$barcode,"."$error\n";
		$result = "Sample ID|$error|$GENDER|$id|||||||\n";
		return($result);
	}
	close(FH);
}

sub evaluateID($$){ 
my ($line,$snpcov)=@_; 
my ($ratioA,$ratioC,$ratioG,$ratioT,$ratioDEL,$reads,$id); 
my $VarCovFLAG = 0; #variant coverage Flag my @VAR = split(/\t/,$line);
my @VAR = @$line;

    #checks coverage
	if ($VAR[1] >= $snpcov){
		
		$VarCovFLAG = 1;
		$ratioA = int(($VAR[2]/$VAR[1])*100);
		$ratioC = int(($VAR[3]/$VAR[1])*100);
		$ratioG = int(($VAR[4]/$VAR[1])*100);
		$ratioT = int(($VAR[5]/$VAR[1])*100);
			
		if ($ratioA >= 90){$id='A';}
		if ($ratioC >= 90){$id='C';}
		if ($ratioG >= 90){$id='G';}
		if ($ratioT >= 90){$id='T';}
		
		if (($ratioA<90) && ($ratioA>10) && ($ratioC<90) && ($ratioC>10)){$id='M';}
		if (($ratioA<90) && ($ratioA>10) && ($ratioG<90) && ($ratioG>10)){$id='R';}
		if (($ratioA<90) && ($ratioA>10) && ($ratioT<90) && ($ratioT>10)){$id='W';}
		if (($ratioC<90) && ($ratioC>10) && ($ratioG<90) && ($ratioG>10)){$id='S';}
		if (($ratioC<90) && ($ratioC>10) && ($ratioT<90) && ($ratioT>10)){$id='Y';}
		if (($ratioG<90) && ($ratioG>10) && ($ratioT<90) && ($ratioT>10)){$id='K';}
	}
	return ($id,$VarCovFLAG);
}
sub evaluateID_RC($$){ 
my ($line,$snpcov)=@_;
my ($ratioA,$ratioC,$ratioG,$ratioT,$ratioDEL,$reads,$id);
my $VarCovFLAG = 0; #variant coverage Flag my @VAR = split(/\t/,$line);
my @VAR = @$line;

    #checks coverage
	if ($VAR[1] >= $snpcov){
		
		$VarCovFLAG = 1;
		$ratioA = int(($VAR[2]/$VAR[1])*100);
		$ratioC = int(($VAR[3]/$VAR[1])*100);
		$ratioG = int(($VAR[4]/$VAR[1])*100);
		$ratioT = int(($VAR[5]/$VAR[1])*100);
		
		if ($ratioA >= 90){$id='T';}
		if ($ratioC >= 90){$id='G';}
		if ($ratioG >= 90){$id='C';}
		if ($ratioT >= 90){$id='A';}
		
		if (($ratioA<90) && ($ratioA>10) && ($ratioC<90) && ($ratioC>10)){$id='K';}
		if (($ratioA<90) && ($ratioA>10) && ($ratioG<90) && ($ratioG>10)){$id='Y';}
		if (($ratioA<90) && ($ratioA>10) && ($ratioT<90) && ($ratioT>10)){$id='W';}
		if (($ratioC<90) && ($ratioC>10) && ($ratioG<90) && ($ratioG>10)){$id='S';}
		if (($ratioC<90) && ($ratioC>10) && ($ratioT<90) && ($ratioT>10)){$id='R';}
		if (($ratioG<90) && ($ratioG>10) && ($ratioT<90) && ($ratioT>10)){$id='M';}
		}
	return ($id,$VarCovFLAG);
}
	
sub evaluateGENDER($$){

my($line,$snpcov)=@_;
my $gender;
my @VAR = @$line;

        #checks coverage
        if ($VAR[1] >= $snpcov){$gender = "M";}
		else {$gender="F";} return($gender);
}
	
	
1;
