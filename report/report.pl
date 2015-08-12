#!/usr/bin/perl -w

## Author: Jeri Dilts
## Date: April 17th, 2014
## Company: American International Biotechnology
## Department: Bioinformatics
## This script checks for duplicate sampleIDs and reports them via email weekly


use strict;
use File::Slurp;
use File::Find;
use CGI;

## prototypes
sub makelinks($$);
sub print_positive_control($$);
sub print_dropouts($$);
sub print_sample_cov($$);
sub link2_filtered_bams($$$$$);
sub print_genotypes($$);
sub print_nt_freq($$$);
sub ParameterExceptionHandling($$);

## variables
my ($plugin_filepath,$ip_address,$run,$plugin_name)= @ARGV;

## main
my $cgi = CGI->new(); #CGI object

open (my $HTML,">","$plugin_filepath/report.html") || die $!;

#-----Building HTML <header>---------------------------
print $HTML $cgi->start_html(

	-title=>'report.html',
	-author=>'Jeri Dilts @ AIBIOTECH',
	-style=>[
			{-src=>'http://code.jquery.com/ui/1.10.3/themes/smoothness/jquery-ui.css'},
			{-src=>'nucfreq_style.css'}],
	-script=>[
			{-type=>'JAVASCRIPT',-src=>'http://code.jquery.com/jquery-1.9.1.js'},
			{-type=>'JAVASCRIPT',-src=>'jquery-ui.js'}
			]
	);
#------------------------------------------------------

#AIBonlinelogo.jpg
print $HTML "<img src=\"logo.jpg\" alt=\"American International Biotechnology\" width=\"1111\" height=\"162\">";

print $HTML "<div id=\"accordion\">";

## execute
makelinks($plugin_filepath,$plugin_name);
print_positive_control($plugin_filepath,$HTML);
print_dropouts($plugin_filepath,$HTML);
my($samples)=print_sample_cov($plugin_filepath,$HTML);

print $HTML $cgi->p("Filtered Bams");
print $HTML "<div>";

my @samples = @$samples;
foreach(@samples){
	my $sample = $_;
	chomp $sample;
	link2_filtered_bams($ip_address,$run,$plugin_name,$HTML,$sample);
}
print $HTML "</div>";

print_genotypes($plugin_filepath,$HTML);
print_nt_freq($plugin_filepath,$HTML,$cgi);


#---------closes <div accordian>-----------------------
print $HTML "</div>";	

#---------javascript-----------------------------------
print $HTML "
<script>
\$( \"#accordion\" ).accordion({heightStyle: 'content'});
</script>
";

#---------footer---------------------------------------
print $HTML "
<footer>
<p>American International Biotechnology</p>
<p>Contact info: <a href=\"mailto:jdilts\@aibiotech.com\">Jeri Dilts</a> and <a href=\"mailto:wbudd\@aibiotech.com\">William Budd</a></p>
<p>Internal #: 6037 or 6044</p>
</footer>
";
	
print $HTML $cgi->end_html;
close($HTML);








## subroutines
sub makelinks($$){

	my ($plugin_filepath,$plugin_name) = @_;
	
	unless (-e "$plugin_filepath/jquery-ui.js"){
	system("ln -s /results/plugins/$plugin_name/style/jquery-ui.js $plugin_filepath/jquery-ui.js");}
	
	unless (-e "$plugin_filepath/nucfreq_style.css"){
	system("ln -s /results/plugins/$plugin_name/style/nucfreq_style.css $plugin_filepath/nucfreq_style.css");}
	
	unless (-e "$plugin_filepath/logo.jpg"){
	system("ln -s /results/plugins/$plugin_name/style/logo.jpg $plugin_filepath/logo.jpg");}
}

sub print_positive_control($$){

	my ($plugin_filepath,$HTML) = @_;
	my $pc_file = $plugin_filepath."/AIB_RESULTS/positive_control_check.txt";
	my @f = read_file($pc_file);

	print $HTML $cgi->p("Positive Control Check");
	print $HTML "<div>"; 
	print $HTML "<p>";
	if (@f) {
		foreach(@f){
			unless(($_=~/^#/) || ($_=~/^$/)){
				print $HTML "$_<br>";
			}
		}
	}
	print $HTML "</p>";
	print $HTML "</div>";
}

sub print_dropouts($$){

	my ($plugin_filepath,$HTML) = @_;
	my $diff = $plugin_filepath."/AIB_RESULTS/dropout_barcodes.txt";
	my @dropouts = read_file($diff);

	print $HTML $cgi->p("Dropouts");
	print $HTML "<div>"; 
	print $HTML "<p>";
	if (@dropouts) {
		foreach(@dropouts){
			unless(($_=~/^#/) || ($_=~/^$/)){
				print $HTML "Barcode $_ dropped out.<br>"
			}
		}
	} else {
		print $HTML "No dropouts.<br>";
	}
	print $HTML "<p>";
	print $HTML "</div>";

}

sub print_sample_cov($$){
	
	my ($plugin_filepath,$HTML) = @_;
	my @sample_cov = read_file("$plugin_filepath"."/AIB_RESULTS/sample_coverage.xls");
	my @samples;
	
	print $HTML $cgi->p("Sample Coverages");
	print $HTML "<div>";
	print $HTML "<p>";
	foreach(@sample_cov){
		unless(($_=~/^#/) || ($_=~/^$/)){
			print $HTML "$_<br>";
			my @line = split(/,/,$_);
			push(@samples, $line[0]);
		}
	}
	print $HTML "</p>";
	print $HTML "</div>";
	return(\@samples);
}

sub link2_filtered_bams($$$$$){

	my ($ip_address,$run,$plugin_name,$HTML,$sample) = @_;
	
	my $bam_link = "http://$ip_address/output/Home/$run/plugin_out/$plugin_name"."_out"."/$sample/$sample"."_rawlib.bam";
	my $bai_link = $bam_link.".bai";
		
	print $HTML "<a href=\"$bam_link\">$bam_link</a><br>";
	print $HTML "<a href=\"$bai_link\">$bai_link</a><br>";
}

sub print_genotypes($$){

	my ($plugin_filepath,$HTML) = @_;
	my @f = read_file("$plugin_filepath"."/AIB_RESULTS/genotype.txt");

	print $HTML $cgi->p("Genotypes");
	print $HTML "<div>"; 
	
	if (@f) {
		print $HTML "<p>";
		foreach(@f){
			unless(($_=~/^#/) || ($_=~/^$/)){
			
				if ($_=~/^Sample ID/){ print $HTML "$_<br><br>";}
				else{print $HTML "$_<br>";}
			}
		}
		print $HTML "</p>";
	}
	print $HTML "</div>";
}

sub print_nt_freq($$$) {

	my ($plugin_filepath,$HTML,$cgi) = @_;
	my @nt_freqs = read_file("$plugin_filepath"."/AIB_RESULTS/nt_freq.xls");
	
	for(my $i = 1; $i<=10;$i++){shift @nt_freqs;} #gets rid of header
	
	my $flag = 1;
	foreach (@nt_freqs){
	
			if (($flag == 1) || ($_ =~/^Chrom/) || ($_=~/^$/)){
			
				if ($flag == 1){
				
					print $HTML $cgi->p($_); #prints barcode/sample patient identifier
					print $HTML "<div>";
					$flag = 0;
				}
				
				if ($_ =~/^Chrom/){		
				
					print $HTML $cgi->start_table({-border=>undef}); #starts table for sample 
					print $HTML $cgi->TR($cgi->th(['Amplicon','Position','Reference','Coverage','A','A%','C','C%','G','G%','T','T%','Deletions','DEL%','Parameters'])); #prints header for sample
					$flag = 0;
				}
				
				if ($_=~/^$/) {	
				
					print $HTML $cgi->end_table();
					print $HTML "</div>"; #end accordian div at the end of the table
					$flag = 1;
				}
			}else{
				my @line = split (/\t/,$_);

				my $amplicon = $line[0];
				my $position = $line[1];		
				my $reference = $line[2];
				my $coverage = $line[3];
				my $Areads = $line[4];
				my $Creads = $line[5];
				my $Greads = $line[6];
				my $Treads = $line[7];
				my $Deletions = $line[8];
					
				if ($coverage != 0){
				
					my $Apercent = sprintf("%.2f",($Areads/$coverage) * 100);
					my $Cpercent = sprintf("%.2f",($Creads/$coverage) * 100);
					my $Gpercent = sprintf("%.2f",($Greads/$coverage) * 100);
					my $Tpercent = sprintf("%.2f",($Treads/$coverage) * 100);
					my $DELpercent = sprintf("%.2f",($Deletions/$coverage) * 100);
						
					my ($parameter) = ParameterExceptionHandling($amplicon,$position);
						
					print $HTML $cgi->TR($cgi->td([$amplicon,$position,$reference,$coverage,$Areads,$Apercent,$Creads,$Cpercent,$Greads,$Gpercent,$Treads,$Tpercent,$Deletions,$DELpercent,$parameter]));
				}
			}
		
	}
}

sub ParameterExceptionHandling($$){

	my ($amplicon,$position) = @_;
	my $flag = 0;
	my $parameter = "80/20";

	if(($amplicon eq "CYP2C19s4") && ($position == 76)){$flag = 1;}
	if($amplicon eq "VKORC1"){$flag =2;}
	if ($flag==1){$parameter = "75/25";}
	if ($flag==2){$parameter = "70/30";}
	
	return $parameter;
}
