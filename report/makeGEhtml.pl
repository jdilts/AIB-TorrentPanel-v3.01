#////////////////////////////////////////////////////////////////////////////////////////////
# Company:	American International Biotechnology, LLC
# Author:	Jeri Dilts
# Date:		December 2013
#////////////////////////////////////////////////////////////////////////////////////////////

#!/usr/bin/perl

use strict;
use warnings;
use CGI;
use File::Slurp;

my ($plugin_filepath,$genotypeTXT) = @ARGV;

#-----PROTOTYPES--------------------------------------
sub makelink($);
sub genotypeHTMLcreator($$);
#-----------------------------------------------------

#-----EXECUTE-----------------------------------------
makelink($plugin_filepath);
genotypeHTMLcreator($plugin_filepath,$genotypeTXT);
#-----------------------------------------------------





#////////////////////////////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////////////////////////////
#-----SUBROUTINES----------------------------------------------------------------------------
#////////////////////////////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////////////////////////////

sub makelink($){

	my $plugin_filepath = shift;
	
	system("ln -s /usr/local/sbin/calls_style.css $plugin_filepath/calls_style.css");
}


sub genotypeHTMLcreator($$){

	my ($plugin_filepath,$genotypeTXT) = @_;
	
	my $cgi = CGI->new(); #CGI object

	open (HTML,">","$plugin_filepath/genotypes.html") || die $!;

	#-----Building HTML <header>---------------------------
	print HTML $cgi->start_html(
	
		-title=>'Calls.html',
		-author=>'Jeri Dilts @ AIBIOTECH',
		-style=>[{-src=>'calls_style.css'}],
		
	);
	#------------------------------------------------------
	
	my @GENOTYPES = read_file($genotypeTXT);
	
	foreach(@GENOTYPES){

		
		unless($_=~/^$/){
			print HTML $cgi->p($_);
		} else{
			print HTML $cgi->br;
			print HTML $cgi->p($_);
		}	
	}
	
	print HTML $cgi->end_html;

}