##############################################################
### printTorrent.pm                                        ###
### April 3,2014                                           ###
### AIBiotech LLC                                          ###
### William Budd; PhD                                      ###
### Version: 2.00                                          ###
### Improvements Needed:                                   ###
##############################################################
package printTorrent;
use strict; 
use warnings;
use Exporter;
use Time::Piece;
our @ISA = qw(Exporter);
our @EXPORT = qw(printHeader printCalls printAccession);

########## Sub-routine prototypes ##############################
sub printHeader($$$);
sub printCalls($$);
sub convertMonth($);
sub printAccession($$);

1;# return statement required for PERL modules!

########### MAIN METHOD ###########################################
sub printHeader($$$){
	my($machine, $version, $file) = @_;
	open(my $OUT,">$file") || die $!;
## Variable declaration 
	my($processDate );	
	$processDate = Time::Piece->new ->mdy;
	print $OUT "##### Ion Torrent| $machine\n";
	print $OUT "##### Assay Type| Ion Torrent Panel\n";
	print $OUT "##### AIBiotech Version| $version\n";
	print $OUT "##### Ion Torrent Protocol| Ampliseq\n";	
	print $OUT "##### Run Date| $processDate\n";
	print $OUT "##### Pipeline process date| $processDate\n";
	print $OUT "##### Input File|$file\n";
	print $OUT "##### Positive Control|\n";
	print $OUT "##### Negative Control|\n\n";
	close ($OUT);
}# end printHeader

sub printAccession($$){
my ($accession, $OUT) = @_;
	print $OUT "$accession\n";
}

sub printCalls($$){
my($callRef, $OUT) = @_;
##variable declarations
my(@calls);
	@calls = @$callRef;
	foreach my $call (@calls){
		print $OUT "$call";
	}
	print "\n";
}# end sub printCalls


sub convertMonth($){
my($monthIn) = @_;
##variables
my($outMon);
	my %month = ( 	"Jan" => "01",
					"Feb" => "02",
					"Mar" => "03",
					"Apr" => "04",
					"May" => "05",
					"Jun" => "06",
					"Jul" => "07",
					"Aug" => "09",
					"Sep" => "10",
					"Oct" => "11",
					"Dec" => "12",
				);
	$outMon= $month{$monthIn};
return $outMon;	
}# end convertMonth sub



