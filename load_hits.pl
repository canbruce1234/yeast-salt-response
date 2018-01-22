#!/usr/bin/perl -w
##########################################################################################
#
#	input: experiment name, file name of hitlist from a chipchip expt.
#	output: insert into "allhits" table expt name, and other data from hitlist
#
##########################################################################################
use DBI;

$expt   = $ARGV[0];    # expt name
$infile = $ARGV[1];    # hitlist

open( IN, $infile ) || die "can't open $infile: $!";

$dbh = DBI->connect( "dbi:mysql:chipchip", "cbruce", "yalejohn" )
  or die "Can't connect to mysql database: $DBI::errstr\n";

while ( defined( $line = <IN> ) ) {
	if ( $line =~ /^#/ ) { next }
	if ( $line =~ /^ / ) { next }
	@fields = split /\t/, $line;
	$chr    = $fields[0];
	$chr    = substr( $chr, 3 );
	$shoulder = $fields[1];
	$start  = $fields[2];
	$peak	= $fields[3];
	$end    = $fields[4];
	$signal = $fields[7];
	$geneW1  = $fields[9];
	$akaW1   = $fields[10];
	$statusW1= $fields[11];
	$distW1 = $fields[12];
	$geneW2  = $fields[13];
	$akaW2   = $fields[14];
	$statusW2= $fields[15];
	$distW2 = $fields[16];
	$geneC1  = $fields[17];
	$akaC1   = $fields[18];
	$statusC1= $fields[19];
	$distC1 = $fields[20];
	$geneC2  = $fields[21];
	$akaC2   = $fields[22];
	$statusC2= $fields[23];
	$distC2 = $fields[24];
        $orf    = $fields[25];

	if ( !defined $geneC1 ) { $geneC1 = "" }
	if ( !defined $geneW1 ) { $geneW1 = "" }
	if ( !defined $akaC1 )  { $akaC1  = "" }
	if ( !defined $akaW1 )  { $akaW1  = "" }
	if ( !defined $geneC2 ) { $geneC2 = "" }
	if ( !defined $geneW2 ) { $geneW2 = "" }
	if ( !defined $akaC2 )  { $akaC2  = "" }
	if ( !defined $akaW2 )  { $akaW2  = "" }
	
	# to handle special case where aka name has a single quote...
	if ($akaC1=~ /'/){$akaC1=~ s/'/\'/g} 
	if ($akaC2=~ /'/){$akaC2=~ s/'/\'/g}
	

#second field is for FDR; will be filled in afterwards
	$stmt1 = $dbh->prepare(
		qq{
    	insert into allhits values (
		'$expt', '', $chr, '$shoulder', $start, $peak, $end, $signal, 
		'$geneW1', "$akaW1", '$statusW1','$distW1', 
		'$geneW2', "$akaW2", '$statusW2','$distW2',
		'$geneC1', "$akaC1", '$statusC1','$distC1', 
		'$geneC2', "$akaC2", '$statusC2','$distC2',$orf )
        }
	  )
	  or die "can't prepare stmt: $DBI::errstr";

	$stmt1->execute() or die "Can't execute SQL statement: $DBI::errstr\n";

}
close IN;

# Now, disconnect from the database
$dbh->disconnect or warn "Disconnection failed: $DBI::errstr\n";

