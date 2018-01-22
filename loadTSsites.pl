#!/usr/bin/perl -w
##########################################################################################
#
#	input: experiment name in hits30min table
#	output: insert into TFsites table chr#, start pos, end pos, normalized
#              geneW1, geneW2, geneC1, geneC2  
#
##########################################################################################
use DBI;

$expt = $ARGV[0];    # expt name
$col  = $ARGV[1];    #col name in table
#$baseline = $ARGV[2];    # baseline shift so histogram mode = 0

$dbh = DBI->connect( "dbi:mysql:chipchip", "cbruce", $password )
  or die "Can't connect to mysql database: $DBI::errstr\n";

$stmt2 = $dbh->prepare(
	qq{
    	select chr, start, peak, end, signal,
    	geneW1, akaW1, statusW1, distW1, 
    	geneW2, akaW2, statusW2, distW2, 
    	geneC1, akaC1, statusC1, distC1, 
    	geneC2, akaC2, statusC2, distC2
    	from hits30min where expt = '$expt'
        }
  )
  or die "can't prepare stmt: $DBI::errstr";

$stmt2->execute()
  or die "Can't execute SQL statement: $DBI::errstr\n";

while ( my $row2 = $stmt2->fetchrow_arrayref ) {
	( $chr, $start, $peak, $end, $signal,
	  $geneW1, $akaW1, $statusW1, $distW1, 
   	  $geneW2, $akaW2, $statusW2, $distW2, 
      $geneC1, $akaC1, $statusC1, $distC1, 
      $geneC2, $akaC2, $statusC2, $distC2, ) 
    	= @$row2;

	$stmt5 = $dbh->prepare(
		qq{
    	select chr, start, end from TFsites
    	where chr=$chr and
    	(start between $start-100 and $end+100
    	 or end between $start-100 and $end+100)
        }
	  )    
	  or die "can't prepare stmt: $DBI::errstr";

	$stmt5->execute()
	  or die "Can't execute SQL statement: $DBI::errstr\n";
	@row5 = $stmt5->fetchrow_array();
	warn "Data fetching terminated early by error: $DBI::errstr\n"
	  if $DBI::err;
	  
	# to handle special case where aka name has a single quote...
	if ($akaC1=~ /'/){$akaC1=~ s/'/''/g; } 
	if ($akaC2=~ /'/){$akaC2=~ s/'/''/g; }
	
	if ( scalar(@row5) == 0 ) {

	
		$stmt3 = $dbh->prepare(
			qq{
    	insert into TFsites 
    		( chr, start, peak, end, $col,
    	      geneW1, akaW1, statusW1, distW1, 
    		  geneW2, akaW2, statusW2, distW2, 
    		  geneC1, akaC1, statusC1, distC1, 
    		  geneC2, akaC2, statusC2, distC2)
    	values ( $chr, $start, $peak, $end, $signal,
    			'$geneW1', '$akaW1', '$statusW1', '$distW1', 
   	  			'$geneW2', '$akaW2', '$statusW2', '$distW2', 
      			'$geneC1', '$akaC1', '$statusC1', '$distC1', 
      			'$geneC2', '$akaC2', '$statusC2', '$distC2' ) 
        }
		  )
		  or die "can't prepare stmt: $DBI::errstr";
		$stmt3->execute()
		  or die "Can't execute SQL statement: $DBI::errstr\n";
	} else {

		$this_start = $row5[1];
		$this_end = $row5[2];
		
		$stmt4 = $dbh->prepare(
			qq{
    	update TFsites set $col=$signal
    	where chr=$chr and start=$this_start and end=$this_end
        }
		  )
		  or die "can't prepare stmt: $DBI::errstr";
		$stmt4->execute()
		  or die "Can't execute SQL statement: $DBI::errstr\n";
	}

}

