#!/usr/bin/perl -w
##########################################################################################
#
#	input: experiment name in hits table
#	output: insert into YAP4times table chr#, start pos, end pos, normalized
#
#
##########################################################################################
use DBI;

$expt = $ARGV[0];    # expt name
$col  = $ARGV[1];    # col name in table

$dbh = DBI->connect( "dbi:mysql:chipchip", "cbruce", $password )
  or die "Can't connect to mysql database: $DBI::errstr\n";

$stmt2 = $dbh->prepare(
	qq{
    	select chr, start, end, signal
    	from allhits where expt = '$expt'
        }
  )
  or die "can't prepare stmt: $DBI::errstr";

$stmt2->execute()
  or die "Can't execute SQL statement: $DBI::errstr\n";

while ( my $row2 = $stmt2->fetchrow_arrayref ) {
	( $chr, $start, $end, $signal) 
    	= @$row2;

	$stmt5 = $dbh->prepare(
		qq{
    	select chr, start, end from yap4chipexpr
    	where chr=$chr and
    	(start between $start-350 and $end+350
    	 or end between $start-350 and $end+350)
        }
	  )    
	  or die "can't prepare stmt: $DBI::errstr";

	$stmt5->execute()
	  or die "Can't execute SQL statement: $DBI::errstr\n";
	@row5 = $stmt5->fetchrow_array();
	warn "Data fetching terminated early by error: $DBI::errstr\n"
	  if $DBI::err;
	  
	if ( scalar(@row5) > 0 ) {
		$this_start = $row5[1];
		$this_end = $row5[2];
		
		$stmt4 = $dbh->prepare(
			qq{
    			update yap4chipexpr set $col=$signal
    			where chr=$chr and start=$this_start and end=$this_end
        		}
		  )  or die "can't prepare stmt: $DBI::errstr";
		$stmt4->execute()
		  or die "Can't execute SQL statement: $DBI::errstr\n";
	}

}

