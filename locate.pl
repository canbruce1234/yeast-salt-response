# This subroutine takes chromosomal coordinates of a ChIP chip "hit"
# and returns gene names, common names(if available), and description
# for both strands where the hit is within 1500 nucleotides upstream
# of the ORF start. It also returns a 1 or 0, depending on whether
# the input coordinates are completely within an ORF (on either strand),
# or not, respectively.
#
# The first subroutine, "locate_run_yeast" does this for yeast.
# A similar routine does it for the encode genes.
#
# to use this subroutine, type the following near the top of your
# calling script:
#	do '/home/cbruce/eclipse/workspace/ChipChip/locate.pl';
use DBI;

sub locate_run_yeast {
	my (@run) = @_;

	$chr2         = $run[0];
	$start_of_run = $run[1];
	$end_of_run   = $run[2];
	$chr          = substr( $chr2, 3 );
	if ( $chr eq "M" ) { $chr = 17 }
	$middle_of_run = int( ( $run[1] + $run[2] ) / 2 );
	@result        = ();

	# print "$chr/$start_of_run/$end_of_run\n";

	my $dbh = DBI->connect( "dbi:mysql:chipchip", "cbruce", $password)
	  or die "Can't connect to mysql database: $DBI::errstr\n";

	# print "$chr2/$chr/$middle_of_run \n";
	my $stmt1 = $dbh->prepare(
		qq{
                        SELECT name, aka, $middle_of_run - start, descr
                          FROM yeast g
	 					 WHERE upstr_end <= $middle_of_run  
           				   AND end > $middle_of_run
                           AND g.strand = 'W'
                           AND g.vvv= '$chr'
                           AND g.type = 'ORF'
                           LIMIT 1
                        }
	  )
	  or die "can't prepare stmt: $DBI::errstr";

	my $stmt2 = $dbh->prepare(
		qq{
                        SELECT name, aka, start - $middle_of_run, descr
                          FROM yeast g
                		 WHERE end <= $middle_of_run  
 						   AND upstr_end > $middle_of_run 
                           AND g.strand = 'C'
                           AND g.vvv= '$chr'
                           AND g.type = 'ORF'
                         LIMIT 1
                        }
	  )
	  or die "can't prepare stmt: $DBI::errstr";

	my $stmt3 = $dbh->prepare(
		qq{
                        SELECT count(*)
                          FROM yeast g
                         WHERE start <  $start_of_run
                           AND end > $end_of_run
                           AND g.strand = 'W'
                           AND g.vvv= '$chr'
                           AND g.type = 'ORF'
                        }
	  )
	  or die "can't prepare stmt: $DBI::errstr";

	my $stmt4 = $dbh->prepare(
		qq{
                        SELECT count(*)
                          FROM yeast g
                         WHERE  start > $end_of_run
                           AND end < $start_of_run
                           AND  g.strand = 'C'
                           AND g.vvv= '$chr'
                           AND g.type = 'ORF'
                        }
	  )
	  or die "can't prepare stmt: $DBI::errstr";

	my $stmt5 = $dbh->prepare(
		qq{
                        SELECT name
                          FROM yeast g
                         WHERE (
                            (upstr_end <= $start_of_run AND end > $end_of_run)
                            OR 
                            (upstr_end > $start_of_run AND end <= $end_of_run)
                            OR
                            (upstr_end between $start_of_run AND $end_of_run)
                            OR
                            (end between $start_of_run AND $end_of_run)
                            )
                           AND g.strand = 'W'
                           AND g.vvv= '$chr'
                           AND g.type = 'ORF'
                        }
	  )
	  or die "can't prepare stmt: $DBI::errstr";

	my $stmt6 = $dbh->prepare(
		qq{
                        SELECT name
                          FROM yeast g
                         WHERE (
                            (upstr_end <= $end_of_run AND end > $start_of_run)
                            OR 
                            (upstr_end > $end_of_run AND end <= $start_of_run)
                            OR
                            (upstr_end between $end_of_run AND $start_of_run)
                            OR
                            (end between $end_of_run AND $start_of_run)
                            )
                           AND g.strand = 'C'
                           AND g.vvv= '$chr'
                           AND g.type = 'ORF'
                        }
	  )
	  or die "can't prepare stmt: $DBI::errstr";

	$stmt1->execute() or die "Can't execute SQL statement: $DBI::errstr\n";
	$stmt2->execute() or die "Can't execute SQL statement: $DBI::errstr\n";
	$stmt3->execute() or die "Can't execute SQL statement: $DBI::errstr\n";
	$stmt4->execute() or die "Can't execute SQL statement: $DBI::errstr\n";
	$stmt5->execute() or die "Can't execute SQL statement: $DBI::errstr\n";
	$stmt6->execute() or die "Can't execute SQL statement: $DBI::errstr\n";

	my @row1 = $stmt1->fetchrow_array();
	warn "Data fetching terminated early by error: $DBI::errstr\n"
	  if $DBI::err;
	if ( scalar(@row1) == 0 ) { @row1 = ( "", "", "", "" ) }

	my @row2 = $stmt2->fetchrow_array();
	warn "Data fetching terminated early by error: $DBI::errstr\n"
	  if $DBI::err;
	if ( scalar(@row2) == 0 ) { @row2 = ( "", "", "", "" ) }

	my @row3 = $stmt3->fetchrow_array();
	warn "Data fetching terminated early by error: $DBI::errstr\n"
	  if $DBI::err;

	my @row4 = $stmt4->fetchrow_array();
	warn "Data fetching terminated early by error: $DBI::errstr\n"
	  if $DBI::err;

	my @row5 = $stmt5->fetchrow_array();
	warn "Data fetching terminated early by error: $DBI::errstr\n"
	  if $DBI::err;

	my @row6 = $stmt6->fetchrow_array();
	warn "Data fetching terminated early by error: $DBI::errstr\n"
	  if $DBI::err;

	if ( $row3[0] != 0 or $row4[0] != 0 ) { $orf = 1; }
	else { $orf = 0 }

	$genes1 = join( ",", @row5 );
	$genes2 = join( ",", @row6 );
	if ( scalar(@row5) > 0 and scalar(@row6) > 0 ) { $sep = "," }
	else { $sep = "" }
	$genes = $genes1 . $sep . $genes2;

	push @result, @row1;
	push @result, @row2;
	push @result, ( $genes, $orf );

	$stmt1->finish;
	$stmt2->finish;
	$stmt3->finish;
	$stmt4->finish;
	$stmt5->finish;
	$stmt6->finish;

	# Now, disconnect from the database
	$dbh->disconnect or warn "Disconnection failed: $DBI::errstr\n";
	return @result;
}

sub locate_run_yeast2 {
	#this subroutine returns an array or arrays
	my (@runs) = @_;

	my $dbh = DBI->connect( "dbi:mysql:chipchip", "cbruce", "yalejohn" )
	  or die "Can't connect to mysql database: $DBI::errstr\n";

	for my $j ( 0 .. $#runs ) {

		$chr2         = $runs[$j][0];
		$start_of_run = $runs[$j][1];
		$end_of_run   = $runs[$j][2];
		$chr          = substr( $chr2, 3 );
		if ( $chr eq "M" ) { $chr = 17 }
		$middle_of_run = int( ( $runs[$j][1] + $runs[$j][2] ) / 2 );
		@result        = ();

		# print "$runs[$j][0]/$runs[$j][1]/$runs[$j][2]\n";
		# print "$chr2/$chr/$middle_of_run \n";
		my $stmt1 = $dbh->prepare(
			qq{
    			SELECT name, aka, $middle_of_run - start
          		  FROM yeast g
	 		 WHERE upstr_end <= $middle_of_run  
           		   AND end > $middle_of_run
			   AND g.strand = 'W' 
		           AND g.vvv= '$chr'
           		   AND g.type = 'ORF'
                           LIMIT 1
                        }
		  )
		  or die "can't prepare stmt: $DBI::errstr";

		my $stmt2 = $dbh->prepare(
			qq{
               		SELECT name, aka, start - $middle_of_run 
                 	  FROM yeast g
                	 WHERE end <= $middle_of_run  
 			   			AND upstr_end > $middle_of_run 
			   			AND g.strand = 'C'  
                  	   AND g.vvv= '$chr'
                  	   AND g.type = 'ORF'
                         LIMIT 1
                        }
		  )
		  or die "can't prepare stmt: $DBI::errstr";

		my $stmt3 = $dbh->prepare(
			qq{
                        SELECT count(*)
                          FROM yeast g
                         WHERE start <  $start_of_run
                           AND end > $end_of_run
                           AND g.strand = 'W'
                           AND g.vvv= '$chr'
                           AND g.type = 'ORF'
                        }
		  )
		  or die "can't prepare stmt: $DBI::errstr";

		my $stmt4 = $dbh->prepare(
			qq{
                        SELECT count(*)
                          FROM yeast g
                         WHERE  start > $end_of_run
                           AND end < $start_of_run
                           AND  g.strand = 'C'
                           AND g.vvv= '$chr'
                           AND g.type = 'ORF'
                        }
		  )
		  or die "can't prepare stmt: $DBI::errstr";

		my $stmt5 = $dbh->prepare(
			qq{
                        SELECT name
                          FROM yeast g
                         WHERE (
                            (upstr_end <= $start_of_run AND end > $end_of_run)
                            OR 
                            (upstr_end > $start_of_run AND end <= $end_of_run)
                            OR
                            (upstr_end between $start_of_run AND $end_of_run)
                            OR
                            (end between $start_of_run AND $end_of_run)
                            )
                           AND g.strand = 'W'
                           AND g.vvv= '$chr'
                           AND g.type = 'ORF'
                        }
		  )
		  or die "can't prepare stmt: $DBI::errstr";

		my $stmt6 = $dbh->prepare(
			qq{
                        SELECT name
                          FROM yeast g
                         WHERE (
                            (upstr_end <= $end_of_run AND end > $start_of_run)
                            OR 
                            (upstr_end > $end_of_run AND end <= $start_of_run)
                            OR
                            (upstr_end between $end_of_run AND $start_of_run)
                            OR
                            (end between $end_of_run AND $start_of_run)
                            )
                           AND g.strand = 'C'
                           AND g.vvv= '$chr'
                           AND g.type = 'ORF'
                        }
		  )
		  or die "can't prepare stmt: $DBI::errstr";

		$stmt1->execute() or die "Can't execute SQL statement: $DBI::errstr\n";
		$stmt2->execute() or die "Can't execute SQL statement: $DBI::errstr\n";
		$stmt3->execute() or die "Can't execute SQL statement: $DBI::errstr\n";
		$stmt4->execute() or die "Can't execute SQL statement: $DBI::errstr\n";
		$stmt5->execute() or die "Can't execute SQL statement: $DBI::errstr\n";
		$stmt6->execute() or die "Can't execute SQL statement: $DBI::errstr\n";

		my @row1 = $stmt1->fetchrow_array();
		warn "Data fetching terminated early by error: $DBI::errstr\n"
		  if $DBI::err;
		if ( scalar(@row1) == 0 ) { @row1 = ( " ", " ", "" ) }

		my @row2 = $stmt2->fetchrow_array();
		warn "Data fetching terminated early by error: $DBI::errstr\n"
		  if $DBI::err;
		if ( scalar(@row2) == 0 ) { @row2 = ( "", "", "" ) }

		my @row3 = $stmt3->fetchrow_array();
		warn "Data fetching terminated early by error: $DBI::errstr\n"
		  if $DBI::err;

		my @row4 = $stmt4->fetchrow_array();
		warn "Data fetching terminated early by error: $DBI::errstr\n"
		  if $DBI::err;

		my @genes=();
		while ( $row5 = $stmt5->fetchrow_arrayref ) {
		   if ( defined $row5 ) {push @genes, @$row5[0]} ;
		}

		warn "Data fetching terminated early by error: $DBI::errstr\n"
		  if $DBI::err;
		  
		#@genes6=();
		while ( $row6 = $stmt6->fetchrow_arrayref ) {
		   if ( defined $row6 ) {push @genes,  @$row6[0]}
		}
		
		warn "Data fetching terminated early by error: $DBI::errstr\n"
		  if $DBI::err;

		if ( $row3[0] != 0 or $row4[0] != 0 ) { $orf = 1; }
		else { $orf = 0 }

		$genes = join ",", @genes;
		# print "genes $genes \n";

		push @result, ( $chr, $start_of_run, $end_of_run );
		push @result, @row1;
		push @result, @row2;
		push @result, ( $genes, $orf );


		$stmt1->finish;
		$stmt2->finish;
		$stmt3->finish;
		$stmt4->finish;
		$stmt5->finish;
		$stmt6->finish;

		push @{ $results[$j] }, @result;

	}

	# Now, disconnect from the database
	$dbh->disconnect or warn "Disconnection failed: $DBI::errstr\n";
	return @results;
}

