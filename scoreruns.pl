#!/usr/bin/perl -w
#################################################################
# Reads four column output from R script that analyzes chipchip data;
# Identifies probes having log2Signal and -log10(p) value above thresholds.
# defines "runs" that are within maxgap1, then eliminates those runs
# that are less than minrun, them joins runs that are within maxgap2
# of each other.
#########################################################
use DBI;
$| = 1;

$ratio_threshold = $ARGV[0];   # threshold
$infile1         = $ARGV[1];   # chipchip expt output file
$outfile         = $ARGV[2];   # hitlist
$outfile2        = $ARGV[3];   # list of genes within hits or genes that are hit

#open( IN1,  $infile1 )     || die "can't open $infile1: $!";
open( OUT,  ">$outfile" )  || die "can't open $outfile: $!";
open( OUT2, ">$outfile2" ) || die "can't open $outfile2: $!";

# set program parameters:

#$ratio_threshold = 0.291;
$probe_length   = 50;          # size of a probe in nucleotides
$prob_threshold = 4
  ; # Minimum p-value for a probe to be part of a run, regardless of signal value.
$maxgap1 =
  120;   # combine runs/probes into runs if separated by less than this distance
$minrun   = 120;    # shortest length for a run
#$maxgap2  = 375;    # combine runs if separated by less than this distance
$maxgap2  = 750;    # combine runs if separated by less than this distance
$peak_sep = 800;    # minimum distance between two peaks
$peak_width = 400;	# width of peak to be reported
$window   =  1000;   			# window size for calculating average slope, 
  					# to determine peak top.

#Read the data file
( $chrRef, $posRef, $probRef, $ratioRef ) = &readfile($infile1);
our @chr   = @$chrRef;
our @pos   = @$posRef;
our @prob  = @$probRef;
our @ratio = @$ratioRef;

@results = &maxgap1( $maxgap1, $prob_threshold, $ratio_threshold );

@results = &minrun( $minrun, @results );

@results = &maxgap2( $maxgap2, @results );

@results = &calc_run(@results);


# get gene info from the chromosomal coordinates of the hits
@results = &locate_run(@results);

# get list of genes that are within hits or genes that are hit

&print_genes(@results);

#get some summary statistics
$avg_ratio = 0;
$avg_prob  = 0;
$sum_count = 0;
$sum_ratio = 0;
for $row (@results) {
	$sum_ratio = $sum_ratio + @{$row}[8] * @{$row}[9];
	$sum_prob  = $sum_prob + @{$row}[8] * @{$row}[10];
	$sum_count = $sum_count + @{$row}[8];
}

if ( $sum_count != 0 ) {
	$avg_ratio = sprintf( "%.2f", $sum_ratio / $sum_count );
	$avg_prob  = sprintf( "%.2f", $sum_prob / $sum_count );
	$avg_run   = sprintf( "%.0f", $sum_count / scalar(@results) );
} else {
	$avg_ratio = "''";
	$avg_prob  = "''";
	$avg_run   = "''";
}
$hit_cnt = scalar(@results);

#sort list by signal strength:
@results = sort { @$b[9] <=> @$a[9] } @results;

# print to stdout:
print
"********************************************************************************\n";
print "$infile1         threshold=$ratio_threshold\n";
print "Tsig\tTprob\tminrun\tmaxgap1\tmaxgap2\tLength\t",
  "mSignal\tmProb\tNorf\tNgenes\tNhit\n";
print
  "$ratio_threshold\t$prob_threshold\t$minrun\t$maxgap1\t$maxgap2\t$avg_run\t",
  "$avg_ratio\t$avg_prob\t$orf_cnt\t$genes_cnt\t$hit_cnt\n";

# print to file:
print OUT "# input file: $infile1\n";
print OUT "# Parameters used:
# --------------------------------------------------------------------------------------------\n";
print OUT
"# Thresholds: log2(signal) = $ratio_threshold, -log10(prob)  = $prob_threshold,
# minrun = $minrun, maxgap1 = $maxgap1, maxgap2 = $maxgap2
# --------------------------------------------------------------------------------------------\n";
print OUT "# SUMMARY\n";
print OUT "# Tsig\tTprob\tminrun\tmaxgap1\tmaxgap2\tmLength\t",
  "mProb\tmSignal\tNprom\tNgenes\tNorf\tNhit\n";
print OUT
"# $ratio_threshold\t$prob_threshold\t$minrun\t$maxgap1\t$maxgap2\t$avg_run\t",
  "$avg_prob\t$avg_ratio\t$prom_cnt\t$genes_cnt\t$orf_cnt\t$hit_cnt\n";
print OUT
"# -------------------------------------------------------------------------------------------------------\n";
print OUT "# Chr\tShoulder\tstart\tpeak\tend\tlength\t#probes\tmSignal\tmProb\t",
  "geneW1\t(aka)\tstatus\tdistW1\tgeneW2\t(aka)\tstatus\tdistW2\t",
  "geneC1\t(aka)\tstatus\tdistC1\tgeneC2\t(aka)\tstatus\tdistC2\torf\n";

for $row (@results) {
	@$row[0] =~ /chr(\d*)/;
	$chrNum = $1;
	$begin  = @$row[4];
	$begin =~ s/ //g;
	$end = @$row[6];
	$end =~ s/ //g;
	$html =
"http://db.yeastgenome.org/cgi-bin/ORFMAP/ORFmap?chr=$chrNum&beg=$begin&end=$end";

	print OUT
"@$row[0]\t@$row[3]\t@$row[4]\t@$row[5]\t@$row[6]\t@$row[7]\t@$row[8]\t@$row[9]\t",
"@$row[10]\t@$row[11]\t@$row[12]\t@$row[13]\t@$row[14]\t@$row[15]\t@$row[16]\t@$row[17]\t",
"@$row[18]\t@$row[19]\t@$row[20]\t@$row[21]\t@$row[22]\t@$row[23]\t@$row[24]\t@$row[25]\t@$row[26]\t@$row[28]\t$html\n";
}

close IN1;
close OUT;
close OUT2;

sub readfile {
	my ($infile1) = @_;

	open( IN1, $infile1 ) || die "can't open $infile1: $!";
	our @chr   = ();
	our @pos   = ();
	our @prob  = ();
	our @ratio = ();
	$i = 0;
	while ( defined( my $line = <IN1> ) ) {
		$line =~ /^(.*?)\t(.*?)\t(.*?)\t(.*?)\n$/;
		$chr = $1;
		if ( substr( $chr, 0, 4 ) eq "mito" ) { $chr = "chr17" }
		$chr[$i]   = $chr;
		$pos[$i]   = $2;
		$prob[$i]  = $3;
		$ratio[$i] = $4;
		$ratio[$i] =~ s/ //g;
		$chr[$i]   =~ s/ //g;
		$i++;
	}
	close IN1;
	return \@chr, \@pos, \@prob, \@ratio;
}

sub maxgap1 {

	my ( $maxgap1, $prob_thresh, $ratio_thresh ) = @_;

	our @chr;
	our @pos;
	our @ratio;
	my $run_flag  = 0;
	my $i         = 0;
	my $chr       = 0;
	my $begin_pos = 0;
	my $begin_idx = 0;
	my $true_run  = 0;
	my @runs1     = ();

	for ( $i = 0 ; $i < scalar(@chr) ; $i++ ) {
		$last_chr = $chr;
		$chr      = $chr[$i];
		$pos      = $pos[$i];
		$prob     = $prob[$i];
		$ratio    = $ratio[$i];

		# $value = $values[$which];

		if ( ( defined $maxgap_end && $pos > $maxgap_end )
			|| $chr ne $last_chr )
		{    # if new probe is on a different chromosome...
			$run_flag = 0;    # there is no more an ongoing run.
			if ( $true_run == 1 )
			{                 # if this run satisfies minimum requirements...

				@temp =
				  ( $last_chr, $begin_idx, $end_idx, $begin_pos, $end_pos );
				push @runs1, [@temp];
				$true_run  = 0;
				$begin_pos = 0;
				$begin_idx = 0;
				undef my $maxgap_end;
			}
		}

		if (   $prob >= $prob_thresh
			&& $ratio >= $ratio_thresh )
		{

			# print"i=$i, pos=$pos, ratio=$ratio\n";
			$end_pos = $pos;    # provisional end of the run is here now.
			$end_idx = $i;
			if ( $run_flag == 0 ) {    # if this is the first probe of a run...
				$maxgap_end = $pos + $maxgap1;
				$begin_pos  = $pos;
				$begin_idx  = $i;
				$run_flag   = 1;
			} elsif ( $pos <= $maxgap_end ) {
				$true_run = 1;    # this is a true run and will be saved later
				$maxgap_end = $pos + $maxgap1;
			}
		}
	}

	#	close IN1;

	# finished while loop, must process last run
	if ( $true_run == 1 ) {    # if this run satisfies minimum requirements...
		@temp = ( $last_chr, $begin_idx, $end_idx, $begin_pos, $end_pos );
		push @runs1, [@temp];
	}
	return @runs1;
}

sub minrun {

	# print "*** start minrun\n";
	my ( $minrun, @runs ) = @_;
	my @runs2 = ();
	for $j ( 0 .. $#runs ) {
		$chr       = $runs[$j][0];
		$begin_idx = $runs[$j][1];
		$end_idx   = $runs[$j][2];
		$begin_pos = $runs[$j][3];
		$end_pos   = $runs[$j][4];

		if ( $end_pos - $begin_pos >= $minrun ) {
			@temp = ( $chr, $begin_idx, $end_idx, $begin_pos, $end_pos );
			push @runs2, [@temp];

			# print join ":", @temp, "\n";
		}
	}
	return @runs2;
}

sub maxgap2 {

	# print"***begin maxgap2...\n";
	my ( $maxgap2, @runs2 ) = @_;
	my @runs3     = ();
	my $chr       = $runs2[0][0];
	my $begin_idx = $runs2[0][1];
	my $end_idx   = $runs2[0][2];
	my $begin_pos = $runs2[0][3];
	my $end_pos   = $runs2[0][4];

	if ( $#runs2 == 0 ) { return $runs2[0]; }

	for my $j ( 1 .. $#runs2 ) {
		$new_chr       = $runs2[$j][0];
		$new_begin_idx = $runs2[$j][1];
		$new_end_idx   = $runs2[$j][2];
		$new_begin_pos = $runs2[$j][3];
		$new_end_pos   = $runs2[$j][4];

	#print "$new_chr $new_begin_idx $new_end_idx $new_begin_pos $new_end_pos\n";

		if ( ( abs( $new_begin_pos - $end_pos ) > $maxgap2 )
			|| $chr ne $last_chr )
		{
			@temp = ( $chr, $begin_idx, $end_idx, $begin_pos, $end_pos );
			push @runs3, [@temp];

			$chr       = $new_chr;
			$begin_idx = $new_begin_idx;
			$end_idx   = $new_end_idx;
			$begin_pos = $new_begin_pos;
			$end_pos   = $new_end_pos;
			next;
		}
		if ( $chr eq $new_chr && abs( $new_begin_pos - $end_pos ) <= $maxgap2 )
		{
			$end_idx = $new_end_idx;
			$end_pos = $new_end_pos;
		}

	}

	# finished for loop, must process last record...
	if ( $#runs2 > 0 ) {
		@temp = ( $chr, $begin_idx, $end_idx, $begin_pos, $end_pos );
		push @runs3, [@temp];

	}
	return @runs3;
}

sub calc_run {

 # re-reads datafile and sums the ratio and probability values of probes
 # within above-threshold runs. returns array of arrays containing these results
 # for each hit.

	(@runs) = @_;

	our @chr;
	our @pos;
	our @ratio;

	if ( scalar @runs == 0 ) { return @runs; }
	our @peaks = ();
	my $j = 0;    # each $j is for one run above threshold
	$run_chr       = $runs[$j][0]; 
	$run_begin_idx = $runs[$j][1];
	$run_end_idx   = $runs[$j][2];
	$run_begin_pos = $runs[$j][3];
	$run_end_pos   = $runs[$j][4];

	$neighb_start =
	  $run_begin_idx - sprintf( "%.0f", ( $window / $probe_length ) / 2 ) - 1;
	if ( $neighb_start < 0 ) { $neighb_start = 0 }
	$neighb_end =
	  $run_end_idx + sprintf( "%.0f", ( $window / $probe_length ) / 2 ) + 1;

# if neighb_start is not on same chromosome as $run_begin_idx, increment it until it is so.
	while ( $chr[$neighb_start] ne $run_chr ) { $neighb_start++ }

# if neighb_end is not on same chroosome as $run_end_idx, decrement it until it is so.

	while ( !defined( $chr[$neighb_end] ) ) {
		$neighb_end--;
	}    # fix for a rare pathological case
	while ( $chr[$neighb_end] ne $run_chr ) { $neighb_end-- }

	# if run is at the very start or very end of the chromosome:
	if ( $neighb_start == $run_begin_idx ) { $run_begin_idx++ }
	if ( $neighb_end == $run_end_idx )     { $run_end_idx-- }

	#print "neighb_start =$neighb_start neighb_end=$neighb_end \n ";

	#for $z ( 0 .. $#runs ) {
	#	$row = $runs[$z];
	#    print "@$row", "\n";
	#}

	$run_flag = 0;

	for ( my $i = 0 ; $i < scalar(@chr) ; $i++ ) {
		$chr = $chr[$i];
		$pos = $pos[$i];

		if ( $chr ne $run_chr && $run_flag == 0 ) {
			next;
		}

		if ( $i < $neighb_start && $run_flag == 0 ) {
			next;
		}

		#		if ( 	$i > $run_begin_idx
		#			&& 	$run_flag == 0
		#			&& chr[$neighb_start] eq $run_chr
		#			) {
		#				# if previous run's array ended after $run_begin_idx
		#				$i= $neighb_start; # resetting the counter here!
		#				 print "resetting counter\n";
		#				$run_flag=1;
		#				next;
		#				}

		if (   $i >= $neighb_start
			&& $i < $neighb_end )
		{

# print"run_begin_idx=$run_begin_idx,run_end_idx=$run_end_idx,i=$i,chr=$chr, pos=$pos\n";

			$run_flag = 1;

			#				$run_flag	= 1;

			# print "not ended... pos{$i}=$pos{$i}\n";
		} elsif (
			( $i = $neighb_end && $run_flag == 1 )

			#				 	||	( $chr ne $run_chr && $run_flag ==1)
		  )
		{
			$run_flag = 0;    # run has ended; will generate summary stats

			# print "*ended... pos{$i}=$pos{$i}\n";

			# for each run, there may be multiple peaks;
			# get idx values of start, peak, and end of each peak.

			( $startRef, $peakRef, $endRef, $peakFlagRef ) =

			  &getpeaks(
				$run_begin_idx ,
				$run_end_idx ,
				$peak_sep, $window, $neighb_start, $neighb_end
			  );
			@start = @$startRef;
			@peak  = @$peakRef;
			@end   = @$endRef;
			@peakFlag = @$peakFlagRef;

			#  print "...STARTs:" , join " ", @start, "\n";
			#  print "...PEAKS:" , join " ", @peak, "\n";
			#  print "...ENDs:" , join " ", @end, "\n";

			#				 $scalar_peak = scalar(@peak);
			#				 print "scalar_peak=$scalar_peak\n";

			for ( $p = 0 ; $p < scalar(@peak) ; $p++ ) {
				$sum_signal = 0;
				$sum_prob   = 0;
				$count      = 0;

	  #		  print " p=$p, start[$p]=$start[$p] end[$p]=$end[$p]  count=$count\n";
	  #		  print "%%% p=$p\n";
				for ( $k = $start[$p] ; $k <= $end[$p] ; $k++ ) {

# print "*** p=$p, k=$k, start[$p]=$start[$p] end[$p]=$end[$p] signalk= $signal{$k} count=$count\n";
					$sum_signal = $sum_signal + $ratio[$k];
					$sum_prob   = $sum_prob + $prob[$k];
					$count      = $count + 1;

				}    #end of $k loop for probes within a peak
				if ( $count < 3 ) { next }

				# print "count=$count, p=$p\n";
				$mean_signal = sprintf( "%.3f", $sum_signal / $count );
				$mean_prob   = sprintf( "%.3f", $sum_prob / $count );
				@temp        = (
					$chr, $start[$p], $end[$p],    #chr, start_idx, end_idx,
					$peakFlag[$p],			#whether a Peak or Shoulder
					$pos[ $start[$p] ],                      # start pos,
					$pos[ $peak[$p] ],                       # peak pos
					$pos[ $end[$p] ] + $probe_length - 1,    # end pos,
					$pos[ $end[$p] ] - $pos[ $start[$p] ] +
					  $probe_length,                         #length
					$count,          # numb of probes in peak
					$mean_signal,    # mean signal
					$mean_prob       # mean prob value
				);

				#					 print join " ", @temp, "\n";
				push @peaks, [@temp];

			}    #end of $p loop for peaks within a run.

			if ( $j == $#runs ) {
				last;    # if end of list of runs, stop reading file
			}

			$j++;        # process next run
			$run_chr       = $runs[$j][0];
			$run_begin_idx = $runs[$j][1];
			$run_end_idx   = $runs[$j][2];
			$run_begin_pos = $runs[$j][3];
			$run_end_pos   = $runs[$j][4];

			#				%pos           = ();
			#				%signal        = ();
			$neighb_start = $run_begin_idx -
			  sprintf( "%.0f", ( $window / $probe_length ) / 2 );
			if ( $neighb_start < 0 ) { $neighb_start = 0 }
			$neighb_end =
			  $run_end_idx + sprintf( "%.0f", ( $window / $probe_length ) / 2 );

# if neighb_start is not on same chromosome as $run_begin_idx, increment it until it is so.
			while ( $chr[$neighb_start] ne $run_chr ) { ++$neighb_start }

# if $neighb_end beyond end of last chromosome, decrement it until it is defined
			while ( !defined $chr[$neighb_end] ) { --$neighb_end }

# if neighb_end is not on same chromosome as $run_end_idx, decrement it until it is so.

			while ( $chr[$neighb_end] ne $run_chr ) { --$neighb_end }

# print "+++$runs[$j][0]|$runs[$j][1]|$runs[$j][2]|$runs[$j][3]|$runs[$j][3]|$runs[$j][4]\n";
			if ( $neighb_start == $run_begin_idx ) { $run_begin_idx++ }
			if ( $neighb_end == $run_end_idx )     { $run_end_idx-- }

		}    # end processing one run
		     #}    # end processing one chromosome
	}    # end reading arrays
	return @peaks;
}

sub locate_run {
	my (@runs) = @_;

	my $dbh = DBI->connect( "dbi:mysql:chipchip", "cbruce", "yalejohn" )
	  or die "Can't connect to mysql database: $DBI::errstr\n";

	for my $j ( 0 .. $#runs ) {

		$chr2 = $runs[$j][0];

		$chr2 =~ s/ //g;
		$chr           = substr( $chr2, 3 );
		$middle_of_run = $runs[$j][5];
		$start_of_run  = $runs[$j][4];
		$end_of_run    = $runs[$j][6];

		my $stmt1 = $dbh->prepare(
			qq{
    			SELECT name, aka, status, start, end, upstr_start, upstr_end
          		  FROM yeast g
	 		 WHERE upstr_end <= $middle_of_run  
#           		   AND start+50 > $middle_of_run
                           AND end > $middle_of_run
			   AND g.strand = 'W' 
		           AND g.vvv= $chr
           		   AND g.type = 'ORF'
           		   LIMIT 1
                   	}
		) or die "can't prepare stmt: $DBI::errstr";

		my $stmt2 = $dbh->prepare(
			qq{
               		SELECT name, aka, status, start, end, upstr_start, upstr_end
                 	  FROM yeast g
#                	 WHERE start-50 <= $middle_of_run
                         WHERE end <= $middle_of_run  
 	        	   AND upstr_end > $middle_of_run 
			   AND g.strand = 'C'  
                  	   AND g.vvv= $chr
                  	   AND g.type = 'ORF'
                  	   LIMIT 1
                 	}
		) or die "can't prepare stmt: $DBI::errstr";

#		Statment to find # hits that are between ORF start +50 and end of ORF, on the W strand
		my $stmt3 = $dbh->prepare(
			qq{
                        SELECT count(*)
                          FROM yeast g
                         WHERE start+50 <  $middle_of_run 
    			   AND end > $middle_of_run
	   		   AND g.strand = 'W'
                           AND g.vvv= $chr
                           AND g.type = 'ORF'
                        }
		) or die "can't prepare stmt: $DBI::errstr";

#		Statment to find # hits that are between ORF start +50 and end of ORF, on the C strand
		my $stmt4 = $dbh->prepare(
			qq {
                        SELECT count(*)
                          FROM yeast g
                         WHERE  start-50 > $middle_of_run
                           AND end < $middle_of_run
			   AND  g.strand = 'C'
                           AND g.vvv= $chr
                           AND g.type = 'ORF'
                        }
		) or die "can't prepare stmt: $DBI::errstr";

		my $stmt5 = $dbh->prepare(
			qq{
                        SELECT name
                          FROM yeast g
                         WHERE (
                            (upstr_end <=  $start_of_run AND end > $end_of_run )
                            OR 
                            (upstr_end >  $start_of_run AND end <= $end_of_run )
                            OR
                            (upstr_end between  $start_of_run AND $end_of_run )
                            OR
                            (end between  $start_of_run AND $end_of_run )
                            )
                           AND g.strand = 'W'
                           AND g.vvv= $chr
                           AND g.type = 'ORF'
                        }
		  )
		  or die "can't prepare stmt: $DBI::errstr";

		my $stmt6 = $dbh->prepare(
			qq{
                        SELECT name
                          FROM yeast g
                         WHERE (
                            (upstr_end <=  $start_of_run AND end > $end_of_run )
                            OR 
                            (upstr_end >  $start_of_run AND end <= $end_of_run )
                            OR
                            (upstr_end between  $start_of_run AND $end_of_run )
                            OR
                            (end between  $start_of_run AND $end_of_run )
                            )
                           AND g.strand = 'C'
                           AND g.vvv= $chr
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

#		Get genes on W strand where hit peak is between end of upstream gene and ORF start+50.		
		@row1 = $stmt1->fetchrow_array();
		warn "Data fetching terminated early by error: $DBI::errstr\n"
		  if $DBI::err;
		if ( scalar(@row1) != 0 ) {
			(
				$geneW1, $akaW1, $statusW1, $startW1, $endW1, $upstr_startW1,
				$upstr_endW1
			  )
			  = @row1;

			$distW1 = $middle_of_run - $startW1;
			if    ( $statusW1 eq "Dubious" )                { $statusW1 = "D" }
			elsif ( $statusW1 eq "Uncharacterized" )        { $statusW1 = "u" }
			elsif ( $statusW1 eq "Verified" )               { $statusW1 = "V" }
			elsif ( $statusW1 eq "Verified|silenced_gene" ) { $statusW1 = "s" }
			else { $statusW1 = "" }

		} else {
			(
				$geneW1, $akaW1, $statusW1, $startW1, $endW1, $upstr_startW1,
				$upstr_endW1
			  )
			  = ( "", "", "", "", "", "", "" );
			$distW1 = "";
		}

#		Get genes on C strand where hit peak is between end of upstream gene and ORF start+50.
		@row2 = $stmt2->fetchrow_array();
		warn "Data fetching terminated early by error: $DBI::errstr\n"
		  if $DBI::err;
		if ( scalar(@row2) != 0 ) {
			(
				$geneC1, $akaC1, $statusC1, $startC1, $endC1, $upstr_startC1,
				$upstr_endC1
			  )
			  = @row2;

			$distC1 = $startC1 - $middle_of_run;
			if    ( $statusC1 eq "Dubious" )                { $statusC1 = "D" }
			elsif ( $statusC1 eq "Uncharacterized" )        { $statusC1 = "u" }
			elsif ( $statusC1 eq "Verified" )               { $statusC1 = "V" }
			elsif ( $statusC1 eq "Verified|silenced_gene" ) { $statusC1 = "s" }
			else { $statusC1 = "" }
		} else {
			(
				$geneC1, $akaC1, $statusC1, $startC1, $endC1, $upstr_startC1,
				$upstr_endC1
			  )
			  = ( "", "", "", "", "", "", "" );
			$distC1 = "";

		}

#		Statement to get the gene downstream from the hit for the W strand.
		my $stmt7 = $dbh->prepare(
			qq{
   				SELECT name, aka, status, start
          		  FROM yeast g
	 		     WHERE upstr_end = $endW1  
			   	   AND g.strand = 'W'
		           AND g.vvv= $chr
           		   AND g.type = 'ORF'
                 LIMIT 1
		  	}
		  )
		  or die "can't prepare stmt: $DBI::errstr";

		if ( $endW1 ne "" ) {
			$stmt7->execute()
			  or die "Can't execute SQL statement: $DBI::errstr\n";
		}

#		Statement to get the gene downstream from the hit for the C strand.
		my $stmt8 = $dbh->prepare(
			qq{
   				SELECT name, aka, status, start 
          		  FROM yeast g
	 		     WHERE upstr_end = $endC1  
			   	   AND g.strand = 'C'
		           AND g.vvv= $chr
           		   AND g.type = 'ORF'
                 LIMIT 1
		  	}
		  )
		  or die "can't prepare stmt: $DBI::errstr";
#		$stmt8->bind_param( 1, $endC1 );
#		$stmt8->bind_param( 2, 'C' );
#		$stmt8->bind_param( 3, $chr );

		if ( $endC1 ne "" ) {
			$stmt8->execute()
			  or die "Can't execute SQL statement 8: $DBI::errstr\n";
		}

		my @row3 = $stmt3->fetchrow_array();
		warn "Data fetching terminated early by error: $DBI::errstr\n"
		  if $DBI::err;

		my @row4 = $stmt4->fetchrow_array();
		warn "Data fetching terminated early by error: $DBI::errstr\n"
		  if $DBI::err;

		my @genes = ();

		while ( $row5 = $stmt5->fetchrow_arrayref ) {
			if ( defined $row5 ) { 	push @genes, @$row5[0] }
		}
		warn "Data fetching terminated early by error: $DBI::errstr\n"
		  if $DBI::err;

		while ( $row6 = $stmt6->fetchrow_arrayref ) {
			if ( defined $row6 ) { push @genes, @$row6[0] }
		}
		warn "Data fetching terminated early by error: $DBI::errstr\n"
		  if $DBI::err;

#		Get downstream gene for W strand 
		if ( $endW1 ne "" ) {
			@row7 = $stmt7->fetchrow_array();
			warn "Data fetching terminated early by error: $DBI::errstr\n"
			  if $DBI::err;
			if ( scalar(@row7) != 0 ) {
				( $geneW2, $akaW2, $statusW2, $startW2 ) = @row7;
				$distW2 = $middle_of_run - $startW2;
				if    ( $statusW2 eq "Dubious" )         { $statusW2 = "D" }
				elsif ( $statusW2 eq "Uncharacterized" ) { $statusW2 = "u" }
				elsif ( $statusW2 eq "Verified" )        { $statusW2 = "V" }
				elsif ( $statusW2 eq "Verified|silenced_gene" ) {$statusW2 = "s";} 
				else {$statusW2 = "";}
			}
		} else {
			( $geneW2, $akaW2, $statusW2, $startW2, $distW2 ) =
			  ( "", "", "", "", "" );
		}

# are there any ORFs on the C stand between genes W1 and W2?
# if yes, geneW2 will not be reported.
#
#         hit     W1              W2
#         |      ====>           ====>
#---------|----------------------------
#  <===   |              <===
#
#
if (!defined($startW2)) {print "startW2 not defined: geneW1= $geneW1 \n";}
if (!defined($endW1)) {print "endW1 not defined: geneW1= $geneW1 \n";}
		my $stmt9 = $dbh->prepare(
			qq{
   				SELECT count(*)
          		  FROM yeast g
	 		     WHERE start < $startW2
			   	   AND end > $endW1
			   	   AND strand = 'C' 
		           AND vvv = $chr
           		   AND g.type = 'ORF'
		  	}
		  )
		  or die "can't prepare stmt: $DBI::errstr";

		if ( $endW1 ne "" && $startW2 ne "" ) {
			$stmt9->execute()
			  or die "Can't execute SQL statement 9: $DBI::errstr\n";

			@row9 = $stmt9->fetchrow_array();
			warn "Data fetching terminated early by error: $DBI::errstr\n"
			  if $DBI::err;
			if ( $row9[0] > 0 ) {
				( $geneW2, $akaW2, $statusW2, $startW2, $distW2 ) =
				  ( "", "", "", "", "" );
			}
		}
		
# are there any ORFs on the opposite strand between hit and gene W1?
# if yes, geneW1 will not be reported.
#
#         hit         W1      W2
#         |          ====>   ====>
#---------|-------------------------
#  <===   | <====          
#

		my $stmt11 = $dbh->prepare(
			qq{
   				SELECT count(*)
          		  FROM yeast g
	 		     WHERE start < $startW1
			   	   AND end > $middle_of_run
			   	   AND strand = 'C' 
		           AND vvv = $chr
           		   AND g.type = 'ORF'
		  	}
		  )
		  or die "can't prepare stmt: $DBI::errstr";

		if ( $startW1 ne "" ) {
			$stmt11->execute()
			  or die "Can't execute SQL statement 11: $DBI::errstr\n";

			@row11 = $stmt11->fetchrow_array();
			warn "Data fetching terminated early by error: $DBI::errstr\n"
			  if $DBI::err;
			if ( $row11[0] > 0 ) {
				( $geneW1, $akaW1, $statusW1, $startW1, $distW1 ) =
				  ( "", "", "", "", "" );
			}
		}

#		Get downstream gene for C strand 
		if ( $endC1 ne "" ) {
			@row8 = $stmt8->fetchrow_array();
			warn "Data fetching terminated early by error: $DBI::errstr\n"
			  if $DBI::err;
			if ( scalar(@row8) != 0 ) {
				( $geneC2, $akaC2, $statusC2, $startC2 ) = @row8;

				$distC2 = $startC2 - $middle_of_run;
				if    ( $statusC2 eq "Dubious" )         { $statusC2 = "D" }
				elsif ( $statusC2 eq "Uncharacterized" ) { $statusC2 = "u" }
				elsif ( $statusC2 eq "Verified" )        { $statusC2 = "V" }
				elsif ( $statusC2 eq "Verified|silenced_gene" ) {
					$statusC2 = "s";
				} else {
					$statusC2 = "";
				}
			} else {
				( $geneC2, $akaC2, $statusC2, $startC2 ) = ( "", "", "", "" );
				$distC2 = "";
			}
		} else {
			( $geneC2, $akaC2, $statusC2, $startC2, $distC2 ) =
			  ( "", "", "", "", "" );
		}

# are there any ORFs on the W strand between geneC1 and gene C2?
# if yes, geneC2 will not be reported.
		my $stmt10 = $dbh->prepare(
			qq{
   				SELECT count(*) 
          		  FROM yeast g
	 		     WHERE start > $startC2 
			   	   AND end < $endC1
			   	   AND strand = 'W' 
		           AND vvv = $chr
           		   AND g.type = 'ORF'
		  	}
		  )
		  or die "can't prepare stmt: $DBI::errstr";

		if ( $endC1 ne "" && $startC2 ne "" ) {
			$stmt10->execute()
			  or die "Can't execute SQL statement 9: $DBI::errstr\n";
			@row10 = $stmt10->fetchrow_array();
			warn
			  "Data fetching terminated early by error: $DBI::errstr\n" 
			  if $DBI::err;
			if ( $row10[0] > 0 ) {
				( $geneC2, $akaC2, $statusC2, $startC2, $distC2 ) =
				  ( "", "", "", "", "" );
			}
		}

		if ( $row3[0] != 0 or $row4[0] != 0 ) { $orf = 1; }
		else { $orf = 0 }

		$genes = join ",", @genes;

# are there any ORFs on the opposite strand between hit and gene C1?
# if yes, geneC1 will not be reported.
#
#                      hit   
#                 ====> |   ====>
#-----------------------|---------
#  <===   <====         |
#   C2      C1
#

		my $stmt12 = $dbh->prepare(
			qq{
   				SELECT count(*)
          		  FROM yeast g
	 		     WHERE start > $startC1
			   	   AND end < $middle_of_run
			   	   AND strand = 'W' 
		           AND vvv = $chr
           		   AND g.type = 'ORF'
		  	}
		  )
		  or die "can't prepare stmt: $DBI::errstr";

		if ( $startC1 ne "" ) {
			$stmt12->execute()
			  or die "Can't execute SQL statement 11: $DBI::errstr\n";

			@row12 = $stmt12->fetchrow_array();
			warn "Data fetching terminated early by error: $DBI::errstr\n"
			  if $DBI::err;
			if ( $row12[0] > 0 ) {
				( $geneC1, $akaC1, $statusC1, $startC1, $distC1 ) =
				  ( "", "", "", "", "" );
                                ( $geneC2, $akaC2, $statusC2, $startC2, $distC2 ) =
                                  ( "", "", "", "", "" );
			}
		}

# are there any ORFs on the W stand between genes C1 and C2?
# if yes, geneC2 will not be reported.
#
#                      hit
#        ====>          |   ====>
#-----------------------|---------
#  <===         <===    |
#   C2            C1
#

                my $stmt13 = $dbh->prepare(
                        qq{
                                SELECT count(*)
                          FROM yeast g
                             WHERE start > $startC2
                                   AND end < $endC1
                                   AND strand = 'W'
                           AND vvv = $chr
                           AND g.type = 'ORF'
                        }
                  )
                  or die "can't prepare stmt: $DBI::errstr";

                if ( $startC1 ne "" &&  $startC2 ne "" ) {
                        $stmt13->execute()
                          or die "Can't execute SQL statement 13: $DBI::errstr\n";

                        @row13 = $stmt13->fetchrow_array();
                        warn "Data fetching terminated early by error: $DBI::errstr\n"
                          if $DBI::err;
                        if ( $row13[0] > 0 ) {
                                ( $geneC2, $akaC2, $statusC2, $startC2, $distC2 ) =
                                  ( "", "", "", "", "" );
                        }
                }




#requirement: if middle of hit not downstream of ORF start+50nt.,
#	and it is not a dubious ORF, don't report next downstream ORF
		if ( $startW1 ne "" ) {
			if ( $middle_of_run < $startW1+50 && $statusW1 ne 'D' ) {
       #                if ( $middle_of_run < $endW1 && $statusW1 ne 'D' ) {
				$geneW2   = "";
				$akaW2    = "";
				$statusW2 = "";
				$distW2   = "";
				}
		}

		#requirement: same as above.
		if ( $startC1 ne "" ) {
			if ( $middle_of_run > $startC1+50 && $statusC1 ne 'D' ) {
       #                if ( $middle_of_run > $endC1 && $statusC1 ne 'D' ) {
				$geneC2   = "";
				$akaC2    = "";
				$statusC2 = "";
				$distC2   = "";
			}
		}

 #requirement: If there is an intervening ORF between a hit and a potential gene
 #			on the other strand, the ORF is assumed to interfere with transcription
 #			of this potential gene.
 #
 #                           hit     yes
 #                    ---->  ||     ----->
 #           ----------------||-----------------
 #           <----
 #              no
 #
 #			That is,
 #			if hit is completely between $startW1 and $upstr_endW1,
 #				and startC1 < $upstrW1, don't report geneC1 and geneC2.
 #			Similarly, if hit is completely between $startC1 and $upstr_endC1,
 #				and startW1 > $upstrC1, don't report geneW1 and geneW2.

		if ( $startC1 ne "" && $upstr_endW1 ne "" && $upstr_startW1 ne "" ) {
			if (   $start_of_run > $upstr_endW1
				&& $startC1 < $upstr_startW1 )
			{
				$geneC1   = "";
				$akaC1    = "";
				$startC1  = "";
				$statusC1 = "";
				$distC1   = "";
				$geneC2   = "";
				$akaC2    = "";
				$statusC2 = "";
				$distC2   = "";
			}
		}

		if ( $upstr_startC1 ne "" && $upstr_endC1 ne "" && $startW1 ne "" ) {
			if (   $end_of_run < $upstr_endC1
				&& $startW1 > $upstr_startC1 )
			{
				$geneW1   = "";
				$akaW1    = "";
				$startW1  = "";
				$statusW1 = "";
				$distW1   = "";
				$geneW2   = "";
				$akaW2    = "";
				$statusW2 = "";
				$distW2   = "";
			}
		}

		push @{ $runs[$j] }, ( $geneW1, $akaW1, $statusW1, $distW1 );
		push @{ $runs[$j] }, ( $geneW2, $akaW2, $statusW2, $distW2 );
		push @{ $runs[$j] }, ( $geneC1, $akaC1, $statusC1, $distC1 );
		push @{ $runs[$j] }, ( $geneC2, $akaC2, $statusC2, $distC2 );
		push @{ $runs[$j] }, ( $genes,  $orf );
	}
	@results = @runs;

	# Now, disconnect from the database
	$dbh->disconnect or warn "Disconnection failed: $DBI::errstr\n";

	return @results;
}

sub DB_output {
	my (
		$prob,  $minrun, $maxgap1,  $maxgap2, $mRun,
		$mProb, $mRatio, $prob_cnt, $hit_cnt
	  )
	  = @_;

	my $dbh = DBI->connect( "dbi:mysql:chipchip", "cbruce", "yalejohn" )
	  or die "Can't connect to mysql database: $DBI::errstr\n";

	my $stmt1 = $dbh->prepare(
		qq{
    	insert into params values (
		'', $prob, $minrun, $maxgap1, $maxgap2, $mRun, $mProb, $mRatio, $prob_cnt, $hit_cnt)
        }
	  )
	  or die "can't prepare stmt: $DBI::errstr";

	$stmt1->execute() or die "Can't execute SQL statement: $DBI::errstr\n";

	# Now, disconnect from the database
	$dbh->disconnect or warn "Disconnection failed: $DBI::errstr\n";
	return @results;
}

sub linreg {
	my ( $Aref, $Bref ) = @_;
	my @A = @$Aref;
	my @B = @$Bref;

	#print join " ", @A, "\n";
	#print join " ", @B, "\n";

	( $sumx, $sumy, $sumx2, $sumy2, $sumxy ) = ( 0, 0, 0, 0, 0 );
	my $n = scalar(@A);    #print "n= $n\n";
	my $i;
	for ( $i = 0 ; $i < scalar(@A) ; $i++ ) {
		$sumx  = $sumx + $A[$i];
		$sumy  = $sumy + $B[$i];
		$sumx2 = $sumx2 + $A[$i]**2;
		$sumy2 = $sumy2 + $B[$i]**2;
		$sumxy = $sumxy + $A[$i] * $B[$i];
	}
	my $m =
	  ( $n * $sumxy - $sumx * $sumy ) / ( $n * $sumx2 - $sumx * $sumx )
	  ;                    # compute slope
	                       #$m = sprintf "%.3f", $m;
	my $b =
	  ( $sumy * $sumx2 - $sumx * $sumxy ) / ( $n * $sumx2 - $sumx * $sumx )
	  ;                    # compute y-intercept
	my $r = ( $sumxy - $sumx * $sumy / $n ) /  # compute correlation coefficient
	  sqrt( ( $sumx2 - $sumx * $sumx / $n ) * ( $sumy2 - $sumy * $sumy / $n ) );

	#	print "m,b,r = $m, $b, $r \n";

	return ( $m, $b, $r );
}

sub getpeaks {

# $startidx: idx# of peak start in original data file
# $endidx: idx# of peak start in original data file
# $peak_sep: minimum separtion between two peaks
# $window: window size used to calculate average slope for locating the peak top
# $Bref:  reference to positions hash (smallest key <= $startidx, largest key >= $endidx)
# $Cref:  reference to signal values hash

	my ( $startidx, $endidx, $peak_sep, $window, $neighb_start, $neighb_end ) =
	  @_;

	our @pos;
	our @ratio;

#if($startidx>91120){
#	print "startidx=$startidx, endidx=$endidx, neighb_start=$neighb_start, neighb_end=$neighb_end \n";
#}
#if (!defined $pos[$startidx]){$startidx =$startidx+1}		# a crude fix!
	if ( $pos[$startidx] == 1 ) { $startidx++ }

	#if ($neighb_start==189963){
	#			print "..chr: ($chr[$startidx])/ pos=($pos[$startidx] )\n";
	#			print "neighb_start=$neighb_start, neighb_end=$neighb_end\n";
	#			print "startidx=$startidx, endidx=$startidx\n";
	#			print "left_idx{$i} = $left_idx{$i} right_idx{$i}= $right_idx{$i}\n";
	#				}

	my $i;
	for ( $i = $startidx - 1 ; $i <= $endidx ; $i++ ) {
		if ( $neighb_start > $startidx ) { $startidx++ }    # a crude fix!

		#loop to calculate slope at idx=$i

# print "startidx= $startidx i = $i neighb_start= $neighb_start, lastj= $lastj\n";
# print "pos{$i}= $pos{$i}\n";

		# define $left_idx and $right_idx, between which average slopes
		# will be calculated.  The distance between the center point
		# and the left_idx or right_idx positions is to be less than
		# the window size.

		for ( $j = $neighb_start ; $j <= $i ; $j++ ) {

			#	print "i=$i, j=$j\n";

			#
			if ( $pos[$i] - $pos[$j] > $window / 2 ) {
				next;
			} else {
				$left_idx{$i} = $j;

				#print "left_idx{$i}=$left_idx{$i}\n";
				last;
			}
		}
		for ( $j = $neighb_end ; $j >= $i ; $j-- ) {

			#	print "i=$i, j=$j\n";
			if ( $pos[$j] - $pos[$i] > $window / 2 ) { next }
			else {
				$right_idx{$i} = $j;

				#				print "right_idx{$i}=$right_idx{$i}\n";
				last;
			}
		}

		# calculate slopes:
		@tempB = ();
		@tempC = ();

		if ( $i < $neighb_start ) { $left_idx{$i} = $neighb_start } # crude fix!

		for ( $j = $left_idx{$i} ; $j <= $right_idx{$i} ; $j++ ) {
			push @tempB, $pos[$j];
			push @tempC, $ratio[$j];
		}

		( $m, undef, undef ) = linreg( \@tempB, \@tempC );
		$slope{$i} = $m;

# print "i=$i left_idx{$i}=$left_idx{$i} right_idx{$i}=$right_idx{$i} slope=$slope{$i};\n";
	}

	# find idx of highest value between startidx and endidx

	for ( $i = $startidx ; $i <= $endidx ; $i++ ) {
		$innerval{$i} = $ratio[$i];
	}
	@sorted_idx = sort { $innerval{$b} <=> $innerval{$a} } keys %innerval;
	%innerval   = ();

	# print join ":", @sorted_idx, "\n";

	$tallest_idx = $sorted_idx[0];

	#print "largestvalue = $tallest_idx\n";

 # will identify start, peak and end indexes for each peak in this range.
 # A start is either a valley preceding a peak, or if there is no such valley,
 # then $startidx, the start of the range.
 # An end is either a valley following a peak, or if there is no such valley,
 # then $endidx, the end of the range.
 # Special case: if index of highest value ($peak_idx) is $startidx,
 # then peak is $startidx, start is $startidx-1.  Similarly, if index of highest
 # value is $endidx, then peak is $endidx, end is $endidx+1.

	@peak      = ();
	@start     = ();
	@end       = ();
        @peakFlag  = ();
	%dist      = ();
	$p         = 0;            # index for peaks
	$start[$p] = $startidx;    # start position of first peak
	$pflag     = 0;            # flag that we are within a peak

	#special case: if highest value is at $startx, that's called a peak
	if ( $tallest_idx == $startidx ) {
		$peak[$p]  = $startidx;
		$start[$p] = $startidx;
#		$peakFlag[$p] = "P";
		
		 #print "******\n";
		 #print "special: start[$p]= $start[$p] peak[$p]=$peak[$p]\n";
		 #print "dist{$p} $dist{$p}\n";
	} 
#       else {
#		$peakFlag[$p] = "S";
#	}
	for ( $i = $startidx ; $i <= $endidx ; $i++ ) {

		# if we are at a valley:
		if ( $slope{ $i - 1 } < 0 && $slope{$i} >= 0 ) {
			if ( !defined $peak[$p] )
			{   # if a valley comes before 1st peak, $startidx called a peak too
				$peak[$p]  = $startidx;
				$start[$p] = $startidx;
				$end[$p]   = $i + 1;

	   #print "valley before first peak: peak[$p]=$peak[$p] end[$p]=$end[$p]\n";
				$p++;
				$start[$p] = $i - 1;
#				$peakFlag[$p] = "S";       # peak flagged by default as a Shoulder 
#	                      				 # unless if corrected later to be a Peak

				#print"start[$p]=$start[$p]\n";

			} else {    # if valley follows a peak
				$end[$p] = $i - 1;

				#print "valley after peak: end[$p]=$end[$p]\n";

 # end-peak pos will be nearest right valley,
 # or position that is at distance $peak_width/2, whichever is less.
 # print "end[$p]>right_idx{ $peak[$p] } $end[$p] > $right_idx{ $peak[$p] } \n";

				if ( $end[$p] > $peak[$p] + $peak_width /(2*$probe_length) ) {
					$end[$p] =
					  sprintf( "%.0f", $peak[$p] + $peak_width /(2* $probe_length) +1 );

			    #print "endp changed, end[$p]=$end[$p] pos:$pos[$start[$p]] \n";
				}
				$p++;
				$start[$p] = $i - 1;
#				$peakFlag[$p] = "S";       # peak flagged by default as a Shoulder 
#	                                       # unless if corrected later to be a Peak

				 #print".start[$p]=$start[$p]\n";
				 #print "valley after peak: start[$p] = $start[$p]\n";

			}
		}

		# if we are at a peak:
		if ( $slope{ $i - 1 } > 0 && $slope{$i} <= 0 ) {

			# if peak is not preceded by a valley, its start is $startidx
			if ( !defined $start[$p] ) {
				$start[$p] = $startidx;

				 #print "no left valley start[$p]= $start[$p]\n";

			}
			$peak[$p] = $i - 1;
			
			$pflag = 1;
		}
	}

	#special case: if highest value is at $endidx, that's called a peak
	if ( $tallest_idx == $endidx ) {
		if ( !defined $start[$p] ) {
			$start[$p] = $startidx;

			 #print "special: start[$p]= start[$p]\n";
		}
		$peak[$p] = $endidx;
		$end[$p]  = $endidx;
		$dist{$p} = 0;
#		$peakFlag[$p] = "P";

		#print "tallest_idx = $tallest_idx endidx=$endidx\n";
		#print "special: peak[$p]= $peak[$p] end[$p]= $end[$p]\n";
		#print "dist{$p} $dist{$p}\n";
	}

	#if a valley never followed last peak, $end is end of range
	if ( !defined( $end[$p] ) ) { $end[$p] = $endidx; }

	#if peak never followed last start site, endof range will be a peak,
	# and following it will be the end:
	if ( !defined( $peak[$p] ) ) {
		$peak[$p] = $endidx; #print "defined peak-p=$peak[$p], ";
		$end[$p]  = $endidx; #print "end-p=$end[$p],";
		$dist{$p} = abs( $tallest_idx - $endidx ); #print "dist p=$dist{$p}\n";
	}

	 #print "....\n";
	# $num_peaks=scalar(@peak);print"num_peaks=$num_peaks\n";
	
	
# Peak resolution if there are more than one peaks:
# If tallest peak is first or last peak of range, assume it is symmetrical, 
# use its outer flank to estimate what the inner flank should be like.  
# Extend calculated inner flank to see where it cuts next peak position.  
# If next peak's height affected by more 10% call it a shoulder, if less, a "Peak".


#for ($p=0;$p<scalar(@peak); $p++){
# print "\n chr=$chr, peak=pos[$peak[$p]] indexes: start[$p]=$start[$p], peak[$p]=$peak[$p], end[$p]=$end[$p]\n";
#}

if (scalar(@peak)==2)
{ #print "Two peaks in run.\n";
  #print "Left peak at chr=$chr, (idx=$peak[0]) position $pos[$peak[0]] \n";
  #print "Right peak at chr=$chr, (idx=$peak[1]) position $pos[$peak[1]] \n";
  $left_outer_10pct_idx= $start[0]+int(0.1* ($peak[0]-$start[0]) );
  $left_outer_90pct_idx= $start[0]+int(0.9*($peak[0]-$start[0]));
  $right_outer_10pct_idx= $peak[1]+int(0.1* ($end[1]-$peak[1]) );
  $right_outer_90pct_idx= $peak[1]+int(0.9*($end[1]-$peak[1]));
  if (($pos[$left_outer_90pct_idx] - $pos[$left_outer_10pct_idx])==0) {$left_peak_slope=0.000001;}
  else {$left_peak_slope = ($ratio[$left_outer_90pct_idx]-$ratio[$left_outer_10pct_idx])/($pos[$left_outer_90pct_idx] - $pos[$left_outer_10pct_idx]);}
  if (($pos[$right_outer_90pct_idx] - $pos[$right_outer_10pct_idx])==0) {$right_peak_slope=0.000001;}
  else{$right_peak_slope = ($ratio[$right_outer_90pct_idx]-$ratio[$right_outer_10pct_idx])/($pos[$right_outer_90pct_idx] - $pos[$right_outer_10pct_idx]);}
  $left_inside_10pct_idx =  $peak[0] + ($peak[0]-$left_outer_10pct_idx); 
  $right_inside_10pct_idx =  $peak[1] - ($right_outer_10pct_idx-$peak[1]);
  $intercept_below_left_peak =  $ratio[$right_inside_10pct_idx] - ($pos[$right_inside_10pct_idx] - $pos[$peak[0]]) * (-$right_peak_slope);
  $intercept_below_right_peak =  $ratio[$left_inside_10pct_idx] + ($pos[$peak[1]] - $pos[$left_inside_10pct_idx]) * (-$left_peak_slope);
  if ($intercept_below_left_peak< 0.1* $ratio[$peak[0]]) 
  {$peakFlag[0]="P"; } 
  else {$peakFlag[0]="M"; }
  if ($intercept_below_right_peak< 0.1* $ratio[$peak[1]]) 
  {$peakFlag[1]="P"; } 
  else {$peakFlag[1]="M"; }
  #print " $peakFlag[0]: height of right peak below left peak= $intercept_below_left_peak, left peak height = $ratio[$peak[0]]\n"; 
  #print " $peakFlag[1]: height of left peak below right peak= $intercept_below_right_peak, right peak height = $ratio[$peak[1]]\n";
}


	#will shrink peaks widths so they are no more than $peak_sep
	for ( $pk = 0 ; $pk < scalar(@peak) ; $pk++ ) {
		$ratio = ( $pos[ $end[$pk] ] - $pos[ $start[$pk] ] ) / $peak_width;

#print "start[$p]=$start[$p] peak[$p]=$peak[$p] end[$p]=$end[$p] dist{$p}= $dist{$p}\n";
	#print"RATIO=$ratio\n";
		if ( $ratio > 1 ) {

#print"==> pos{start[$pk]}=$pos{$start[$pk]} pos{peak[$pk]}=$pos{$peak[$pk]} pos{end[$pk]}=$pos{$end[$pk]}\n";
			$start[$pk] = sprintf( "%.0f",
				$peak[$pk] - ( $peak[$pk] - $start[$pk] ) / $ratio );
			$end[$pk] = sprintf( "%.0f",
				$peak[$pk] + ( $end[$pk] - $peak[$pk] ) / $ratio );

# print"==> pos{start[$pk]}=$pos{$start[$pk]} pos{peak[$pk]}=$pos{$peak[$pk]} pos{end[$pk]}=$pos{$end[$pk]}\n";
		}
	}
	# next, we will combine any peaks
	# that are separated from it by less than the resolution limit.
	# The nearest neighbor within resolution limit was previously
	# identified as $left_idx and right_idx values while setting the
	# boundaries of the window used to calculate the slope of the
	# regression line.

	my @final_peaks  = $peak[0];
	my @final_starts = $start[0];
	my @final_ends   = $end[0];
	my $last_peak    = $peak[0];
	my $last_start   = $start[0];
	my $last_end     = $end[0];

	if ( scalar(@peak) > 1 ) {
		for ( $i = 0 ; $i < scalar(@peak) ; $i++ ) {
                        $current_peak = $peak[$i];

			if ( abs( $pos[$current_peak] - $pos[$last_peak] ) < $peak_sep )

			  #			&& 	($end[$i]-$start[$i])>2)
			{
				pop @final_peaks;  
				pop @final_starts; 
				pop @final_ends; 

				$current_peak = int( ( $current_peak + $last_peak ) / 2 );
				push @final_peaks,  $current_peak;
				push @final_starts, $last_start;
				push @final_ends,   $end[$i];
				$last_peak  = $current_peak;
				$last_start = $last_start;
				$last_end   = $end[$i];
			}

			#		elsif ( abs($current_peak - $last_peak) > 2 )
			#					&& ($end[$i]-$start[$i]>2))
			else {

				push @final_peaks,  $current_peak;
				push @final_starts, $start[$i];
				push @final_ends,   $end[$i];
				$last_peak  = $current_peak;
				$last_start = $start[$i];
				$last_end   = $end[$i];
                        }
		}
	}

        @sorted_peaks = sort { $ratio[$b] <=> $ratio[$a] } @final_peaks;
	$tallest_idx = $sorted_peaks[0];


        for ($p=0; $p < scalar(@final_peaks); $p++) {

#                        peak start will be nearest left-valley, or start idx of region,
#                       # or position that is at distance $peak_width/2, whichever is less.

                        if ( $final_starts[$p] < $final_peaks[$p] - $peak_width / (2 * $probe_length) ) {
                                $final_starts[$p] =
                                  sprintf( "%.0f", $final_peaks[$p] - $peak_width /(2* $probe_length) -1 );
                        }

             if ($final_peaks[$p] == $tallest_idx) {$peakFlag[$p]="P" }
        
             # if peak's start and end positions have values no more than 50%
             # of peak value, i.e., peak is a well resolved shoulder, then it
             # will be flagged as a Peak instead of a Shoulder.

             elsif ($ratio[$final_starts[$p]] < 0.5* $ratio[$final_peaks[$p]]
               &&
               $ratio[$final_ends[$p]] < 0.5* $ratio[$final_peaks[$p]]
               )
             {$peakFlag[$p] = "P";}
	     else 
             {$peakFlag[$p] = "S";}
        } 

	return \@final_starts, \@final_peaks, \@final_ends, \@peakFlag;

}

sub print_genes {

 # re-reads datafile and sums the ratio and probability values of probes
 # within above-threshold runs. returns array of arrays containing these results
 # for each hit.

	my (@results) = @_;

	print OUT2 "# Parameters used:
# --------------------------------------------------------------------------------------------\n";
	print OUT2
"# Thresholds: log10(signal) = $ratio_threshold, -log2(prob)  = $prob_threshold,
# minrun = $minrun, maxgap1 = $maxgap1, maxgap2 = $maxgap2
# --------------------------------------------------------------------------------------------\n";

	$prom_cnt  = 0;
	$orf_cnt   = 0;
	$genes_cnt = 0;
	@genes     = ();

	for $row (@results) {
		if ( ${$row}[11] ne "" || ${$row}[14] ne "" ) { ++$prom_cnt }
		$orf_cnt = $orf_cnt + ${$row}[28];

		# need to eliminate duplicates from the genes list:
		%seen        = ();
		@uniq        = ();
		@these_genes = split /,/, ${$row}[27];
		push( @genes, @these_genes );

		foreach $item (@genes) {
			if ( $item eq "" ) { next }
			push( @uniq, $item ) unless $seen{$item}++;

		}
	}

	%seen = ();
	@uniq = ();
	foreach $item (@genes) {
		unless ( $seen{$item} ) {
			$seen{$item}++;
			push( @uniq, $item );
			$genes_cnt++;
			print OUT2 $item, "\n";
		}
	}

}
