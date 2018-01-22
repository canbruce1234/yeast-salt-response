#!/usr/bin/perl -w
###############################################################
# Compare a sorted gene listed from a microarray experiment
# with a coordinate list from a chipchip experiment.  Converts
# the gene names to a coordinate list ( end of upstream gene to
# end of given gene ), then looks for overlaps with the coordinates
# in the chipchip hit list.  For each region pair, looks for
# coordinates of common sequence between the two regions.  If
# this exists, returns:
#	1) sum of ranks for the overlapping region.
#       2) ORF name, common name, strand, start and end from expression expt.
#	3) chipchip hit start, end, length and midpoint;
#	4) location of chipchip hit midpoint relative to ORF start from expn expt
#
#	THIS VERSION INTERSECTS A TF-TF INTERSECT FILE WITH A RANKED GENE NAME LIST
#
###############################################################
use DBI;

my $dbh = DBI->connect( "dbi:mysql:chipchip", "cbruce", "yalejohn" )
  or die "Can't connect to mysql database: $DBI::errstr\n";

$infile1 = $ARGV[0];    # chipchip expt output file
$infile2 = $ARGV[1];    # expression expt output file
$outfile = $ARGV[2];    # common genes files

open( IN1, $infile1 )    || die "can't open $infile1: $!";
open( IN2, $infile2 )    || die "can't open $infile2: $!";
open( OUT, ">$outfile" ) || die "can't open $outfile: $!";

$chipchip_rank = 0;

while ( defined( $line = <IN1> ) ) {
	if ( $line =~ /^#/ )      { next }
	if ( $line =~ /^\s*?\n/ ) { next }
	$line =~ s/ chr/chr/g;      #
	$line =~ s/mito/chr17/g;    #
	@fields = split /\t/, $line;
	$list1_chr = $fields[0];
	$list1_chr =~ s/ //g;
	$list1_begin = $fields[1];
	$list1_begin =~ s/ //g;
	$list1_peak = $fields[2];
	$list1_peak =~ s/ //g;
	$list1_end = $fields[3];
	$list1_end =~ s/ //g;
	$list1_chr = substr( $list1_chr, 3 );

	$chipchip_rank++;
	@temp = ( $list1_chr, $list1_begin, $list1_end, $chipchip_rank, $list1_peak );
	push @list1, [@temp];
}
close IN1;

$rank  = 0;
@list2 = ();
while ( defined( $line = <IN2> ) ) {
	chomp $line;
	if ( $line =~ /^#/ )      { next }
	if ( $line =~ /^\s*?\n/ ) { next }
	@fields = split /	/, $line;

	$name = $fields[0];
	$signal = sprintf( "%.2f", $fields[1] );
	$signal =~ s/ //g;
	$name =~ s/"//g;
	$name =~ s/ //g;

	$stmt1 = $dbh->prepare(
		qq{	 SELECT vvv, start, end, upstr_end, strand, name, aka, descr
               FROM yeast g
	 		  WHERE name = '$name' or aka = '$name' or zzz like '%$name%'
	 		  }
	  )
	  or die "can't prepare stmt: $DBI::errstr";

	$stmt1->execute()
	  or die "Can't execute SQL statement: $DBI::errstr\n";

	my @row1 = $stmt1->fetchrow_array();
	warn "Data fetching terminated early by error: $DBI::errstr\n" if $DBI::err;
if (scalar(@row1)==0){print "no output for $name !!!! \n"}
	( $chr, $start, $end, $upstr_end, $strand, $name, $aka, $descr ) = @row1;

	if ( !defined $upstr_end ) {

		#print "$name has no upstr_end\n";
		$upstr_end = $end;
	}

	$rank++;
	@temp = (
		$chr,  $start, $end,   $upstr_end, $strand,
		$name, $aka,   $descr, $signal,    $rank
	);

# print join(@temp),"\n";

	push @list2, [@temp];

}
close IN2;
$stmt1->finish;

FOR_J: for $j ( 0 .. $#list1 ) {

	$a_chr   = $list1[$j][0];
	$a_begin = $list1[$j][1];
	$a_end   = $list1[$j][2];
	$a_rank  = $list1[$j][3];
	$a_peak  = $list1[$j][4];

  FOR_K: for $k ( 0 .. $#list2 ) {
		$b_chr       = $list2[$k][0];
		$b_begin     = $list2[$k][1];
		$b_end       = $list2[$k][2];
		$b_upstr_end = $list2[$k][3];
		$b_strand    = $list2[$k][4];
		$b_name      = $list2[$k][5];
		$b_aka       = $list2[$k][6];
		$b_descr     = $list2[$k][7];
		$b_signal    = $list2[$k][8];
		$b_rank      = $list2[$k][9];

		# need to switch upstr_end with end for Crick strand genes
		if ( $b_strand eq "C" ) {
			$temp        = $b_upstr_end;
			$b_upstr_end = $b_end;
			$b_end       = $temp;
		}

		if ( $b_chr ne $a_chr ) { next }

		if ( $b_upstr_end > $a_end ) {
			next FOR_K;
		}    # skip loop if no overlap

		if ( $b_end <= $a_begin ) {
			next FOR_K;
		}    # skip loop if no overlap

		$sum_rank = $a_rank + $b_rank;
		if ( $b_strand eq "C" ) {
			$temp        = $b_upstr_end;
			$b_upstr_end = $b_end;
			$b_end       = $temp;
		}

		$hit_length = $a_end - $a_begin + 1;
		#$hit_midpoint = sprintf( "%.0f", ( $a_begin + $a_end ) / 2 );
		if ( $b_strand eq "W" ) {
		#	$rel_pos = $hit_midpoint - $b_begin;
			$rel_pos = $a_peak - $b_begin;
		} else {
		#	$rel_pos = $b_begin - $hit_midpoint;
			$rel_pos = $b_begin - $a_peak;
		}

		push @overlap,
		  [
			(
				$sum_rank,    $b_rank,   $a_rank,  $b_name,     $b_aka,
				$b_descr,     $b_signal, $rel_pos, $hit_length, $b_chr,
				$b_upstr_end, $b_begin,  $b_end,   $a_begin,    $a_end
			)
		  ];

	}    # end FOR_K

}    # end FOR_J

#sort @uniq_overlap by sum_rank:
@overlap = sort { @$a[0] <=> @$b[0] } @overlap;

$size1       = scalar(@list1);
$size2       = scalar(@list2);
$common_size = scalar(@overlap);

print OUT "# SUMMARY\n";
print OUT "# file1 ($size1 hits): $infile1 \n";
print OUT "# file2 ($size2 genes): $infile2 \n";
print OUT "# output ($common_size in common): $outfile\n";

print "# file1 ($size1 hits):\t$infile1 \n";
print "# file2 ($size2 genes):\t$infile2 \n";
print "# output ($common_size in common):\t$outfile\n";

print OUT "#", "-" x 110, "\n";
print OUT
  "# sum\texpr\thit\t\t\texpr\thit\thit\t      |\tupstr\t\t      |\thit\thit\n";
print OUT
"# rank\trank\trank\tgene\taka\tsignal\tpeak\tlength\tchr   |\tend\tstart\tend   |\tstart\tend\n";
print OUT "#", "-" x 110, "\n";

print  "#", "-" x 110, "\n";
print  "# sum\texpr\thit\t\t\texpr\thit\thit\t      |\tupstr\t\t      |\thit\thit\n";
print  
"# rank\trank\trank\tgene\taka\tsignal\tpeak\tlength\tchr   |\tend\tstart\tend   |\tstart\tend\n";
print "#", "-" x 110, "\n";

for $row (@overlap) {
	(
		$sum_rank,    $b_rank,   $a_rank,  $b_name,     $b_aka,
		undef,        $b_signal, $rel_pos, $hit_length, $b_chr,
		$b_upstr_end, $b_begin,  $b_end,   $a_begin,    $a_end
	  )
	  = @$row;
	print "$sum_rank\t$b_rank\t$a_rank\t$b_name\t$b_aka\t$b_signal\t",
          "$rel_pos\t$hit_length\t$b_chr\t$b_upstr_end\t$b_begin\t$b_end\t$a_begin\t$a_end\n" ;
	print OUT "$sum_rank\t$b_rank\t$a_rank\t$b_name\t$b_aka\t$b_signal\t",
          "$rel_pos\t$hit_length\t$b_chr\t$b_upstr_end\t$b_begin\t$b_end\t$a_begin\t$a_end\n" ;
}
print "#", "-" x 100, "\n";
for $row (@overlap) {
	print OUT "# @$row[3] @$row[4]: @$row[5]\n";
}

print OUT "#", "-" x 100, "\n";
for $row (@overlap) {
	print "@$row[3] @$row[4]: @$row[5]\n";
}



close OUT;

$dbh->disconnect or warn "Disconnection failed: $DBI::errstr\n";

