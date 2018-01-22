#!/usr/bin/perl -w
#######################################################
# reads a chipchip hitlist, reads hit lengths in order,
# then generates random hits in the yeast genome having
# the same hit lengths. Produces a hitlist with hit locations.
#
#
########################################################

do '/home/cbruce/eclipse/workspace/ChipChip/locate.pl';

$infile1  = $ARGV[0];    # hitlist1
$outfile  = $ARGV[1];    # begin and end points of hits; nearby promoters
$outfile2 = $ARGV[2];    # genes and scores within hits

open( IN1,  $infile1 )     || die "can't open $infile1: $!";
open( OUT,  ">$outfile" )  || die "can't open $outfile: $!";
open( OUT2, ">$outfile2" ) || die "can't open $outfile2: $!";

while ( defined( $line = <IN1> ) ) {
	if ( $line =~ /^#/ ) { next }
	if ( $line =~ /^ / ) { next }
	$line =~ s/ chr/chr/g;      #
	$line =~ s/mito/chr17/g;    #
	@fields = split /\t/, $line;
	$length = $fields[4];
	push @lengths, $length;

}

%chr_length = (
	"chr1"  => "230208",
	"chr2"  => "813178",
	"chr3"  => "316616",
	"chr4"  => "1531916",
	"chr5"  => "576869",
	"chr6"  => "270148",
	"chr7"  => "1090946",
	"chr8"  => "562642",
	"chr9"  => "439885",
	"chr10" => "745666",
	"chr11" => "666454",
	"chr12" => "1078174",
	"chr13" => "924429",
	"chr14" => "784331",
	"chr15" => "1091287",
	"chr16" => "948062",
	"chr17" => "85779",
);

$cumm_sum{"chr1"}  = $chr_length{"chr1"};
$cumm_sum{"chr2"}  = $chr_length{"chr2"} + $cumm_sum{"chr1"};
$cumm_sum{"chr3"}  = $chr_length{"chr3"} + $cumm_sum{"chr2"};
$cumm_sum{"chr4"}  = $chr_length{"chr4"} + $cumm_sum{"chr3"};
$cumm_sum{"chr5"}  = $chr_length{"chr5"} + $cumm_sum{"chr4"};
$cumm_sum{"chr6"}  = $chr_length{"chr6"} + $cumm_sum{"chr5"};
$cumm_sum{"chr7"}  = $chr_length{"chr7"} + $cumm_sum{"chr6"};
$cumm_sum{"chr8"}  = $chr_length{"chr8"} + $cumm_sum{"chr7"};
$cumm_sum{"chr9"}  = $chr_length{"chr9"} + $cumm_sum{"chr8"};
$cumm_sum{"chr10"} = $chr_length{"chr10"} + $cumm_sum{"chr9"};
$cumm_sum{"chr11"} = $chr_length{"chr11"} + $cumm_sum{"chr10"};
$cumm_sum{"chr12"} = $chr_length{"chr12"} + $cumm_sum{"chr11"};
$cumm_sum{"chr13"} = $chr_length{"chr13"} + $cumm_sum{"chr12"};
$cumm_sum{"chr14"} = $chr_length{"chr14"} + $cumm_sum{"chr13"};
$cumm_sum{"chr15"} = $chr_length{"chr15"} + $cumm_sum{"chr14"};
$cumm_sum{"chr16"} = $chr_length{"chr16"} + $cumm_sum{"chr15"};
$cumm_sum{"chr17"} = $chr_length{"chr17"} + $cumm_sum{"chr16"};

$total_length = $cumm_sum{"chr17"};

@random_hits = ();
print OUT
"# Chr	Start	Mid	End	mLength	#probes mProb	mSignal	geneW	akaW	distW	geneC akaC	distC	orf\n";
LABEL1:
foreach $length (@lengths) {
	$start = sprintf( "%.0f", rand($total_length) );
	( $start_chr, $start_pos ) = locate_random_hit($start);
	( $end_chr,   $end_pos )   = locate_random_hit( $start + $length );
	if ( $start_chr ne $end_chr ) { redo LABEL1 }
	foreach $hit (@random_hits) {
		if ( @$hit[0] ne $start_chr ) { next }
		if ( @$hit[1] >= $start_pos && @$hit[1] <= $end_pos ) { redo LABEL1 }
		if ( @$hit[2] >= $start_pos && @$hit[2] <= $end_pos ) { redo LABEL1 }
		push @random_hits, [ $start_chr, $start_pos, $end_pos ];
	}

	( $geneW, $akaW, $distW, undef, $geneC, $akaC, $distC, undef, $genes, $orf )
	  = &locate_run_yeast( $chr, $start_pos, $end_pos );

	@temp = ( $start_chr, $start_pos, $end_pos );

	push @temp_results, [@temp];
}

@results = locate_run_yeast2(@temp_results);

$genes_cnt = 0;

foreach $row (@results) {
	(
		$chr,   $start_pos, $end_pos, $geneW, $akaW, $distW,
		$geneC, $akaC,      $distC,   $genes, $orf
	  )
	  = (
		@$row[0], @$row[1], @$row[2], @$row[3], @$row[4], @$row[5],
		@$row[6], @$row[7], @$row[8], @$row[9], @$row[10]
	  );

	# print "*** $genes", "\n";
	$mid_pos = int( ( $start_pos + $end_pos ) / 2 );

	$length = $end_pos - $start_pos;
	$chr    = "chr" . $chr;

	if ( $akaW ne "" ) { $akaW = "(" . $akaW . ")" }
	if ( $akaC ne "" ) { $akaC = "(" . $akaC . ")" }

	@these_genes = split /,/, $genes;
	push( @genes, @these_genes );

	print OUT
"$chr\t$start_pos\t$mid_pos\t$end_pos \t$length\t0\t0\t0\t$geneW\t$akaW\t$distW\t$geneC\t$akaC\t$distC\t$orf\n";
}    

%seen = ();
@uniq = ();
foreach $item (@genes) {
	unless ( $seen{$item} ) {
		$seen{$item}++;
		push( @uniq, $item );
		$genes_cnt ++;
		print OUT2 $item, "\n";
	}
}
print OUT2 "\n";


close IN1;
close OUT;

close OUT2;

sub locate_random_hit {
	my ($rnd) = (@_);

	if ( $rnd <= $cumm_sum{"chr1"} ) { $chr = "chr1"; $pos = $rnd; }
	elsif ( $rnd > $cumm_sum{"chr1"} && $rnd <= $cumm_sum{"chr2"} ) {
		$chr = "chr2";
		$pos = $rnd - $cumm_sum{"chr1"};
	} elsif ( $rnd > $cumm_sum{"chr2"} && $rnd <= $cumm_sum{"chr3"} ) {
		$chr = "chr3";
		$pos = $rnd - $cumm_sum{"chr2"};
	} elsif ( $rnd > $cumm_sum{"chr3"} && $rnd <= $cumm_sum{"chr4"} ) {
		$chr = "chr4";
		$pos = $rnd - $cumm_sum{"chr3"};
	} elsif ( $rnd > $cumm_sum{"chr4"} && $rnd <= $cumm_sum{"chr5"} ) {
		$chr = "chr5";
		$pos = $rnd - $cumm_sum{"chr4"};
	} elsif ( $rnd > $cumm_sum{"chr5"} && $rnd <= $cumm_sum{"chr6"} ) {
		$chr = "chr6";
		$pos = $rnd - $cumm_sum{"chr5"};
	} elsif ( $rnd > $cumm_sum{"chr6"} && $rnd <= $cumm_sum{"chr7"} ) {
		$chr = "chr7";
		$pos = $rnd - $cumm_sum{"chr6"};
	} elsif ( $rnd > $cumm_sum{"chr7"} && $rnd <= $cumm_sum{"chr8"} ) {
		$chr = "chr8";
		$pos = $rnd - $cumm_sum{"chr7"};
	} elsif ( $rnd > $cumm_sum{"chr8"} && $rnd <= $cumm_sum{"chr9"} ) {
		$chr = "chr9";
		$pos = $rnd - $cumm_sum{"chr8"};
	} elsif ( $rnd > $cumm_sum{"chr9"} && $rnd <= $cumm_sum{"chr10"} ) {
		$chr = "chr10";
		$pos = $rnd - $cumm_sum{"chr9"};
	} elsif ( $rnd > $cumm_sum{"chr10"} && $rnd <= $cumm_sum{"chr11"} ) {
		$chr = "chr11";
		$pos = $rnd - $cumm_sum{"chr10"};
	} elsif ( $rnd > $cumm_sum{"chr11"} && $rnd <= $cumm_sum{"chr12"} ) {
		$chr = "chr12";
		$pos = $rnd - $cumm_sum{"chr11"};
	} elsif ( $rnd > $cumm_sum{"chr12"} && $rnd <= $cumm_sum{"chr13"} ) {
		$chr = "chr13";
		$pos = $rnd - $cumm_sum{"chr12"};
	} elsif ( $rnd > $cumm_sum{"chr13"} && $rnd <= $cumm_sum{"chr14"} ) {
		$chr = "chr14";
		$pos = $rnd - $cumm_sum{"chr13"};
	} elsif ( $rnd > $cumm_sum{"chr14"} && $rnd <= $cumm_sum{"chr15"} ) {
		$chr = "chr15";
		$pos = $rnd - $cumm_sum{"chr14"};
	} elsif ( $rnd > $cumm_sum{"chr15"} && $rnd <= $cumm_sum{"chr16"} ) {
		$chr = "chr16";
		$pos = $rnd - $cumm_sum{"chr15"};
	} elsif ( $rnd > $cumm_sum{"chr16"} && $rnd <= $cumm_sum{"chr17"} ) {
		$chr = "chr17";
		$pos = $rnd - $cumm_sum{"chr16"};
	}

	return ( $chr, $pos );
}

