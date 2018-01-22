#!/usr/bin/perl -w
##############################################################
#
# Reads in two data files (for Red and Green channels) from a
# chipchip experiment.  Gathers the chromosome numbers, sequence
# position of probe, red and green signals and outputs them in
# an output file with 6 cols, to be used by R script for analysis:
# columns are chr, pos, x, y, Red, Green
# 
##############################################################

$infile1 = $ARGV[0];
$infile2 = $ARGV[1];
$outfile = $ARGV[2];

open( IN1, $infile1 )    || die "can't open $infile1: $!";
open( IN2, $infile2 )    || die "can't open $infile2: $!";
open( OUT, ">$outfile" ) || die "can't open $outfile: $!";

@Red   = ();
@Green = ();

my %map = (
	"NC_001133" => "chr1",
	"NC_001134" => "chr2",
	"NC_001135" => "chr3",
	"NC_001136" => "chr4",
	"NC_001137" => "chr5",
	"NC_001138" => "chr6",
	"NC_001139" => "chr7",
	"NC_001140" => "chr8",
	"NC_001141" => "chr9",
	"NC_001142" => "chr10",
	"NC_001143" => "chr11",
	"NC_001144" => "chr12",
	"NC_001145" => "chr13",
	"NC_001146" => "chr14",
	"NC_001147" => "chr15",
	"NC_001148" => "chr16",
	"NC_001224" => "chr17",
);

while ( defined( $line = <IN1> ) ) {
	if ( $line =~ /^#/ )        { next }
	if ( $line =~ /^IMAGE_ID/ ) { next }    # skip column headers

	$line =~
/^(.*?)\t(.*?)\t(.*?)\s(.*?)\s(.*?)\s(.*?)\s(.*?)\s(.*?)\s(.*?)\s(.*?)\s(.*?)/;

	next if ( $2 eq "RANDOM" ) ;      # skip random data
	$probe  = $4;
	$pos    = $5;
        $x = $6;
        $y = $7;
	$signal = $9;
	$probe =~ /(NC_001...).*/;
	$ncbi = $1;
	$chr  = $map{$ncbi};
	push @Red, [ $chr, $pos, ,$x, $y, $signal ];
}
close IN1;

while ( defined( $line = <IN2> ) ) {
	if ( $line =~ /^#/ )        { next }
	if ( $line =~ /^IMAGE_ID/ ) { next }    # skip column headers
	$line =~
/^(.*?)\t(.*?)\t(.*?)\s(.*?)\s(.*?)\s(.*?)\s(.*?)\s(.*?)\s(.*?)\s(.*?)\s(.*?)/;

	next if ( $2 eq "RANDOM" );    # skip random data
	$probe  = $4;
	$pos    = $5;
        $x = $6;
        $y = $7;
	$signal = $9;
	$probe =~ /(NC_001...).*/;
	$ncbi = $1;
	$chr  = $map{$ncbi};
	push @Green, [ $chr, $pos, $signal ];
}
close IN2;

if ( scalar(@Green) != scalar(@Red) ) {
	print "Arrays are of unequal length. \n";
	exit;
}
for ( $i = 0 ; $i < scalar(@Red) ; $i++ ) {

	 if($Red[$i][0] eq $Green[$i][0] && $Red[$i][1] == $Green[$i][1]) {
	print OUT "$Red[$i][0]\t$Red[$i][1]\t$Red[$i][2]\t$Red[$i][3]\t$Red[$i][4]\t$Green[$i][2]\t\n";

	 }
	else {
		print
		  "data mismatch: Red:  $Red[$i][0]\t$Red[$i][1]\t$Red[$i][2]\n
								Green:$Green[$i][0]\t$Green[$i][1]\t$Green[$i][2]\n";
		exit;
	}
}

close OUT;

