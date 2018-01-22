#!/usr/bin/perl -w
#################################################################
# Reads file with chr, pos, red, green
# Produces output file of sgr format:
# chr, pos, log(red/green). 
# 
#########################################################

$infile1 = $ARGV[0];
$outfile1 = $ARGV[1];   # ratio.out

open( IN1, $infile1 ) || die "can't open $infile1: $!";
open( OUT1, ">$outfile1" ) || die "can't open $outfile1: $!";

	while ( defined( my $line = <IN1> ) ) {
		$line =~ /^(.*?)\t(.*?)\t(.*?)\t(.*?)\n$/;
		$chr = $1;
		$pos   = $2;
		$red  = $3;
		$green = $4;
		$ratio= log($red/$green)/log(2);
	print OUT1 "$chr\t$pos\t$ratio\n";	
	}
	close IN1;
	close OUT1;

