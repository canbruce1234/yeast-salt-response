#!/usr/bin/perl -w
#################################################################
$infile         = $ARGV[0];   # chipchip expt output file
$outfile         = $ARGV[1];   # hitlist

open( IN, $infile) || die "can't open $infile: $!";
open( OUT, ">$outfile") || die "can't open $outfile: $!";

$newline = "";
while ( defined( my $line = <IN> ) ) {
	if ( $line =~ /^>/){
		if ($newline ne "") {
			print OUT $newline . "\n";
			$newline = "";
		}
		$line = substr($line, 0, length($line)-1);
		print OUT $line, "\n";
		next;
	}
	
	$line = substr($line, 0, length($line)-2);
	$newline=$newline . $line;
}
close IN;
close OUT;

