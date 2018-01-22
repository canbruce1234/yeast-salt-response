#!/usr/bin/perl -w
#################################################################
# Reads four column output from R script that analyzes chipchip data;
# Splits data into two files of sgr format, having columns
# chr, pos, (ratio/prob). 
# 
# 8/23/05 (CB): Created.
# 9/5/05 CB: simplified analyze.pl 
#########################################################

$infile1 = $ARGV[0];
$outfile1 = $ARGV[1];   # rename.out

open( IN1, $infile1 ) || die "can't open $infile1: $!";
open( OUT1, ">$outfile1" ) || die "can't open $outfile1: $!";

my %map = (
        "chr1" => "chrI",
	"chr2" => "chrII",
	"chr3" => "chrIII",
	"chr4" => "chrIV",
	"chr5" => "chrV",
	"chr6" => "chrVI",
	"chr7" => "chrVII",
	"chr8" => "chrVIII",
	"chr9" => "chrIX",
	"chr10" => "chrX",
	"chr11" => "chrXI",
	"chr12" => "chrXII",
	"chr13" => "chrXIII",
	"chr14" => "chrXIV",
	"chr15" => "chrXV",
	"chr16" => "chrXVI",
	"mito" => "chrM",
);
	while ( defined( my $line = <IN1> ) ) {
		$line =~ /^(.*?) \t (.*?) \t (.*?) \t (.*?) \n$/;
		$chr = $1;
		$pos   = $2;
		$prob  = $3;
		$ratio = $4;
        # print ",$chr,$pos,$prob,$ratio,\n";

	print OUT1 "$map{$chr}\t$pos\t$prob\t$ratio\n";	
	}
	close IN1;
	close OUT1;

