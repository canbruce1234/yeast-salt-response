#!/usr/bin/perl -w

$infile1 = $ARGV[0];
$outfile1 = $ARGV[1];   # rename.out

open( IN1, $infile1 ) || die "can't open $infile1: $!";
open( OUT1, ">$outfile1" ) || die "can't open $outfile1: $!";

my %map = (
        "chrI" => "chr1",
	"chrII" => "chr2",
	"chrIII" =>  "chr3",
	"chrIV" => "chr4",
	"chrV" => "chr5",
	"chrVI" => "chr6",
	"chrVII" => "chr7",
	"chrVIII" => "chr8",
	"chrIX" => "chr9",
	"chrX" => "chr10",
	"chrXI" => "chr11",
	"chrXII" => "chr12",
	"chrXIII" =>"chr13",
	"chrXIV" => "chr14",
	"chrXV" => "chr15",
	"chrXVI" => "chr16",
	"mito" => "chr17",
);

	while ( defined( my $line = <IN1> ) ) {
		chomp $line;
		@items = split /\t/, $line;
		$chr = $items[0];
		$chr =~ s/\s?//g;
		$chr = $map{$chr};
		$pos   = $items[1];
		$pos =~ s/\s?//g;
		$red  = $items[2];
		$red =~ s/\s?//g;
		$green = $items[3];
		$green =~ s/\s?//g;

	print OUT1 "$chr\t$pos\t$red\t$green\n";	
	}
	close IN1;
	close OUT1;
