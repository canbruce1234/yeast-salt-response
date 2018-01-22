#!/usr/bin/perl -w
###############################
# Makes a gff file using the p-values
#
################################
$infile1 = $ARGV[0];
$outfile = $ARGV[1];

open( IN1, $infile1 )    || die "can't open $infile1: $!";
open( OUT, ">$outfile" ) || die "can't open $outfile: $!";

@temp = ();
$name = $infile1;
$name =~ /^(.*?)\./;
$array = "p-value";
while ( defined( my $line = <IN1> ) ) {
chomp $line;
@items = split /\t/, $line;
$chr = $items[0];
$pos = $items[1];
$pos =~ s/ //g;
$prob = $items[2];
$prob =~ s/ //g;
$chr =~ s/mito/chr17/g;
$chr =~ s/\s?//g;
$end = $pos+49;
push @temp,[$chr, $array, $pos, $end, $prob];
}
@temp=
	sort{ substr (@$a[0],3) <=> substr(@$b[0],3) 
             || @$a[2] <=> @$b[2]} @temp;
foreach $row (@temp) {
  $chr = @$row[0];
  $array = @$row[1];
  $pos = @$row[2];
  $end = @$row[3];
  $prob = @$row[4];
  print OUT "$chr\tnimblegen\t$array\t$pos\t$end\t$prob\t.\t.\t.\n";
} 

 
 
 
