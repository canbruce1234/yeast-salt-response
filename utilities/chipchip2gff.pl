#!/usr/bin/perl -w
$infile1 = $ARGV[0];
$outfile = $ARGV[1];
$name = $outfile;
$name =~ /^(.*?)\./;
$array = $1;

open( IN1, $infile1 )    || die "can't open $infile1: $!";
open( OUT, ">$outfile" ) || die "can't open $outfile: $!";

@temp = ();
while ( defined( my $line = <IN1> ) ) {
chomp $line;
@items = split /\t/, $line;
$chr = $items[0];
$pos = $items[1];
$signal = $items[3];
$chr =~ s/mito/chr17/g;
$chr =~ s/ chr/chr/g;
$end = $pos+49;
push @temp,[$chr, $array, $pos, $end, $signal];
}
@temp=
	sort{ substr (@$a[0],3) <=> substr(@$b[0],3) 
             || @$a[2] <=> @$b[2]} @temp;
foreach $row (@temp) {
  $chr = @$row[0];
  $array = @$row[1];
  $pos = @$row[2];
  $end = @$row[3];
  $signal = @$row[4];
  print OUT "$chr\t$array\tnimblegen\t$pos\t$end\t$signal\t.\t.\t.\n";
}
 
