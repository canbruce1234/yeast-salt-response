#!/usr/bin/perl -w
$infile1 = $ARGV[0];
$infile2 = $ARGV[1];
$outfile = $ARGV[2];

open( IN1, $infile1 ) || die "can't open $infile1: $!";
open( IN2, $infile2 ) || die "can't open $infile2: $!";
open( OUT, ">$outfile" ) || die "can't open $outfile: $!";

@masks = ();
while ( defined( $line = <IN1> ) ) {
	if ( $line =~ /^#/ ) { next }
	chomp $line;
	$line =~ /^(.*?)\t(.*?)\t(.*?)\t(.*?)$/;
	$chr = $1;
	if ( $2 < $3 ) {
		$start = $2;
		$end   = $3;
	}
	else {
		$start = $3;
		$end   = $2;
	}
	push @masks, [ $chr, $start, $end ];
}
close IN1;

@sorted_masks =
  sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @masks;

#for $i (0..$#sorted_masks){
#	print "@{$sorted_masks[$i]}","\n";
#}


while ( defined( $line = <IN2> ) ) {
	if ( $line =~ /^#/ ) { next }
	chomp $line;
	$line =~ /^(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)$/;
	$chr   = substr( $1, 3 );
	$pos   = $2;
	$x     = $3;
	$y     = $4;
	$red   = $5;
	$green = $6;
	push @data, [ $chr, $pos, $x, $y, $red, $green ];

	#	print "$chr, $pos, $x, $y, $red, $green \n";
}
close IN2;

@sorted_data =
  sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @data;

#for $i (0..$#sorted_data){
#	print "@{$sorted_data[$i]}","\n";
#}

$mask_flag = 0;
$i=0;
$imax=$#sorted_masks;
LINEj: for $j ( 0 .. $#sorted_data ) {
#  		print "i=$i, j=$j new mask:", "@{$sorted_masks[$i]}", "\n";
  		if (   $sorted_data[$j][0] == $sorted_masks[$i][0]
			&& $sorted_data[$j][1] > $sorted_masks[$i][1]
			&& $sorted_data[$j][1] > $sorted_masks[$i][2] ){
				$i++;
			}
		if (   $sorted_data[$j][0] == $sorted_masks[$i][0]
			&& $sorted_data[$j][1] > $sorted_masks[$i][1]
			&& $sorted_data[$j][1] < $sorted_masks[$i][2] )
		{
#			print "i=$i, j=$j ", "@{$sorted_masks[$i]}","*","@{$sorted_data[$j]}",  "\n";
			$mask_flag = 1;
			next LINEj;
		}
		else {
			if ( $mask_flag == 1 ) {
#				print "@{$sorted_data[$j]}",  "\n";
				$mask_flag = 0;
				if ($i<$imax){
					$i++;
#					print "i=$i, new mask:", "@{$sorted_masks[$i]}", "\n";
				}
#				else
#				{last};
				next LINEj;
			}
			else {
				print OUT "chr",join ("\t", @{$sorted_data[$j]}),  "\n";
				next LINEj;
			}
		}
}

