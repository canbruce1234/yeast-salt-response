#!/usr/bin/perl -w
###############################################################
#	Takes two hitlists from a chipchip expt.  Writes a
#	a file with the number of overlaps for the first N hits in
#	each list having the highest signals.
#
###############################################################
do '/home/cbruce/eclipse/workspace/ChipChip/intersect.pl';
do '/home/cbruce/eclipse/workspace/ChipChip/locate.pl';

$infile1 = $ARGV[0];    # hitlist1
$infile2 = $ARGV[1];    # hitlist2
$outfile = $ARGV[2];    # data file to graph ranks with overlaps

open( IN1, $infile1 )    || die "can't open $infile1: $!";
open( IN2, $infile2 )    || die "can't open $infile2: $!";
open( OUT, ">$outfile" ) || die "can't open $outfile: $!";

while ( defined( $line1 = <IN1> ) ) {
	chomp $line1;
	if ( $line1 =~ /^#/ )  { next }
	if ( $line1 =~ /^\n/ ) { next }
	@items = split /	/, $line1;
	push @list1, $items[8];
	push @list1, $items[11];
}

close IN1;
print "\n";
while ( defined( $line2 = <IN2> ) ) {	
	chomp $line2;
	if ( $line2 =~ /^#/ ) { next }
	if ( $line2 =~ /^ / ) { next }
	@items = split /	/, $line2;
	push @list2, $items[8];
	push @list2, $items[11];
	
}
close IN2;

# foreach $row (@list1) {print "@$row[0]\t@$row[1]\t@$row[2]\t@$row[3]\n";}

# Compare two lists of equal length N.
# N ($maxlength) will be the smallest of the lengths of the two lists.
if ( scalar(@list1) > scalar(@list2) ) { $maxlength = scalar(@list2) }
else { $maxlength = scalar(@list1) }
print "maxlength = $maxlength\n";

 #$maxlength=41;
$window = 10;
$rank1        = 0;
#$last = 0;
for ( $rank2 = 9 ; $rank2 < $maxlength  ; $rank2 = $rank2 + 10 ) {
	for ( $idx = $rank1 ; $idx < $rank2 + 1 ; $idx++ ) {
		$sublist1[ $idx - $rank1 ] = $list1[$idx];
	}

	for ( $idx = $rank1 ; $idx < $rank2 + 1 ; $idx++ ) {
		$sublist2[ $idx - $rank1 ] = $list2[$idx];
	}

	@union = @isect = ();
	%union = %isect = ();

	foreach $e (@sublist1) { $union{$e} = 1 ; 	
		#print"$e|";
		}
	#print"\n";

	foreach $e (@sublist2) {
		#print"$e!";
		if ( $union{$e} ) { $isect{$e} = 1 }
		$union { $e } = 1;
	}

#print "\n";
	@union = keys %union;
	@isect = keys %isect;
	
	$isect_cnt{$rank2} = scalar(@isect);
	$rank = $rank2 +1;
	if ($rank2-$window > 0) {$delta = $isect_cnt{$rank2} - $isect_cnt{$rank2-$window} }
	else {$delta = 0 }
	print OUT "$rank $isect_cnt{$rank2} $delta\n"; 
	
}



close OUT;

