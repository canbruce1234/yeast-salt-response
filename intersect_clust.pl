#!/usr/bin/perl -w
##############################################################
# Takes a chipchip hitlists and a clustered hitlist
# and returns overlapping hits.
#
##############################################################
do '/home/cbruce/eclipse/workspace/ChipChip/intersect2.pl';

 $infile1 = $ARGV[0];	# a whole hitlist1
 $infile2 = $ARGV[1];	# a clustered hitlist (chr, start, end, geneW1, akaW1, statusW1, distW1, etc)
 $outfile = $ARGV[2];	# begin and end points of common hits
 
open( IN1, $infile1 )    || die "can't open $infile1: $!";
open( IN2, $infile2 )    || die "can't open $infile2: $!";
open( OUT, ">$outfile" ) || die "can't open $outfile: $!";

while ( defined( $line = <IN1> ))  {
	$list1_chr=""; 
	$list1_begin=""; 
	$list1_end=""; 
	$list1_geneW1=""; 
	$list1_geneW2=""; 
	$list1_geneC1=""; 
	$list1_geneC2="";
	if ( $line =~ /^#/ )        { next }
	if ( $line =~ /^\n/ )        { next }
	@fields = split /\t/, $line;
	$list1_chr   = $fields[0] ;
	$list1_chr = substr( $list1_chr, 3 );
	$list1_begin = $fields[1];
	$list1_end   = $fields[3];
	$list1_geneW1   	 = $fields[8];
	$list1_geneW2		 = $fields[12];	
	$list1_geneC1      = $fields[16];
	$list1_geneC2      = $fields[20];

	@temp = ($list1_chr, $list1_begin, $list1_end, $list1_geneW1, $list1_geneW2, $list1_geneC1, $list1_geneC2);
	push @list1, [@temp];
 }

close IN1;

while ( defined( $line = <IN2> ))  {
	$list2_chr=""; 
	$list2_begin=""; 
	$list2_end=""; 
	$list2_geneW1=""; 
	$list2_geneW2=""; 
	$list2_geneC1=""; 
	$list2_geneC2="";
	if ( $line =~ /^#/ )        { next }
	if ( $line =~ /^ / )        { next }
	$line =~ s/ //g;
	#print $line, "\n";
	@fields = split /	/, $line;
	$list2_chr   = $fields[0];
	$list2_begin = $fields[1];
	$list2_end   = $fields[2];
	$list2_geneW1   	 = $fields[8];
	$list2_geneW2		 = $fields[11];	
	$list2_geneC1      = $fields[15];
	$list2_geneC2      = $fields[19];
	
	@temp = ($list2_chr, $list2_begin, $list2_end, $list2_geneW1, $list2_geneW2, $list2_geneC1, $list2_geneC2);
	push @list2, [@temp];
}
close IN2;


@overlap = &intersect2(\@list1, \@list2, );	
print OUT "#\tTF\t\tYAP4_cluster\n";
print OUT "#Chr	Start	End	Start	End  	gene_W1 gene_W2 gene_C1 gene_C2\n";
 for $row (@overlap) {
 	 #print split "\t", @$row;
	 print OUT  "@$row[0]\t@$row[1]\t@$row[2]\t@$row[3]\t@$row[4]\t@$row[5]\t@$row[6]\t@$row[7]\t@$row[8]\n";
 }
close OUT;
