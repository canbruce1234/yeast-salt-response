#!/usr/bin/perl -w
# do '/home/cbruce/eclipse/workspace/ChipChip/intersect.pl';

sub intersect {

my ( $min_overlap, $aref, $bref ) = @_;
@a = @$aref;
@b = @$bref;

# foreach $row (@a) {print "A:@$row[0],@$row[1],@$row[2], @$row[4]\n";}
# foreach $row (@b) {print "B:@$row[0],@$row[1],@$row[2], @$row[4]\n";}
# print"***************\n";

FOR_J: for $j ( 0 .. $#a ) {

	$a_chr       = $a[$j][0];
	$a_begin = $a[$j][1];
	$a_end   = $a[$j][2];
	$a_rank = $a[$j][4];

  FOR_K: for $k ( 0 .. $#b ) {
		$b_chr = $b[$k][0];
		$b_begin = $b[$k][1];
		$b_end   = $b[$k][2];
		$b_rank = $b[$k][4];
		if ( $b_chr ne $a_chr ) {next }

		if ( $b_begin > $a_end ){	
			next FOR_J;	}    # skip loop if no overlap
		if ( $b_end < $a_begin ){	
			next FOR_K;	}    # skip loop if no overlap
		     #
		     #		a |---|
		     #		b   |---|
		     #
		if (   $a_begin < $b_begin
			&& $a_end <= $b_end
			&& $a_end > $b_begin
			&& $a_end - $b_begin > $min_overlap )
		{
			$c_chr        = $b_chr;
			$c_begin      = $a_begin;
			$c_end        = $b_end;
			$overlap	  = $a_end - $b_begin+1;
			$sum_rank     = $a_rank + $b_rank;
		}

		#
		#		a   |---|
		#       b |---|
		#
		if (   $a_begin > $b_begin
			&& $a_end > $b_end
			&& $b_end > $a_begin
			&& $b_end - $a_begin > $min_overlap )
		{
			$c_chr        = $b_chr;
			$c_begin      = $b_begin;
			$c_end        = $a_end;
			$overlap	  = $b_end - $a_begin+1;
			$sum_rank     = $a_rank + $b_rank;
		}

		#
		#	   a   |---|
		#      b |-------|
		#
		if (   $a_begin > $b_begin
			&& $a_end <= $b_end
			&& $a_end - $a_begin > $min_overlap )
		{   
			$c_chr        = $b_chr;
			$c_begin      = $b_begin;
			$c_end        = $b_end;
			$overlap	  = $a_end - $a_begin+1;
			$sum_rank     = $a_rank + $b_rank;
		}

		#
		#	  a |-------|
		#     b   |---|
		#
		if (   $a_begin <= $b_begin
			&& $a_end > $b_end
			&& $b_end - $b_begin > $min_overlap )
		{
			$c_chr        = $b_chr;
			$c_begin      = $a_begin;
			$c_end        = $b_end;
			$overlap	  = $b_end - $b_begin+1;
			$sum_rank     = $a_rank + $b_rank;
		}
		#
		#	  a |-------|
		#     b |-------|
		#
		if (   $a_begin == $b_begin
			&& $a_end == $b_end )
		{
			$c_chr        = $b_chr;
			$c_begin      = $a_begin;
			$c_end        = $a_end;
			$overlap	  = $a_end - $a_begin+1;
			$sum_rank     = $a_rank + $b_rank;
		}

		if (defined $c_begin )
			{
			push @c, [( $c_chr, $c_begin, $c_end, $overlap, $c_end - $c_begin+1 , $sum_rank, $a_rank, $b_rank)];
			 # print "overlap: $c_chr, $c_begin, $c_end, $overlap\n";
			undef $c_chr; undef $c_begin; undef $c_end ;
			}
	}    # end FOR_K

}    # end FOR_J


# need to eliminate duplicates from the list:
%seen = ();
@uniq = ();
foreach $item (@c) {
	@$item[0] =~ s/ //g;
	@$item[1] =~ s/ //g;
	@$item[2] =~ s/ //g;
	$position_key = join '|', @$item[0] , @$item[1] ,  @$item[2];
	# print "$position_key = @$item[0] * @$item[1] * @$item[2]\n";
	push(@uniq, $item) unless $seen{$position_key}++;
	}
return @uniq;
}
