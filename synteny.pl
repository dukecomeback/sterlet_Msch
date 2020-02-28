#!/usr/bin/perl -w

my $usage=<<EOF;
--------------------------------
I keep ortholog pairs that are in synteny. If you are using me in -ohno model, I will make sure that the pair are not in the same chromosome to avoid tandem duplication cluster.

Here is how I do it:
	sort -V the second last column and then walk through lines and check the last column:
	1) initiating: when there are \$egg (default 5) genes in a row not melted, initiate the chain
	2) fusing mechanism: when the chain is interupt by a gap of \$gap (default 15) genes, melt the chain  

Usage: perl $0 all.pep.bla.best.RBH.pos.sort (-egg 5) (-gap 15) (-ohno) >all.pep.bla.best.RBH.pos.sort.ohno
	-egg	how many un-melted genes is needed for iniciating the chain  
	-gap 	how many gene is needed to melt the chain
	-ohno	ohnolog model to identify ohnolog pairs, which should be located in different chromosomes

Demo of the file "all.pep.bla.best.RBH.pos.sort" (spe_chr_loci):
rem_g1004.t1	ste_g42369.t1	rem_246_1	ste_42_50
rem_g1008.t1	ste_g42365.t1	rem_246_4	ste_42_48
rem_g1011.t1	ste_g56254.t1	rem_247_1	ste_U_2050
rem_g1015.t1	ste_g35233.t1	rem_249_3	ste_33_193
rem_g1016.t1	ste_g35234.t1	rem_249_4	ste_33_194

                                            		Du Kang 2019-5-6
--------------------------------
EOF
	# about how to use the script

@ARGV or die $usage;

$egg=5;
$gap=15;
$ohno=0;
foreach $i (0..@ARGV-1) {
	$egg=$ARGV[$i+1] if $ARGV[$i] eq "-egg";
	$gap=$ARGV[$i+1] if $ARGV[$i] eq "-gap";
	$ohno=1 if $ARGV[$i]=~/-ohno/;
}
	# set parameters

$log="$ARGV[0].synteny.log";
unlink $log if -e $log;
open LOG, ">>$log" or die $!;
	# logs output here

$perl=q#perl -lane 'print "$F[-2]\t$_"'#;
$in="cat $ARGV[0] |$perl |sort -V |cut -f2- |";
	# sort the second last column

open IN, $in or die $!;
@lines=<IN>;
$ini=0;
foreach $i (0..@lines-1){
	$l=$i+1;
	@col0=split /\s+/, $lines[$i];

	$col0[-1]=~/.*?_(.*?)_(.*?)(_|$)/;
	$chain0="$1";
	$anchor_chain=$chain0 if $i==0;
	$pos0=$2;

	$col0[-2]=~/.*?_(.*?)_/;
	$chr0="$1";

	if ($ini==1 and $chain0 eq $anchor_chain and abs($pos0-$pre_pos0)>$gap) {
		$ini=0;
		print LOG "line $l: gap size large then $gap, melt the chain: $pos0\t$pre_pos0\n\n";

	} elsif ($ini==1 and $chain0 ne $anchor_chain){
		print LOG "line $l: chain interupted, start the melting mechanism\n";

		$heating=1;
		$reconecte=0;
		foreach $n ($i+1..@lines-1){
			$l1=$n+1;
			@col1=split /\s+/, $lines[$n];

			$col1[-1]=~/.*?_(.*?)_(.*?)(_|$)/;
			$chain1="$1";
			$pos1=$2;

			if($chain1 ne $anchor_chain){
				$heating++;		# heat up to melt the chain. When heat achieve the $gap, melt the chain
			}elsif(abs($pos1-$pre_pos0) >$gap){
				$heating=$gap;
			}else{
				$reconecte=1;
			}

			if ($reconecte==1){
				print LOG "line $l: chain reconected on line $l1, stop the melting mechanism: $pos1\t$pre_pos0\n\n";
				last;
			}elsif ($heating==$gap) {
				$ini=0;
				print LOG "line $l: gap $gap achived on line $l1, melt the chain\n\n";
				last;
			} 

		}
	} 

	if (($ohno==1 and $ini==0 and $chain0 ne $chr0) or ($ohno==0 and $ini==0)) {
		print LOG "line $l: start the initiating processing\n";

		$heating=0;
		$count=1,
		$anchor_chain=$chain0;
		$pre_pos1=$pos0;
		foreach $n ($i+1..@lines-1){
			$l1=$n+1;
			@col1=split /\s+/, $lines[$n];

			$col1[-1]=~/.*?_(.*?)_(.*?)(_|$)/;
			$chain1="$1";
			$pos1=$2;
			
			if($chain1 ne $anchor_chain){
				$heating++;
			}elsif(abs($pos1-$pre_pos1) >$gap){
				$heating=$gap;
			}else{
				$count++;
			}

			if ($count==$egg){
				$ini=1;
				print LOG "line $l: initiated\n\n";
				last;
			} elsif ($heating==$gap) {
				print LOG "line $l: melted the initiating processing on line $l1, initiating failed\n\n";
				last;
			}

			$pre_pos1=$pos1 if $chain1 eq $anchor_chain;
		}
	}

	if ($ohno==1){
		print $lines[$i] if $ini==1 and $chain0 ne $chr0 and $chain0 eq $anchor_chain;	# need the two chromosome not to be the same
	}else{
		print $lines[$i] if $ini==1 and $chain0 eq $anchor_chain;		
	}

	$pre_pos0=$pos0 if $chain0 eq $anchor_chain;
}
