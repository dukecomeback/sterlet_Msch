#!/usr/bin/perl -w

my $usage=<<EOF;
------------------------------------------------------------------------
Feed me with gff output from genewise (-gff), I'll give you non-overlapped best results.
This is what I'm gonna do:
1) sort the genewise output file in order of scaffold name (strand included) and decreasing genewise score.
2) go though the file line by line and print if it have no overlap with printed results.
Usage: perl $0 genewise.gtf >genewise.gtf.best
										Du Kang 2017-07-14
EOF

$ARGV[0] or die $usage;
`cp $ARGV[0] tmp20170713`;
system (q#cat tmp20170713 | perl -ne 's/\n/&/;s/(\/\/)/\n$1/; print "\/\/&" if $.==1;print'|sort -k1,1 -k7,7 -k6nr,6 >tmp20170714#);


use Set::IntervalTree;
$flag='initial';
open IN, "tmp20170714" or die $!;
while (<IN>) {
	next if /^\/\/\&$/;
	chomp;
	@_=split; 
	$_[0]=~/\&(.*)/;
	$scaffold=$1;
	if($_[6] eq '+'){
		$current=$scaffold."for";
		$low=$_[3];
		$high=$_[4];
	}else{
		$current=$scaffold."rev";
                $low=$_[4];
                $high=$_[3];
	}
	if ($current ne $flag) {
		$flag=$current;
		$tree = Set::IntervalTree->new;
		$object=1;
		$tree->insert($object,$low,$high);
		s/&/\n/g; 
		print;
	}else {
		$results = $tree->fetch($low,$high);
		if (! scalar(@$results)){
			$object++;
			$tree->insert($object,$low,$high);
			s/&/\n/g; 
			print;	
		}
	}
}
print "\/\/";

`rm tmp20170713 tmp20170714`;
