#!/usr/bin/perl -w
my $usage=<<E;
----------------------------------------
This is for concatenated aligned fastas.
Usage: perl $0 aa.mus.tri.all(obtained by cat all fasta files) >aa.align.concatenated
										Du Kang 2017-1-30
---------------------------------------
E
open IN, $ARGV[0] or die $usage;
while (<IN>) {
	chomp;
	if (/>/) {
		s/>//;
# in 180809		/([a-zA-Z_]*)\d*/;
		/(.*?)_/;
		$name = $1;
	}else{
		$seq{$name} .= $_;
	}
}
foreach $key (keys %seq) {
	print ">$key\n$seq{$key}\n";
}
