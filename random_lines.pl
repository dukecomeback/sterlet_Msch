#!/usr/bin/perl -w

my $usage=<<EOF;
--------------------------------
I randomly sample n(10) line from a file

Usage: perl $0 infile sample_size >infile.sample
                                            		Du Kang 2019-06-04
--------------------------------
EOF
    # about how to use the script

@ARGV==2 or die $usage;

$size=$ARGV[1];

`wc $ARGV[0]` =~ /(\d+?)\s+/;
$max=$1-1;

@keep=sample_no_return($max, $size);

open IN, $ARGV[0] or die $!;
while (<IN>) {
	print if  grep { $_ eq $. } @keep;
}



sub sample_no_return {
	# Usage: sample_no_return (sampling_range, sample_size)
	my @sample=();
	my $range=$_[0];
	my $size=$_[1];
	while (@sample <$size) {
		my $random=int(rand($range));
		push @sample, $random unless grep { $_ eq $random } @sample;
	}
	return @sample;
}
