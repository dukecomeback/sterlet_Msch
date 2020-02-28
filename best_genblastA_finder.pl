#!/usr/bin/perl -w

my $usage=<<EOF;
------------------------------------------------------------------------
Feed me with genblastA output, I'll give you non-overlapped best hit.
This is what I'm gonna do:
1) sort the file in order of scaffold name (strand included) and decreasing score.
2) go though the file line by line and print it if it have no overlap with printed results.
Usage: perl $0 genblastA.out >genblastA.out
                                                                                Du Kang 2017-09-23
EOF

$ARGV[0] or die $usage;
`cp $ARGV[0] tmp20170923`;
system (q#perl -i -ne 'chomp;if(/rank:/){print "\n$_"} if(/HSP.*query:\((.*)-(.*)\)/){print "\t$1\t$2"}' tmp20170923#);
system (q#perl -i -ne 'chomp;@F=split /\t/, $_; print "$F[0]\t";use List::Util qw/max min/; shift @F; $max=max @F; $min=min @F;print "$min\t$max\n"' tmp20170923#);
system (q#perl -i -lane '/.*\|(.*):(\d+)\.\.(\d+)\|(.*)\|.*\|.*:(.*)\|.*/;print "$1\t$4\t$2\t$3\t$5\t$_"' tmp20170923#);
system (q#sort -k1,1 -k2,2 -k5,5nr tmp20170923 -o tmp20170923#);


use Set::IntervalTree;
$flag='initial';
open IN, "tmp20170923" or die $!;
while (<IN>) {
	next if /^\s+$/;
        @_=split;
        $current= $_[1] eq "+" ? $_[0]."for" : $_[0]."rev";
        $low=$_[2];
        $high=$_[3];
	s/$_[0]\t[+|-]\t$_[2]\t$_[3]\t$_[4]\t//;
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

`rm tmp20170923`;

