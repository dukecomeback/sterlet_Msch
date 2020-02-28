#!/usr/bin/perl -w

my $usage=<<EOF;
------------------------------------------------------------------------
Feed me with gtf output from exonerate, I'll give you non-overlapped best results.
This is what I'm gonna do:
1) sort the exonerate output file in order of scaffold name (strand included) and decreasing exonerate score.
2) go though the file line by line and print if it have no overlap with printed results.
Usage: perl $0 exonerate.out > exonerate.out.best
                                                                                Du Kang 2017-07-14
EOF

$ARGV[0] or die $usage;
`cp $ARGV[0] tmp20170713`;
system (q#cat tmp20170713 | perl -ne 'next if /Command|Hostname/;s/\n/&/;s/(\# --- START OF GFF DUMP)/\n$1/; print'|perl -aF'\t' -lne 'next if $.==1;/attributes&\#&(.*?)\s*exonerate/;print "$1\t$F[5]\t".$_'|sort -t $'\t' -k1,1 -k9,9 -k2nr,2 >tmp20170714#);

use Set::IntervalTree;
$flag='initial';
open IN, "tmp20170714" or die $!;
while (<IN>) {
        chomp;
        @_=split /\t/;
        $current= $_[8] eq "+"? $_[0]."for" : $_[0]."rev";
        $low=$_[5];
        $high=$_[6];
        if ($current ne $flag) {
                $flag=$current;
                $tree = Set::IntervalTree->new;
                $object=1;
                $tree->insert($object,$low,$high);
		s/$_[0]\t$_[1]\t//;
                s/&/\n/g;
                print;
        }else {
                $results = $tree->fetch($low,$high);
                if (! scalar(@$results)){
                        $object++;
                        $tree->insert($object,$low,$high);
			s/$_[0]\t$_[1]\t//;
                        s/&/\n/g;
                        print;
                }
        }
}

`rm tmp20170713 tmp20170714`;
