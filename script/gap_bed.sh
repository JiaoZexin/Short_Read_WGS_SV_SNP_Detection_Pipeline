#!/bin/bash

# Usage: bash gap_bed.sh ref.fa output.txt
ref_fa=$1
output_txt=$2

perl -ne 'chomp;
    if(/^>(.*)/){
        $head = $1; $i = 0; next
    };
    @a = split("",$_);
    foreach(@a){
        $i++;
        if($_ eq "N" && $s==0){
            $z = $i-1;
            print "$head\t$z";
            $s=1
        }
        elsif($s==1 && $_ ne "N"){
            $j=$i-1;
            print "\t$j\n";
            $s=0
        }
    }' $ref_fa > $output_txt
