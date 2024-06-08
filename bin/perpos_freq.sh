#! /usr/bin/env bash

# calculate per base allele frequency from the output of `perbase base-depth`
# then filter for potential heterozygous positions

echo "POS DEPTH A C G T INS DEL"
awk 'NR > 1 { print $2, $3, $4/$3, $5/$3, $6/$3, $7/$3, $9/$3, $10/$3 }' $1 | \
# get suspicious positions only (heterozygous)
awk 'NF>1 {for(i=3; i <= NF; i+=1) if ($i > 0.3 && $i < 0.7) print }' | uniq