#! /usr/bin/env bash

# requires fasterplot to be installed and in env
# arg1 is path to csv file with user,sample,barcode,dna_size
# arg2 is path to fastq_pass
# output is csv supplemented with maxbin [integer]

barcode_idx=$(head -1 $1 | sed 's/,/\n/g' | nl | grep 'barcode' | cut -f 1)

echo "$(head -n 1 $1),obs_size" > samplesheet-validated.csv

while IFS="," read line; do
    barcode=$(echo $line | cut -f $barcode_idx -d,)
    obs_size=$(cat $2/$barcode/*.fastq.gz | seqkit seq -M 49999 -g | fasterplot -l - |  grep "# maxbin:" | cut -f2 )
    
    echo -e "$line,$obs_size" >> samplesheet-validated.csv


done < <(tail -n +2 $1)