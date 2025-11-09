#! /usr/bin/env bash

# requires fasterplot to be installed and in env
# arg1 is path to csv file with user,sample,barcode,dna_size
# arg2 is path to fastq_pass
# output is csv supplemented with maxbin [integer] and nreads [integer]

INPUT_FILE="$1"
FASTQ_DIR="$2"
OUTPUT_FILE="00-samplesheet-validated.csv"

# 1. Check if 'obs_size' column already exists
HEADER=$(head -n 1 "$INPUT_FILE")
if echo "$HEADER" | grep -q '\bobs_size\b'; then
    # If 'obs_size' is present, just return the filename 
    echo "The column 'obs_size' is already present in the input file."
    # We copy the original file to the output name, assuming it's the valid one.
    cp "$INPUT_FILE" "$OUTPUT_FILE"
    exit 0
fi

# If 'obs_size' is NOT present, proceed with calculation and addition

# 2. Find the index of the 'barcode' column
barcode_idx=$(echo "$HEADER" | sed 's/,/\n/g' | nl | grep 'barcode' | awk '{print $1}')

# 3. Create the new header for the output file
echo "$HEADER,obs_size,nreads" > "$OUTPUT_FILE"

# 4. Process data rows
# Read data rows (tail -n +2) and process them
while IFS="," read line; do
    
    # Extract the barcode using the determined index
    barcode=$(echo "$line" | cut -f "$barcode_idx" -d,)
    
    # Calculate obs_size (maxbin)
    obs_size=$(cat "$FASTQ_DIR/$barcode"/*.fastq.gz | \
               seqkit seq -M 49499 -g | \
               fasterplot -l - | \
               grep "# maxbin:" | \
               cut -f2 | \
               tr -d ' ') # Remove potential leading/trailing spaces
    nreads=$(cat "$FASTQ_DIR/$barcode"/*.fastq.gz | faster2 -ts - | cut -f 2 | tr -d ' ')
    
    # Append the original line and the calculated obs_size to the output file
    echo "$line,$obs_size,$nreads" >> "$OUTPUT_FILE"

done < <(tail -n +2 "$INPUT_FILE")

echo "Successfully generated $OUTPUT_FILE with the new 'obs_size' column."