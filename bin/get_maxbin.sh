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
echo "Barcode index: $barcode_idx"

# 3. Create the new header for the output file
echo "$HEADER,obs_size,nreads" > "$OUTPUT_FILE"

# 4. Process data rows
# Read data rows (tail -n +2) and process them
while IFS="," read line; do
    
    # Extract the barcode using the determined index
    barcode=$(echo "$line" | cut -f "$barcode_idx" -d, | tr -d '"')
    echo "Working on $FASTQ_DIR/$barcode/"
    
    # Detect BAM vs fastq.gz and compute obs_size / nreads accordingly
    if [ ! -d "$FASTQ_DIR/$barcode" ]; then
        echo "Warning: directory $FASTQ_DIR/$barcode does not exist, skipping" >&2
        obs_size=""
        nreads=0
    elif ls "$FASTQ_DIR/$barcode"/*.bam 1>/dev/null 2>&1; then
        # BAM mode: many BAM files per barcode are expected, so merge with
        # samtools merge (-u = uncompressed, faster since it's a pipe) and
        # stream straight into samtools fastq via stdin - no tmp file required.
        bam_files=("$FASTQ_DIR/$barcode"/*.bam)

        read_count=$(samtools merge -u -o - "${bam_files[@]}" | samtools view -c -)
        if [ "$read_count" -eq 0 ]; then
            echo "Warning: no reads found for barcode $barcode, skipping conversion" >&2
            obs_size=""
            nreads=0
        else
            obs_size=$(samtools merge -u -o - "${bam_files[@]}" | samtools fastq - 2>/dev/null | \
                       seqkit seq -M 49499 -g | \
                       fasterplot -l - | \
                       grep "# maxbin:" | \
                       cut -f2 | \
                       tr -d ' ')
            nreads=$(samtools merge -u -o - "${bam_files[@]}" | samtools fastq - 2>/dev/null | faster2 -ts - | cut -f 2 | tr -d ' ')
        fi
    elif ls "$FASTQ_DIR/$barcode"/*.fastq.gz 1>/dev/null 2>&1; then
        fq_files=("$FASTQ_DIR/$barcode"/*.fastq.gz)

        read_count=$(( $(gunzip -c "${fq_files[@]}" | wc -l) / 4 ))
        if [ "$read_count" -eq 0 ]; then
            echo "Warning: no reads found for barcode $barcode, skipping conversion" >&2
            obs_size=""
            nreads=0
        else
            obs_size=$(cat "${fq_files[@]}" | \
                       seqkit seq -M 49499 -g | \
                       fasterplot -l - | \
                       grep "# maxbin:" | \
                       cut -f2 | \
                       tr -d ' ') # Remove potential leading/trailing spaces
            nreads=$(cat "${fq_files[@]}" | faster2 -ts - | cut -f 2 | tr -d ' ')
        fi
    else
        # Directory exists but has no BAM or fastq.gz files in it
        echo "Warning: no BAM or fastq.gz files found in $FASTQ_DIR/$barcode, skipping" >&2
        obs_size=""
        nreads=0
    fi
    
    # Append the original line and the calculated obs_size to the output file
    echo "$line,$obs_size,$nreads" >> "$OUTPUT_FILE"

done < <(tail -n +2 "$INPUT_FILE")

echo "Successfully generated $OUTPUT_FILE with the new 'obs_size' column."