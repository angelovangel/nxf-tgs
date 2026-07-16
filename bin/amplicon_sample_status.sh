#!/usr/bin/env bash

# Generate amplicon_sample_status.txt for wf-amplicon runs.
# Pass/fail is determined by presence of the sample alias in all-consensus-seqs.fasta.
# The file is always written, even when the FASTA is absent (all samples failed).
#
# arg1: path to per-user samplesheet CSV (header: alias,barcode,approx_size)
# arg2: path to 02-assembly dir (default: 02-assembly)
#
# Output: 02-assembly/amplicon_sample_status.txt  (sample,pass_fail,length)

SAMPLESHEET="$1"
ASSEMBLY_DIR="${2:-02-assembly}"
FASTA="${ASSEMBLY_DIR}/all-consensus-seqs.fasta"
OUT="${ASSEMBLY_DIR}/amplicon_sample_status.txt"

# Build a lookup file: "alias,pass,<length>" for each sequence in the FASTA
> _pass.txt
if [ -f "$FASTA" ]; then
    awk '/^>/{
            if (seq) print name",pass,"length(seq)
            sub(/^>/, ""); name=$1; seq=""; next
         }
         { seq = seq $0 }
         END { if (seq) print name",pass,"length(seq) }' \
        "$FASTA" > _pass.txt
fi

# Write output - iterate over all aliases in the samplesheet as source of truth
printf "sample,pass_fail,length\n" > "$OUT"
while IFS= read -r s; do
    line=$(grep "^${s}," _pass.txt || true)
    if [ -n "$line" ]; then
        echo "$line" >> "$OUT"
    else
        echo "${s},fail,N/A" >> "$OUT"
    fi
done < <(tail -n +2 "$SAMPLESHEET" | cut -d',' -f1)

echo "Wrote amplicon sample status to ${OUT}"
