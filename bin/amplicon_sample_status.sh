#!/usr/bin/env bash

# Generate amplicon_sample_status.txt for wf-amplicon runs.
# Pass/fail is determined by presence of the sample alias in all-consensus-seqs.fasta.
# The file is always written, even when the FASTA is absent (all samples failed).
#
# arg1: path to per-user samplesheet CSV (header: alias,barcode,approx_size)
# arg2: path to 02-assembly dir (default: 02-assembly)
#
# Output: 02-assembly/amplicon_sample_status.txt  (sample,pass_fail,length)
#         02-assembly/amplicon.assembly_stats.tsv   (sample_name,mean_quality) - as fastcat output in wf-clone-validation

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

# Mean consensus quality per sample from <sample>/consensus/consensus.fastq
# as in fastcat: mean of per-base error probabilities, converted back to phred
STATS="${ASSEMBLY_DIR}/amplicon.assembly_stats.tsv"
printf "sample_name\tmean_quality\n" > "$STATS"
for fq in "${ASSEMBLY_DIR}"/*/consensus/consensus.fastq; do
    [ -s "$fq" ] || continue
    sample=$(basename "$(dirname "$(dirname "$fq")")")
    awk -v s="$sample" '
        BEGIN { for (i = 33; i < 127; i++) q[sprintf("%c", i)] = i - 33 }
        NR % 4 == 0 { for (i = 1; i <= length($0); i++) { e += 10^(-q[substr($0, i, 1)] / 10); n++ } }
        END { if (n > 0) printf "%s\t%.2f\n", s, -10 * log(e / n) / log(10) }
    ' "$fq" >> "$STATS"
done

echo "Wrote amplicon assembly quality to ${STATS}"
