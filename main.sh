#!/usr/bin/env bash
set -e

# ==============================================================================
# Bash Orchestration of NXF-TGS Pipeline
# ==============================================================================

function show_help {
    echo "Usage: ./main.sh [OPTIONS]"
    echo "  --fastq <path>            : path to raw fastq_pass data"
    echo "  --samplesheet <path>      : path to csv or excel samplesheet"
    echo "  --pipeline <name>         : epi2me workflow (wf-clone-validation, wf-bacterial-genomes, wf-amplicon, report-only)"
    echo "  --assembly_args <args>    : additional arguments for assembly workflow"
    echo "  --assembly_profile <prof> : profile to use (standard, singularity, test) default: standard"
    echo "  --outdir <dir>            : output directory, default: output"
    echo "  --cpus <num>              : number of cpus to use, default: 4"
    echo "  --help                    : show this help message"
    exit 0
}

# Default parameters
PIPELINE="report-only"
ASSEMBLY_PROFILE="standard"
OUTDIR="output"
CPUS=4
ASSEMBLY_ARGS=""

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --fastq) FASTQ="$2"; shift ;;
        --samplesheet) SAMPLESHEET="$2"; shift ;;
        --pipeline) PIPELINE="$2"; shift ;;
        --assembly_args) ASSEMBLY_ARGS="$2"; shift ;;
        --assembly_profile) ASSEMBLY_PROFILE="$2"; shift ;;
        --outdir) OUTDIR="$2"; shift ;;
        --cpus) CPUS="$2"; shift ;;
        --help) show_help ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [[ -z "$FASTQ" || -z "$SAMPLESHEET" ]]; then
    echo "Error: --fastq and --samplesheet are required."
    show_help
fi

# Determine pipeline revision
PIPELINE_VER="N/A"
case "$PIPELINE" in
    "wf-clone-validation") PIPELINE_VER="v1.8.4" ;;
    "wf-amplicon") PIPELINE_VER="v1.2.2" ;;
    "wf-bacterial-genomes") PIPELINE_VER="v2.0.2" ;;
esac

# Ensure paths are absolute for Docker mounts
[[ "$FASTQ" != /* ]] && FASTQ="$(pwd)/$FASTQ"
[[ "$SAMPLESHEET" != /* ]] && SAMPLESHEET="$(pwd)/$SAMPLESHEET"
[[ "$OUTDIR" != /* ]] && OUTDIR="$(pwd)/$OUTDIR"

# Tools container
DOCKER_IMAGE="aangeloo/nxf-tgs:latest"
IGV_DOCKER_IMAGE="aangeloo/igv-reports-image:latest"

# Setup directories
mkdir -p "$OUTDIR"
TMPDIR=$(mktemp -d)
echo "Using temporary directory: $TMPDIR"
# trap 'rm -rf -- "$TMPDIR"' EXIT

# 1. Read and Validate Samplesheet
SAMPLESHEET_CSV="$SAMPLESHEET"
if [[ "$SAMPLESHEET" == *.xlsx ]]; then
    echo "Converting Excel to CSV..."
    docker run --rm -v "$(pwd):/work" -v "$FASTQ":"$FASTQ" -v "$(dirname "$SAMPLESHEET")":"$(dirname "$SAMPLESHEET")" -v "$OUTDIR":"$OUTDIR" -w "$OUTDIR" "$DOCKER_IMAGE" bash -c "export PATH=\"/work/bin:\$PATH\"; convert_excel.R '$SAMPLESHEET'"
    SAMPLESHEET_CSV="$OUTDIR/$(basename "${SAMPLESHEET%.xlsx}.csv")"
fi

echo "Validating Samplesheet..."
docker run --rm -v "$(pwd):/work" -v "$FASTQ":"$FASTQ" -v "$(dirname "$SAMPLESHEET")":"$(dirname "$SAMPLESHEET")" -v "$OUTDIR":"$OUTDIR" -w "$OUTDIR" "$DOCKER_IMAGE" bash -c "
    export PATH=\"/work/bin:\$PATH\"
    validate_samplesheet.R '$SAMPLESHEET_CSV' '$FASTQ'
    if [ '$PIPELINE' != 'wf-bacterial-genome' ]; then
        get_maxbin.sh samplesheet-validated.csv '$FASTQ'
    else
        mv samplesheet-validated.csv 00-samplesheet-validated.csv
    fi
"

# Parse users and samples
USERS=$(awk -F, 'NR>1 && $0 ~ /OK/ {print $1}' "$OUTDIR/00-samplesheet-validated.csv" | sort | uniq | tr -d '\r')

# 2. Merge Reads & Generate Reports per User
for USER in $USERS; do
    echo "Processing user: $USER"
    USER_OUT="$OUTDIR/$USER"
    mkdir -p "$USER_OUT/01-fastq"

    # Get user specific samples
    awk -F, -v u="$USER" 'NR>1 && $1==u && $0~/OK/ {print $2","$4","$3}' "$OUTDIR/00-samplesheet-validated.csv" > "$TMPDIR/${USER}_samples.csv"
    
    # Merge fastq files
    while IFS=, read -r sample barcode size; do
        if [[ "$barcode" == barcode* ]]; then
            echo "  Merging reads for $sample ($barcode)..."
            docker run --rm -v "$(pwd):/work" -v "$FASTQ":"$FASTQ" -v "$(dirname "$SAMPLESHEET")":"$(dirname "$SAMPLESHEET")" -v "$OUTDIR":"$OUTDIR" -w /work "$DOCKER_IMAGE" bash -c "
                export PATH=\"/work/bin:\$PATH\"
                shopt -s extglob
                find '$FASTQ/${barcode}/' -type f ! -name '*.gz' -exec pigz {} \;
                cat '$FASTQ/${barcode}/'@(*.fastq|*.fq).gz > '${USER_OUT}/01-fastq/${sample}.fastq.gz'
            "
        fi
    done < <(tr -d '\r' < "$TMPDIR/${USER}_samples.csv")

    # Generate Reports
    echo "  Generating FASTQ reports..."
    docker run --rm -v "$(pwd):/work" -v "$FASTQ":"$FASTQ" -v "$(dirname "$SAMPLESHEET")":"$(dirname "$SAMPLESHEET")" -v "$OUTDIR":"$OUTDIR" -w /work "$DOCKER_IMAGE" bash -c "
        export PATH=\"/work/bin:\$PATH\"
        FASTQ_FILES=\$(ls ${USER_OUT}/01-fastq/*.fastq.gz | tr '\n' ' ')
        echo -e 'file\treads\tbases\tn_bases\tmin_len\tmax_len\tN50\tGC_percent\tQ20_percent' > '${USER_OUT}/01-${USER}-faster-report.tsv'
        faster2 -ts \$FASTQ_FILES >> '${USER_OUT}/01-${USER}-faster-report.tsv'
        
        # HTML Report wrapper
        cd '${USER_OUT}' && faster-report.R -p '01-fastq' \
            --outfile '01-${USER}-faster-report' \
            --user '${USER}' --rundate 'NA' --flowcell 'NA' --basecall 'NA'
    "
done

# 3. Assembly
if [[ "$PIPELINE" != "report-only" ]]; then
    echo "Starting Assembly with pipeline: $PIPELINE"
    
    # NOTE: The epi2me-labs pipelines themselves are Nextflow pipelines. 
    # If 'no nextflow' strictly means avoiding Nextflow entirely, 
    # this step would require manually implementing the entire epi2me pipeline in Bash.
    # We use Nextflow here specifically to run the epi2me-labs pipeline, as it was treated as a black-box in main.nf.
    
    for USER in $USERS; do
        echo "  Running assembly for $USER..."
        USER_OUT="$OUTDIR/$USER"
        
        # Create user specific samplesheet
        echo "alias,barcode,approx_size" > "$TMPDIR/${USER}_samplesheet.csv"
        awk -F, -v u="$USER" 'NR>1 && $1==u && $0~/OK/ {print $2","$4","$3}' "$OUTDIR/00-samplesheet-validated.csv" >> "$TMPDIR/${USER}_samplesheet.csv"

        mkdir -p "$USER_OUT/02-assembly"
        
        REV_ARG=""
        if [[ "$PIPELINE_VER" != "N/A" ]]; then
            REV_ARG="-r $PIPELINE_VER"
        fi
        
        (
            cd "$USER_OUT/02-assembly"
            
            echo "executor { name = 'local'; cpus = ${CPUS} }" > child.config
            
            NXF_VER="24.10.9" nextflow run epi2me-labs/${PIPELINE} \
                --fastq "$FASTQ" \
                --sample_sheet "$TMPDIR/${USER}_samplesheet.csv" \
                --out_dir "$USER_OUT/02-assembly" \
                $ASSEMBLY_ARGS \
                -profile $ASSEMBLY_PROFILE \
                -c child.config \
                $REV_ARG
        )

        if [[ "$PIPELINE" == "wf-clone-validation" ]]; then
            # Fix annotations.bed
            docker run --rm -v "$(pwd):/work" -v "$FASTQ":"$FASTQ" -v "$(dirname "$SAMPLESHEET")":"$(dirname "$SAMPLESHEET")" -v "$OUTDIR":"$OUTDIR" -w /work "$DOCKER_IMAGE" bash -c "
                export PATH=\"/work/bin:\$PATH\"
                cd '${USER_OUT}/02-assembly' && make_bed.R feature_table.txt
            "
        fi
        
        # Ensure assembly stats exists
        if ! ls ${USER_OUT}/02-assembly/*.assembly_stats.tsv 1> /dev/null 2>&1; then
            touch "${USER_OUT}/02-assembly/empty.assembly_stats.tsv"
            echo -e "sample_name\tmean_quality" >> "${USER_OUT}/02-assembly/empty.assembly_stats.tsv"
        fi
    done
fi

# 4. Post-Assembly Workflows (Mapping, IGV, Summary)
if [[ "$PIPELINE" == "wf-clone-validation" ]]; then
    for USER in $USERS; do
        USER_OUT="$OUTDIR/$USER"
        
        # Sample Status
        docker run --rm -v "$(pwd):/work" -v "$FASTQ":"$FASTQ" -v "$(dirname "$SAMPLESHEET")":"$(dirname "$SAMPLESHEET")" -v "$OUTDIR":"$OUTDIR" -w /work "$DOCKER_IMAGE" bash -c "
            export PATH=\"/work/bin:\$PATH\"
            cd '${USER_OUT}/02-assembly'
            sample_status.R '${USER}' 'sample_status.txt' '${OUTDIR}/00-samplesheet-validated.csv'
            mv sample-status.csv '../00-${USER}-sample_status.csv' || true
        "
        
        # Mapping for each sample
        mkdir -p "$USER_OUT/03-mapping"
        mkdir -p "$USER_OUT/04-igv-reports"
        
        for FASTA in ${USER_OUT}/02-assembly/*.final.fasta; do
            [[ -f "$FASTA" ]] || continue
            SAMPLE=$(basename "$FASTA" .final.fasta)
            FASTQ_FILE="${USER_OUT}/01-fastq/${SAMPLE}.fastq.gz"
            
            echo "  Mapping $SAMPLE for $USER..."
            docker run --rm -v "$(pwd):/work" -v "$FASTQ":"$FASTQ" -v "$(dirname "$SAMPLESHEET")":"$(dirname "$SAMPLESHEET")" -v "$OUTDIR":"$OUTDIR" -w /work "$DOCKER_IMAGE" bash -c "
                export PATH=\"/work/bin:\$PATH\"
                minimap2 -ax lr:hq '$FASTA' '$FASTQ_FILE' > '${USER_OUT}/03-mapping/mapping.sam'
                samtools view -S -b -T '$FASTA' '${USER_OUT}/03-mapping/mapping.sam' | \
                samtools sort -o '${USER_OUT}/03-mapping/${SAMPLE}.bam' -
                samtools index '${USER_OUT}/03-mapping/${SAMPLE}.bam'
                rm '${USER_OUT}/03-mapping/mapping.sam'
                
                perbase base-depth --threads 4 '${USER_OUT}/03-mapping/${SAMPLE}.bam' -F 260 > '${USER_OUT}/03-mapping/${SAMPLE}.perbase.tsv'
                perpos_freq.sh '${USER_OUT}/03-mapping/${SAMPLE}.perbase.tsv' > '${USER_OUT}/03-mapping/${SAMPLE}.problems.tsv'
            "
            
            # IGV Reports
            echo "  Generating IGV report for $SAMPLE..."
            docker run --rm -v "$(pwd):/work" -v "$FASTQ":"$FASTQ" -v "$(dirname "$SAMPLESHEET")":"$(dirname "$SAMPLESHEET")" -v "$OUTDIR":"$OUTDIR" -w /work "$IGV_DOCKER_IMAGE" bash -c "
                export PATH=\"/work/bin:\$PATH\"
                samtools faidx '$FASTA'
                LEN=\$(awk '{print \$2}' '${FASTA}.fai' | head -n 1)
                HEADER=\$(grep '>' '$FASTA' | cut -c 2-)
                
                echo -e \"\$HEADER\t0\t\$LEN\tHET\" > '${USER_OUT}/04-igv-reports/bedfile.bed'
                
                create_report '${USER_OUT}/04-igv-reports/bedfile.bed' \
                    --fasta '$FASTA' \
                    --tracks '${USER_OUT}/02-assembly/${SAMPLE}.annotations2.bed' '${USER_OUT}/03-mapping/${SAMPLE}.bam' \
                    --output '${USER_OUT}/04-igv-reports/${SAMPLE}.igvreport.html' \
                    --flanking 200
            "
        done
        
        # Mapping Summary
        docker run --rm -v "$(pwd):/work" -v "$FASTQ":"$FASTQ" -v "$(dirname "$SAMPLESHEET")":"$(dirname "$SAMPLESHEET")" -v "$OUTDIR":"$OUTDIR" -w /work "$DOCKER_IMAGE" bash -c "
            export PATH=\"/work/bin:\$PATH\"
            cd '${USER_OUT}/03-mapping'
            mapping_counts.R '.*mapping-counts.txt' '${USER}'
            mv mapping-summary.csv '../00-${USER}-mapping-summary.csv' || true
        "
    done
    
    # Global Summaries
    docker run --rm -v "$(pwd):/work" -v "$FASTQ":"$FASTQ" -v "$(dirname "$SAMPLESHEET")":"$(dirname "$SAMPLESHEET")" -v "$OUTDIR":"$OUTDIR" -w /work "$DOCKER_IMAGE" bash -c "
        export PATH=\"/work/bin:\$PATH\"
        mkdir -p '${OUTDIR}/summary_temp'
        cp '${OUTDIR}'/*/00-*-sample_status.csv '${OUTDIR}/summary_temp/' 2>/dev/null || true
        cp '${OUTDIR}'/*/00-*-mapping-summary.csv '${OUTDIR}/summary_temp/' 2>/dev/null || true
        cd '${OUTDIR}/summary_temp'
        sample_summary.R '.*sample_status.csv' 'bash-run' 'HEAD' '$PIPELINE' '$PIPELINE_VER' || true
        mapping_summary.R '.*mapping-summary.csv' 'bash-run' || true
        mv *.html .. 2>/dev/null || true
        cd .. && rm -rf summary_temp
    "
fi

echo "Pipeline finished successfully!"
