
if (params.help) {
    helpMessage()
    exit(0)
}

def helpMessage() {
    log.info """\
    ========================================================================================================================
    NXF - TGS ONT PIPELINE
    process (per user) raw fastq_pass folder - merge/rename, generate report, assembly (plasmid, amplicon, bacterial genome)
    ========================================================================================================================
    Usage:
    -----------------------------------
    fastq         : path to raw fastq_pass data
    samplesheet   : path to csv or excel with (at least) columns sample, barcode, user
    pipeline      : epi2me workflow to use - can be wf-clone-validation, wf-bacterial-genomes, wf-amplicon, report-only
    assembly_args : additional command-line arguments passed to the assembly workflow
    outdir        : where to save results, default is 'output'
    """
    .stripIndent(true)
}

log.info """\
    ========================================================================================================================
    NXF - TGS ONT PIPELINE
    process (per user) raw fastq_pass folder - merge/rename, generate reports, assembly (plasmid, amplicon, bacterial genome)
    ========================================================================================================================
    fastq           : ${params.fastq}
    samplesheet     : ${params.samplesheet}
    pipeline        : ${params.pipeline}
    assembly_args   : ${params.assembly_args}
    outdir          : ${params.outdir}
    """
    .stripIndent(true)

fastq_pass_ch = Channel.fromPath(params.fastq, type: 'dir', checkIfExists: true)
samplesheet_ch = Channel.fromPath(params.samplesheet)
wf_versions = Channel.from(['wf-clone-validation', 'v1.5.0'], ['wf-amplicon', 'v1.1.3'], ['wf-bacterial-genomes', 'v1.4.0'])
wf_ver = Channel.from(params.pipeline).join(wf_versions)

// takes in csv, checks for duplicate barcodes, unique sample names per user
// emits *-checked.csv if all ok, exits with error if not
// also add observed peak size (as seen by fasterplot)
process VALIDATE_SAMPLESHEET {
    container 'aangeloo/nxf-tgs:latest'
    publishDir "$params.outdir", mode: 'copy', pattern: '00-samplesheet-validated.csv'

    input: 
    path(csv)
    path(fastq_pass)

    output:
    path("00-samplesheet-validated.csv")

    script:
    """
    validate_samplesheet.R $csv $fastq_pass
    get_maxbin.sh samplesheet-validated.csv $fastq_pass

    """
}

process READEXCEL {
    container 'aangeloo/nxf-tgs:latest'

    input:
    path(excelfile)

    output:
    path("*.csv")

    script:
    """
    convert_excel.R $excelfile
    """
}

process MERGE_READS {
    container 'aangeloo/nxf-tgs:latest'
    tag "$user - $samplename"
    errorStrategy 'ignore' //because some barcodes defined in the samplesheet might be missing in the data
    publishDir "$params.outdir/$user/01-fastq", mode: 'copy', pattern: '*.fastq.gz'

    input:
    tuple val(samplename), val(barcode), val(user), path(mypath)
    
    output: 
    tuple val(user), path('*.fastq.gz'), emit: merged_fastq_ch
    
    script:
    """
    #https://www.gnu.org/savannah-checkouts/gnu/bash/manual/bash.html#Pattern-Matching
    shopt -s extglob

    find ${mypath}/${barcode}/ -type f ! -name "*.gz" -exec pigz {} \\;
    cat ${mypath}/${barcode}/@(*.fastq|*.fq).gz > ${samplename}.fastq.gz
    """
}

process REPORT {
    container 'aangeloo/nxf-tgs:latest'
    tag "$user"
    //publishDir "$params.outdir/$user", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(user), path(fastqfiles)
    
    output:
    path('*faster-report.tsv')
    
    script:
    """
    echo "file\treads\tbases\tn_bases\tmin_len\tmax_len\tN50\tGC_percent\tQ20_percent" > 01-${user}-faster-report.tsv
    # parallel -k faster -ts ::: $fastqfiles >> 00-${user}-faster-report.tsv
    faster2 -ts $fastqfiles >> 01-${user}-faster-report.tsv
    """
}

process HTMLREPORT {
    container 'aangeloo/nxf-tgs:latest'
    tag "$user"
    errorStrategy 'retry'
    maxRetries 3
    publishDir "$params.outdir/$user", mode: 'copy', pattern: '*.html' // used for sharing with the user
    //publishDir "$params.outdir", mode: 'copy', pattern: '*.html' // v5.1.14 of epi2me-labs looks for html reports recursively, so not needed
    
    input:
    tuple val(user), path(fastqpath)

    output:
    path('*.html')

    script:
    """
    # get run info from fastq header to use in faster-report
    # it is the same for all users
    FLOWCELL=\$(gzip -cd ${fastqpath[0]} | head -n 1 | grep -oE "flow_cell_id=.*" | cut -d" " -f1 | cut -d= -f2)
    RUNDATE=\$(gzip -cd ${fastqpath[0]} | head -n 1 | grep -oE "start_time=.*" | cut -d" " -f1 | cut -d= -f2 | cut -dT -f1)
    BC_MODEL=\$(gzip -cd ${fastqpath[0]} | head -n 1 | grep -oE "model_version_id=.*" | cut -d" " -f1 | cut -d= -f2 | cut -dT -f1)

    if [[ -z "\$FLOWCELL" ]]; then
        FLOWCELL="NA"
    fi

    if [[ -z "\$RUNDATE" ]]; then
        RUNDATE="NA"
    fi
    
    if [[ -z "\$BC_MODEL" ]]; then
        BC_MODEL="NA"
    fi

    faster-report.R -p . \
        --outfile 01-${user}-faster-report \
        --user ${user} \
        --rundate \$RUNDATE \
        --flowcell \$FLOWCELL \
        --basecall \$BC_MODEL
    """
}

// no need to run this in docker as it is already dockerized
process ASSEMBLY {
    tag "$user"
    errorStrategy 'retry'
    publishDir (
        "$params.outdir/$user", 
        mode: "copy", 
        pattern: "02-assembly/**{txt,fasta,fai,fastq,gbk,bam,bai,json,annotations.bed}" //wf html report is handled separately
    )
    publishDir ( 
        "$params.outdir/$user", 
        mode: "copy", 
        pattern: "02-assembly/*html", 
        saveAs: { fn -> "02-${user}-${file(fn).baseName}.html" } // rename wf-report to add username 
    ) 
    // [user, /path/to/samplesheet.csv, /path/to/fastq_pass, version]
    input:
    tuple val(user), path(samplesheet), path(fastq_pass), val(ver)
    
    output:
    //path "output/*report.html"
    path "**"
    // this is not output by wf-amplicon and bacterial genome, so no mapping and IGV report there
    tuple val(user), path("02-assembly/*.final.fasta"), path("02-assembly/*.annotations2.bed"), optional: true, emit: fasta_ch
    tuple val(user), path("02-assembly/sample_status.txt"), optional: true, emit: sample_status_ch
    tuple val(user), path("02-assembly/*.assembly_stats.tsv"), optional: true, emit: assembly_stats_ch
    
    script:
    def assembly_args = params.assembly_args ?: ''
    """
    nextflow run epi2me-labs/${params.pipeline} \
        --fastq $fastq_pass \
        --sample_sheet $samplesheet \
        --out_dir '02-assembly' \
        ${assembly_args} \
        -r $ver
    # fix annotations.bed
    # this has to be moved out in IGV process to be able to run in docker because of the R libraries. Or use base R!
    if [ ${params.pipeline} = 'wf-clone-validation' ]; then
        # if feature_table is empty do not run make_bed.R - implemented in make_bed.R
        # feature_counts=\$(wc -l < 02-assembly/feature_table.txt)
        cd 02-assembly && make_bed.R feature_table.txt && cd ..
    fi

    # put an empty assembly_stats.tsv in case no sample for this user works, to avoid having null in assembly_stats_ch
    if [ ! -f 02-assembly/*.assembly_stats.tsv ]; then
        touch 02-assembly/empty.assembly_stats.tsv
        echo -e "sample_name\tmean_quality" >> 02-assembly/empty.assembly_stats.tsv
    fi
    """
}

process SAMPLE_STATUS {
    container 'aangeloo/nxf-tgs:latest'
    tag "$user"
    errorStrategy 'ignore'

    publishDir(
        "$params.outdir/$user", 
        mode: "copy", 
        pattern: "*.csv",
        saveAs: { fn -> "00-${user}-${file(fn).baseName}.csv" } 
    )

    //
    // stage assembly_stats.tsv files but don't pass as arguments to the Rscript
    input:
    tuple val(user), path(sample_status), path(samplesheet_validated), path(assembly_stats)

    output:
    path("*.csv"), emit: merged_sample_status_ch

    script:
    // *assembly_stats.tsv is not there for failing samples, if no samples succeed for a user then assembly_stats_ch is not emitted
    // if remainder = true, then assembly_stats_ch is null
    
    """
    sample_status.R ${user} ${sample_status} ${samplesheet_validated}
    """
}

// https://www.nextflow.io/docs/latest/process.html#multiple-input-files
process SAMPLE_SUMMARY {
    container 'aangeloo/nxf-tgs:latest'
    errorStrategy 'ignore'
    publishDir "$params.outdir", mode: "copy", pattern: "*.html"

    input:
    path "sample-status*.csv"

    output:
    path("*.html")

    script:
    """
    sample_summary.R "*.csv" $workflow.runName
    """
}

process MAPPING {
    container 'aangeloo/nxf-tgs:latest'
    tag "$user - $sample"
    publishDir "$params.outdir/$user/03-mapping", mode: 'copy'

    //[user, sample, [final.fasta, annotations.bed], fastq.gz]
    input:
    tuple val(user), val(sample), path(finalfasta), path(fastq)

    output:
    path "*.{bam,bai,tsv}"
    tuple val(user), val(sample), path("*.{bam,bai,problems.tsv}"), emit: bam_ch

    script:
    """
    minimap2 -ax lr:hq ${finalfasta[0]} $fastq > mapping.sam

    samtools view -S -b -T ${finalfasta[0]} mapping.sam | \
    samtools sort -o ${sample}.bam -
    samtools index ${sample}.bam
    rm mapping.sam

    perbase base-depth --threads 4 ${sample}.bam -F 260 > ${sample}.perbase.tsv
    perpos_freq.sh ${sample}.perbase.tsv > ${sample}.problems.tsv
    """

}

process IGV_REPORTS {
    container 'aangeloo/nxf-tgs:latest'
    // https://training.nextflow.io/basic_training/debugging/#dynamic-resources-allocation
    errorStrategy { task.exitStatus == 1 ? 'retry' : 'ignore' }
    tag "$user - $sample"
    publishDir "$params.outdir/$user/04-igv-reports", mode: 'copy'
    //[user, sample, [bam, bam.bai, problems.tsv], [final.fasta, annotations2.bed]]
    input:
    tuple val(user), val(sample), path(mapping), path(fasta)

    output:
    path "*.igvreport.html"

    script:
    """
    len=\$(faster2 -l ${fasta[0]})
    header=\$(grep ">" ${fasta[0]} | cut -c 2-)
    # dynamic calculation for subsampling, subsample for > 500 alignments
    count=\$(samtools view -c ${mapping[0]})
    subsample=\$(echo \$count | awk '{if (\$1 <500) {print 1} else {print 500/\$1}}')
    
    # construct bed file
    echo -e "\$header\t0\t\$len\tsubsampled alignments (\$subsample)" > bedfile.bed
    awk -v OFS='\t' -v chr=\$header 'NR>1 {print chr, \$1-1, \$1, "HET"}' ${mapping[2]} >> bedfile.bed

    create_report \
        bedfile.bed \
        --fasta ${fasta[0]} \
        --tracks ${fasta[1]} ${mapping[0]} \
        --output ${sample}.igvreport.html \
        --flanking 200 \
        --subsample \$subsample
    """
}

workflow prep_samplesheet {
    main:
    if (params.samplesheet.endsWith(".csv")) {
        VALIDATE_SAMPLESHEET(samplesheet_ch, fastq_pass_ch) 
        .tap {validated_samplesheet_ch }
        .splitCsv(header: true)
        .filter{ it -> it.barcode =~ /^barcode*/ }
        .filter{it -> it.validate =~ /OK/ }
        .tap { internal_ch }
        .map { row -> tuple(row.sample, row.barcode, row.user) } 
        .combine(fastq_pass_ch) 
        //.view()
        .set { prepped_samplesheet_ch } 
    } else if (params.samplesheet.endsWith(".xlsx")) {
        VALIDATE_SAMPLESHEET(READEXCEL(samplesheet_ch), fastq_pass_ch) \
        .splitCsv(header: true)
        .filter{it -> it.barcode =~ /^barcode*/}
        .filter{it -> it.validate =~ /OK/ }
        .tap { internal_ch }
        .map { row -> tuple(row.sample, row.barcode, row.user) }
        .combine(fastq_pass_ch)
        .set { prepped_samplesheet_ch } 
    } else {
        exit 'Please provide either a .csv or a .xlsx samplesheet'
    }
    // generate sample sheets per user and save as files
    assembly_ch = internal_ch
        .collectFile(keepHeader: true, storeDir: "${workflow.workDir}/samplesheets"){ row ->
            sample      = row.sample
            barcode     = row.barcode
            size        = row.dna_size
            user        = row.user
            ["${user}_samplesheet.csv", "alias,barcode,approx_size\n${sample},${barcode},${size}\n"]
        }
        .map { file -> 
            def key = file.name.toString().tokenize('_').get(0)
            return tuple(key, file)
        }
    
    emit:
    prepped_samplesheet_ch
    assembly_ch
    validated_samplesheet_ch
}

workflow merge_reads {
    //prep_samplesheet()
    prep_samplesheet().prepped_samplesheet_ch \
    | MERGE_READS 
}

workflow report {
    prep_samplesheet().prepped_samplesheet_ch \
    | MERGE_READS \
    | groupTuple \
    | (REPORT & HTMLREPORT)

    emit:
    fastq_ch = MERGE_READS.out.merged_fastq_ch
    new_assembly_ch = prep_samplesheet.out.assembly_ch
    validated_samplesheet_ch2 = prep_samplesheet.out.validated_samplesheet_ch
}

//barcode,alias,approx_size are needed by epi2me/wf
workflow {
    report() 
    if (params.pipeline != 'report-only') {
        report.out.new_assembly_ch
        .combine(fastq_pass_ch)
        .combine(wf_ver.flatten().last())
        //.join(assembly_versions, by: [0,3]) \
        | ASSEMBLY

        report.out.fastq_ch
        .map{ it -> [ it[0], it.toString().split("/").last().split("\\.")[0], it[1] ] }
        .set { fastq_ch }
    
        ASSEMBLY.out.fasta_ch
        .transpose()
        .map{ it -> [ it[0], it.toString().split("/").last().split("\\.")[0], it[1..2] ] }
        .join(fastq_ch, by:[0,1])
        //.view()
        .set { mapping_ch }
    }

    if (params.pipeline == 'wf-clone-validation') {
        ASSEMBLY.out.sample_status_ch
        .combine(report.out.validated_samplesheet_ch2)
        .join(ASSEMBLY.out.assembly_stats_ch, remainder: true)
        //.view()
        | SAMPLE_STATUS
        
        SAMPLE_STATUS.out.merged_sample_status_ch
        .collect()
        | SAMPLE_SUMMARY

        MAPPING(mapping_ch)
        
        MAPPING.out.bam_ch
        .join( mapping_ch, by: [0,1] )
        .map{ it -> it[0..3] } 
        //[user, sample, [bam, bam.bai, problems.tsv], [final.fasta, annotations2.bed]]
        //.view()
        | IGV_REPORTS
    }
}
