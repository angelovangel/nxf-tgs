params.fastq = "fastq_pass"
params.samplesheet = "samplesheet.csv" //user, sample, barcode NEEDED
params.outdir = "${workflow.launchDir}/results-tgs"
params.pipeline = "wf-clone-validation" //can be wf-clone-validation, wf-bacterial-genomes, wf-amplicon
params.help = ""

if (params.help) {
    helpMessage()
    exit(0)
}

def helpMessage() {
    log.info """\
    ===================================
    NXF - TGS ONT PIPELINE
    process (per user) raw fastq_pass folder - merge/rename, generate report, assembly (plasmid, amplicon, bacterial genome)
    ===================================
    Usage:
    -----------------------------------
    fastq       : path to raw fastq_pass data
    samplesheet : path to csv with (at least) columns sample, barcode, user
    assembly    : epi2me workflow to use - can be wf-clone-validation, wf-bacterial-genomes, wf-amplicon
    outdir      : where to save results, default is results-tgs
    """
    .stripIndent(true)
}

log.info """\
    ===================================
    NXF - TGS ONT PIPELINE
    process (per user) raw fastq_pass folder - merge/rename, generate report, assembly (plasmid, amplicon, bacterial genome)
    ===================================
    fastq       : ${params.fastq}
    samplesheet : ${params.samplesheet}
    assembly    : ${params.pipeline}
    outdir      : ${params.outdir}
    """
    .stripIndent(true)

fastq_pass_ch = Channel.fromPath(params.fastq, type: 'dir', checkIfExists: true)
samplesheet_ch = Channel.fromPath(params.samplesheet).splitCsv(header: true)
wf_samplesheet_ch = samplesheet_ch
    .collectFile(keepHeader: true, storeDir: "${workflow.workDir}/samplesheets"){ row ->
            sample      = row.sample
            barcode     = row.barcode
            size        = row.dna_size
            user        = row.user
            ["${user}_samplesheet.csv", "alias,barcode,approx_size\n${sample},${barcode},${size}\n"]
    }

process MERGE_READS {
    container 'aangeloo/nxf-tgs:latest'
    tag "$user - $samplename"
    errorStrategy 'ignore' //because some barcodes defined in the samplesheet might be missing in the data
    publishDir "$params.outdir/$user/01-fastq", mode: 'copy', pattern: '*.fastq.gz'

    input:
    tuple val(samplename), val(barcode), val(user), path(mypath)
    
    output: 
    tuple val(user), path('*.fastq.gz')
    
    script:
    """
    cat ${mypath}/${barcode}/*.fastq.gz > ${samplename}.fastq.gz
    """
}

process REPORT {
    container 'aangeloo/nxf-tgs:latest'
    tag "$user"
    publishDir "$params.outdir/$user", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(user), path(fastqfiles)
    
    output:
    path('faster-report*')
    
    script:
    """
    echo "file\treads\tbases\tn_bases\tmin_len\tmax_len\tN50\tGC_percent\tQ20_percent" > faster-report-${user}.tsv
    #parallel -k faster -ts ::: $fastqfiles >> faster-report-${user}.tsv
    faster2 -ts $fastqfiles >> faster-report-${user}.tsv
    """
}

process HTMLREPORT {
    container 'aangeloo/nxf-tgs:latest'
    tag "$user"
    publishDir "$params.outdir/$user", mode: 'copy', pattern: '*.html'
    
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
    if [[ -z "\$FLOWCELL" ]]; then
        FLOWCELL="NA"
    fi

    if [[ -z "\$RUNDATE" ]]; then
        RUNDATE="NA"
    fi
    
    faster-report.R -p . --outfile faster-report-${user} --user ${user} --rundate \$RUNDATE --flowcell \$FLOWCELL
    """
}

// no need to run this in docker as it is already dockerized
process ASSEMBLY {
    tag "$user"
    publishDir "$params.outdir/$user", mode: "copy", pattern: "02-assembly/*{html,txt,fasta,fastq,gbk,bed,json}"
    //publishDir "$params.outdir/$user", mode: "copy", pattern: "02-assembly/*html", saveAs: { filename -> filename.getName }

    input:
    tuple val(user), path(samplesheet), path(fastq_pass) // input is [user, /path/to/samplesheet.csv, /path/to/fastq_pass]
    
    output:
    //path "output/*report.html"
    path "**"

    script:
    """
    nextflow run epi2me-labs/${params.pipeline} --fastq $fastq_pass --sample_sheet $samplesheet --out_dir '02-assembly'
    """
}

workflow merge_reads {
    samplesheet_ch
    .map { row -> tuple(row.sample, row.barcode, row.user) }
    .combine(fastq_pass_ch)
    //.view()
    | MERGE_READS
}

workflow report {
    samplesheet_ch
    .map { row -> tuple(row.sample, row.barcode, row.user) }
    .combine(fastq_pass_ch)
    //.view()
    | MERGE_READS  | groupTuple | (REPORT & HTMLREPORT)
}



//barcode,alias,approx_size are needed by epi2me/wf
workflow {
    report()
    // get user, wf-samplesheet and combine with fastq_pass
    wf_samplesheet_ch
        .map { file -> 
            def key = file.name.toString().tokenize('_').get(0)
            return tuple(key, file)
        }
        .combine(fastq_pass_ch) | ASSEMBLY
}
