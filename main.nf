params.fastq = "fastq_pass"
params.samplesheet = "samplesheet.csv" //user, sample, barcode NEEDED
params.outdir = "${workflow.launchDir}/results-tgs"
params.pipeline = "wf-clone-validation"


log.info """\
    ===================================
    NXF - TGS ONT PIPELINE
    process raw fastq_pass folder - merge/rename, generate report, assembly (plasmid, amplicon, bacterial genome)
    ===================================
    fastq_pass :  ${params.fastq}
    samplesheet:  ${params.samplesheet}
    outdir     :  ${params.outdir}
    container  :  ${workflow.container}

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
    tag "$user - $samplename"
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

process ASSEMBLY {
    tag "$user"
     publishDir "$params.outdir/$user", mode: 'copy', pattern: 'output/*'

    input:
    tuple val(user), path(samplesheet), path(fastq_pass) // input is [user, /path/to/samplesheet.csv, /path/to/fastq_pass]
    
    output:
    path 'output/*'

    script:
    """
    nextflow run epi2me-labs/${params.pipeline} --fastq $fastq_pass --sample_sheet $samplesheet
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
    | MERGE_READS  | groupTuple | REPORT
}

// workflow {
//     //report()

// }

//barcode,alias,approx_size are needed by epi2me/wf
workflow {
    report()
    // get user, wf-samplesheet and combine with fastq_pass
    wf_samplesheet_ch
        .map { file -> 
            def key = file.name.toString().tokenize('_').get(0)
            return tuple(key, file)
        }
        .combine(fastq_pass_ch) | view //ASSEMBLY
}
