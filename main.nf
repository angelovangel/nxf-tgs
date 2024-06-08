params.fastq = "fastq_pass"
params.samplesheet = "samplesheet.csv" //user, sample, barcode NEEDED
params.outdir = "${workflow.launchDir}/output"
params.pipeline = "wf-clone-validation" //can be wf-clone-validation, wf-bacterial-genomes, wf-amplicon
params.assembly_args = null
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
    samplesheet : path to csv or excel with (at least) columns sample, barcode, user
    assembly    : epi2me workflow to use - can be wf-clone-validation, wf-bacterial-genomes, wf-amplicon
    outdir      : where to save results, default is output
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
samplesheet_ch = Channel.fromPath(params.samplesheet)
wf_versions = Channel.from(['wf-clone-validation', 'v1.2.0'], ['wf-amplicon', 'v1.1.0'], ['wf-bacterial-genomes', 'v1.2.0'])
wf_ver = Channel.from(params.pipeline).join(wf_versions)

// takes in csv, checks for duplicate barcodes, unique sample names per user
// emits *-checked.csv if all ok, exits with error if not
process VALIDATE_SAMPLESHEET {
    container 'aangeloo/nxf-tgs:latest'
    publishDir "$params.outdir", mode: 'copy', pattern: 'validated-samplesheet.csv'

    input: 
    path(csv)

    output:
    path("validated-samplesheet.csv")

    script:
    """
    validate_samplesheet.R $csv
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
    path('*faster-report.tsv')
    
    script:
    """
    echo "file\treads\tbases\tn_bases\tmin_len\tmax_len\tN50\tGC_percent\tQ20_percent" > faster-report-${user}.tsv
    #parallel -k faster -ts ::: $fastqfiles >> faster-report-${user}.tsv
    faster2 -ts $fastqfiles >> ${user}-faster-report.tsv
    """
}

process HTMLREPORT {
    container 'aangeloo/nxf-tgs:latest'
    tag "$user"
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
    if [[ -z "\$FLOWCELL" ]]; then
        FLOWCELL="NA"
    fi

    if [[ -z "\$RUNDATE" ]]; then
        RUNDATE="NA"
    fi
    
    faster-report.R -p . --outfile ${user}-faster-report --user ${user} --rundate \$RUNDATE --flowcell \$FLOWCELL
    """
}

// no need to run this in docker as it is already dockerized
process ASSEMBLY {
    tag "$user"
    publishDir (
        "$params.outdir/$user", 
        mode: "copy", 
        pattern: "02-assembly/**{txt,fasta,fastq,gbk,bed,json,bam,bai}" //wf html report is handled separately
        )
    publishDir ( 
        "$params.outdir/$user", 
        mode: "copy", 
        pattern: "02-assembly/*html", 
        saveAs: { fn -> "${user}-${file(fn).baseName}.html" } // rename wf-report to add username 
        ) 

    input:
    tuple val(user), path(samplesheet), path(fastq_pass), val(ver) // input is [user, /path/to/samplesheet.csv, /path/to/fastq_pass]
    
    output:
    //path "output/*report.html"
    path "**"
    tuple val(user), path("02-assembly/*.final.fasta"), emit: fasta_ch

    script:
    def assembly_args = params.assembly_args ?: ''
    """
    nextflow run epi2me-labs/${params.pipeline} \
        --fastq $fastq_pass \
        --sample_sheet $samplesheet \
        --out_dir '02-assembly' \
        ${assembly_args} \
        -r $ver
    """
}

process MAPPING {
    container 'aangeloo/nxf-tgs:latest'
    tag "$user - $sample"
    publishDir "$params.outdir/$user/03-mapping", mode: 'copy'

    input:
    tuple val(user), val(sample), path(finalfasta), path(fastq)

    output:
    path "*.{bam,bai,tsv}"

    script:
    """
    minimap2 -ax lr:hq $finalfasta $fastq > mapping.sam

    samtools view -S -b -T $finalfasta mapping.sam | \
    samtools sort -o ${sample}.bam -
    samtools index ${sample}.bam
    rm mapping.sam

    perbase base-depth ${sample}.bam -F 260 > ${sample}.perbase.tsv
    perpos_freq.sh ${sample}.perbase.tsv > ${sample}.problems.tsv
    """

}

workflow prep_samplesheet {
    main:
    if (params.samplesheet.endsWith(".csv")) {
        VALIDATE_SAMPLESHEET(samplesheet_ch) 
        .splitCsv(header: true)
        .filter{it -> it.barcode =~ /^barcode*/}
        .tap { internal_ch }
        .map { row -> tuple(row.sample, row.barcode, row.user) } 
        .combine(fastq_pass_ch) 
        //| view()
        .set { prepped_samplesheet_ch } 
    } else if (params.samplesheet.endsWith(".xlsx")) {
        VALIDATE_SAMPLESHEET(READEXCEL(samplesheet_ch)) \
        .splitCsv(header: true)
        .filter{it -> it.barcode =~ /^barcode*/}
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
}

//barcode,alias,approx_size are needed by epi2me/wf
workflow {
    report()                       //this is first time prep_samplesheet is called

    prep_samplesheet().assembly_ch //this is second time prep_samplesheet is called
    .combine(fastq_pass_ch)
    .combine(wf_ver.flatten().last())
    //.join(assembly_versions, by: [0,3]) \
    | ASSEMBLY

    report.out.fastq_ch
    .map{ it -> [ it[0], it.toString().split("/").last().split("\\.")[0], it[1] ] }
    .set { fastq_ch }
    
    ASSEMBLY.out.fasta_ch
    .transpose()
    .map{ it -> [ it[0], it.toString().split("/").last().split("\\.")[0], it[1] ] }
    .join(fastq_ch, by:[0,1])
    .set { mapping_ch }

    MAPPING(mapping_ch)

}
