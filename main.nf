
wfVersionMap = ['wf-clone-validation': 'v1.8.4', 'wf-amplicon': 'v1.2.2', 'wf-bacterial-genomes': 'v2.0.2']

def writePipelineLog() {
    def outdir = new File(params.outdir ?: 'output')
    outdir.mkdirs()

    def logFile = new File(outdir, 'pipeline.log')
    def cmdLine = workflow.commandLine ?: 'N/A'
    def hostname = java.net.InetAddress.localHost.hostName ?: 'unknown'
    def osName = System.getProperty('os.name') ?: 'unknown'
    def osVersion = System.getProperty('os.version') ?: 'unknown'
    def osArch = System.getProperty('os.arch') ?: 'unknown'
    def javaVersion = System.getProperty('java.version') ?: 'unknown'
    def nextflowVersion = nextflow.version ?: 'unknown'
    def dockerVersion = 'unknown'
    try {
        dockerVersion = 'docker -v'.execute().text.trim()
    } catch (Exception e) {}
    def processors = Runtime.runtime.availableProcessors()
    def pipelineVersion = wfVersionMap[params.pipeline] ?: 'N/A'
    def lines = []

    lines << "============================================================"
    lines << "NXF-TGS PIPELINE RUN ENVIRONMENT"
    lines << "============================================================"
    lines << "${new java.util.Date()}"
    lines << ""
    lines << "command_line: "
    lines << "------------------------------------------------------------"
    lines << "${cmdLine}"
    lines << ""
    lines << "machine_info:"
    lines << "------------------------------------------------------------"
    lines << "hostname           : ${hostname}"
    lines << "user               : ${System.getProperty('user.name') ?: 'unknown'}"
    lines << "available_cores    : ${processors}"
    lines << "java_version       : ${javaVersion}"
    lines << "nextflow_version   : ${nextflowVersion}"
    lines << "docker_version     : ${dockerVersion}"
    lines << ""
    lines << "os_info:"
    lines << "------------------------------------------------------------"
    lines << "name               : ${osName}"
    lines << "version            : ${osVersion}"
    lines << "architecture       : ${osArch}"
    lines << ""
    lines << "pipeline:"
    lines << "------------------------------------------------------------"
    lines << "nxf_tgs_version    : ${workflow.manifest.version ?: 'unknown'}"
    lines << "pipeline           : ${params.pipeline}"
    lines << "pipeline_version   : ${pipelineVersion}"
    lines << ""
    lines << "parameters:"
    lines << "------------------------------------------------------------"

    params.keySet().sort().each { key ->
        def value = params[key] != null ? params[key].toString() : 'null'
        lines << "${key.padRight(18)} : ${value}"
    }

    lines << "------------------------------------------------------------"

    logFile.text = lines.join('\n') + '\n'
    log.info "Wrote pipeline metadata to ${logFile}"
}

writePipelineLog()

if (params.help) {
    helpMessage()
    exit(0)
}

def helpMessage() {
    log.info """\
    ========================================================================================================================
    NXF - TGS ONT PIPELINE
    process (per user) raw fastq_pass/bam_pass folder - merge/rename, generate report, assembly (plasmid, amplicon, bacterial genome)
    ========================================================================================================================
    Usage:
    -----------------------------------
    reads            : path to raw data folder (fastq_pass or bam_pass, as output by MinKNOW), auto-detects bam or fastq input based on file extensions in the barcode directories
    samplesheet      : path to csv or excel with (at least) columns user, sample, barcode, dna_size
    pipeline         : epi2me workflow to use - can be wf-clone-validation, wf-bacterial-genomes, wf-amplicon, report-only
    assembly_args    : additional command-line arguments passed to the assembly workflow
    assembly_profile : profile to use for the assembly workflow (standard, singularity, test), default is 'standard'
    outdir           : where to save results, default is 'output'
    nxf_ver          : nextflow version to use for the internal nextflow run, default is 24.04.2
    cpus             : number of cpus to use, default is 4 (this is the executor local total cpus, also passed to assembly workflow! - minimum is 4) 
    -----------------------------------
    """
    .stripIndent(true)
}

log.info """\
    ========================================================================================================================
    NXF - TGS ONT PIPELINE
    process (per user) raw fastq_pass/bam_pass folder - merge/rename, generate reports, assembly (plasmid, amplicon, bacterial genome)
    ========================================================================================================================
    reads           : ${params.reads}
    samplesheet     : ${params.samplesheet}
    pipeline        : ${params.pipeline}
    assembly_args   : ${params.assembly_args}
    assembly_profile: ${params.assembly_profile}
    outdir          : ${params.outdir}
    nxf_ver         : ${params.nxf_ver}
    cpus            : ${params.cpus}
    -----------------------------------------------------------------------------------------------------------------------
    """
    .stripIndent(true)

reads_dir_ch = Channel.fromPath(params.reads, type: 'dir', checkIfExists: true)
samplesheet_ch = Channel.fromPath(params.samplesheet, type: 'file', checkIfExists: true)
wf_versions = Channel.from(wfVersionMap.collect { k, v -> [k, v] })
wf_ver = Channel.from(params.pipeline).join(wf_versions)

// Auto-detect bam input once at startup by checking file extensions in barcode subdirectories
isBamInput = file(params.reads).listFiles()?.any { subdir ->
    subdir.isDirectory() && subdir.listFiles()?.any { it.name.endsWith('.bam') }
} ?: false

// has to be a value channel
reference_ch = Channel.value( file( "${projectDir}/assets/mg1655.fasta" ) )

// takes in csv, checks for duplicate barcodes, unique sample names per user
// emits *-checked.csv if all ok, exits with error if not
// also add observed peak size (as seen by fasterplot) and nreads per barcode if plasmid workflow
process VALIDATE_SAMPLESHEET {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    publishDir "$params.outdir", mode: 'copy', pattern: '00-samplesheet-validated.csv'

    input: 
    path(csv)
    path(reads_dir)

    output:
    path("00-samplesheet-validated.csv")

    script:
    """
    validate_samplesheet.R $csv $reads_dir
    
    if [ ${params.pipeline} != 'wf-bacterial-genomes' ]; then
        get_maxbin.sh samplesheet-validated.csv $reads_dir
    else
        mv samplesheet-validated.csv 00-samplesheet-validated.csv
    fi

    """
}

process READEXCEL {
    container 'docker.io/aangeloo/nxf-tgs:latest'

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
    container 'docker.io/aangeloo/nxf-tgs:latest'
    tag "$user - $samplename"
    // missing barcodes are filtered out upstream in prep_samplesheet, so real errors surface here
    publishDir "$params.outdir/$user/01-fastq", mode: 'copy', pattern: '*.fastq.gz'
    publishDir "$params.outdir/$user/01-bam", mode: 'copy', pattern: '*.bam', enabled: isBamInput

    input:
    tuple val(samplename), val(barcode), val(user), path(mypath)
    
    output: 
    tuple val(user), path('*.fastq.gz'), emit: merged_fastq_ch
    tuple val(user), path('*.bam'), optional: true, emit: merged_bam_ch
    
    script:
    // Auto-detect bam or fastq based on file extensions in the barcode directory
    """
    shopt -s extglob nullglob

    bam_files=( ${mypath}/${barcode}/*.bam )
    if [ \${#bam_files[@]} -gt 0 ]; then
        # bam mode: merge bams per barcode (goes to HTMLREPORT), also derive fastq for mapping
        samtools merge -@ ${task.cpus} -o ${samplename}.bam ${mypath}/${barcode}/*.bam
        samtools fastq -T '*' ${samplename}.bam | pigz > ${samplename}.fastq.gz
    else
        #https://www.gnu.org/savannah-checkouts/gnu/bash/manual/bash.html#Pattern-Matching
        # never modify the raw data dir - ASSEMBLY reads it at the same time (in-place pigz caused duplicated reads)
        gz_files=( ${mypath}/${barcode}/@(*.fastq|*.fq).gz )
        plain_files=( ${mypath}/${barcode}/@(*.fastq|*.fq) )

        if [ \${#gz_files[@]} -eq 0 ] && [ \${#plain_files[@]} -eq 0 ]; then
            echo "ERROR: No fastq/fq or bam files found in ${mypath}/${barcode}/." >&2
            exit 1
        fi

        # concatenated gzip streams are a valid gzip file
        > ${samplename}.fastq.gz
        if [ \${#gz_files[@]} -gt 0 ]; then
            cat "\${gz_files[@]}" >> ${samplename}.fastq.gz
        fi
        if [ \${#plain_files[@]} -gt 0 ]; then
            cat "\${plain_files[@]}" | pigz -p ${task.cpus} >> ${samplename}.fastq.gz
        fi
    fi
    """
}

process REPORT {
    container 'docker.io/aangeloo/nxf-tgs:latest'
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
    //container 'docker.io/aangeloo/faster-report:latest'
    tag "$user"
    maxForks 1 // This ensures sequential execution for this process, to prevent faster-report.knit.md which is in the script dir to be mixed up
    errorStrategy 'retry'
    maxRetries 3
    publishDir "${params.outdir}/$user", 
        mode: 'copy', 
        saveAs: { filename -> file(filename).getName() }
    
    input:
    tuple val(user), path(readsfiles)

    output:
    path('output/*.html')

    script:
    """
    nextflow run angelovangel/faster-report \
    --reads . \
    --outfile 01-${user}-rawreads-report.html \
    --user ${user} \
    --simgel
    """
}

// wf-amplicon does not accept --bam, so convert bam_pass to a fastq_pass-like dir (same barcode subdirs)
process BAM_TO_FASTQ {
    container 'docker.io/aangeloo/nxf-tgs:latest'

    input:
    path(reads_dir)

    output:
    path('fastq_pass')

    script:
    """
    shopt -s nullglob
    for dir in ${reads_dir}/*/; do
        barcode=\$(basename \$dir)
        bam_files=( \$dir*.bam )
        [ \${#bam_files[@]} -eq 0 ] && continue
        mkdir -p fastq_pass/\$barcode
        samtools cat "\${bam_files[@]}" | samtools fastq -@ ${task.cpus} -T '*' - | pigz > fastq_pass/\$barcode/\${barcode}.fastq.gz
    done
    """
}

// no need to run this in docker as it is already dockerized
process ASSEMBLY {
    tag "$user - ${params.pipeline} $ver"
    errorStrategy 'retry'
    publishDir (
        "$params.outdir/$user", 
        mode: "copy", 
        pattern: "02-assembly/**{fasta,fasta.gz,fai,fastq,gbk,gbff,gff3,bam,bai,flye_stats.tsv}" //wf html report is handled separately
    )
    publishDir ( 
        "$params.outdir/$user", 
        mode: "copy", 
        pattern: "02-assembly/*html", 
        saveAs: { fn -> "03-${user}-${file(fn).baseName}.html" } // rename wf-report to add username 
    ) 
    // [user, /path/to/samplesheet.csv, /path/to/reads_dir, version]
    input:
    tuple val(user), path(samplesheet), path(reads_dir), val(ver)
    
    output:
    //path "output/*report.html"
    path "**"
    // this is not output by wf-amplicon and bacterial genome, so no mapping and IGV report there
    tuple val(user), path("02-assembly/*.final.fasta"), path("02-assembly/*.annotations2.bed"), optional: true, emit: assembly_fasta_ch
    tuple val(user), path("02-assembly/sample_status.txt"), optional: true, emit: sample_status_ch
    tuple val(user), path("02-assembly/*.assembly_stats.tsv"), optional: true, emit: assembly_stats_ch
    tuple val(user), path("02-assembly/amplicon_sample_status.txt"), optional: true, emit: amplicon_status_ch
    tuple val(user), path("02-assembly/results.json"), optional: true, emit: bacterial_results_ch
    tuple val(user), path("02-assembly/*/alignments/*.aligned.sorted.bam"), optional: true, emit: amplicon_bam_ch
    
    script:
    def assembly_args = params.assembly_args ?: ''
    def custom_configs = workflow.configFiles.findAll { !it.name.endsWith('nextflow.config') }
    def append_configs = custom_configs ? custom_configs.collect { "cat ${it} >> child.config" }.join('\n    ') : ''
    // cap memory requests of the wf-pipeline processes (low_mem profile)
    def resource_limits = params.max_memory ? "process { resourceLimits = [memory: '${params.max_memory}', cpus: ${params.cpus}] }" : ''
    // bakta occasionally hangs on single samples - give up on annotation for those instead of blocking the run
    def bakta_limit = params.pipeline == 'wf-bacterial-genomes' ? "process { withName: 'runBakta' { time = '${params.bakta_time}'; errorStrategy = 'ignore' } }" : ''
    // the epi2me wf runs in a stable dir per session and user, so that a retry of this task resumes it (-resume)
    def child_dir = "${workflow.workDir}/epi2me/${workflow.sessionId}-${user}"
    // wf-amplicon gets fastq converted by BAM_TO_FASTQ
    def reads_flag = isBamInput && params.pipeline != 'wf-amplicon' ? '--bam' : '--fastq'
    """
    # do this in this shell, or better set it up in the calling shell!
    export NXF_SINGULARITY_CACHEDIR="\$HOME/singularity-cache"
    
    echo "executor {
      name = 'local'
      cpus = ${params.cpus}
    }" > child.config
    echo "${resource_limits}" >> child.config
    echo "${bakta_limit}" >> child.config
    
    ${append_configs}

    # absolute, resolved input paths - they must be identical between attempts for -resume to work
    TASK_DIR=\$PWD
    READS=\$(readlink -f $reads_dir)
    SHEET=\$(readlink -f $samplesheet)

    mkdir -p ${child_dir} && cd ${child_dir}
    NXF_VER=${params.nxf_ver} nextflow run epi2me-labs/${params.pipeline} \
        ${reads_flag} \$READS \
        --sample_sheet \$SHEET \
        --out_dir \$TASK_DIR/02-assembly \
        ${assembly_args} \
        -c \$TASK_DIR/child.config \
        -w ${child_dir}/work \
        -resume \
        -r $ver \
        -profile ${params.assembly_profile}
    cd \$TASK_DIR
        
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

    # for wf-amplicon: derive pass/fail from all-consensus-seqs.fasta
    if [ ${params.pipeline} = 'wf-amplicon' ]; then
        amplicon_sample_status.sh $samplesheet 02-assembly
    fi
    """
}

// wf-amplicon: mean depth and primary ratio from reads aligned to the final consensus
// primary_ratio = primary / (primary + secondary + supplementary), as in the wf-amplicon de novo QC
process AMPLICON_ALN_STATS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    tag "$user"
    errorStrategy 'ignore'

    input:
    tuple val(user), path(bams)

    output:
    tuple val(user), path("amplicon.alignment_stats.tsv")

    script:
    """
    printf "sample_name\tmean_depth\tprimary_ratio\n" > amplicon.alignment_stats.tsv
    for bam in ${bams}; do
        sample=\$(basename \$bam .aligned.sorted.bam)
        primary=\$(samtools view -c -F 0x904 \$bam)
        other=\$(( \$(samtools view -c -f 0x100 \$bam) + \$(samtools view -c -f 0x800 \$bam) ))
        depth=\$(samtools depth -a \$bam | awk '{ s += \$3 } END { if (NR > 0) printf "%.0f", s / NR; else print "NA" }')
        ratio=\$(awk -v p=\$primary -v o=\$other 'BEGIN { if (p + o > 0) printf "%.3f", p / (p + o); else print "NA" }')
        printf "%s\t%s\t%s\n" \$sample \$depth \$ratio >> amplicon.alignment_stats.tsv
    done
    """
}

process SAMPLE_STATUS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    tag "$user"
    errorStrategy 'ignore'

    // publishDir(
    //     "$params.outdir/$user", 
    //     mode: "copy", 
    //     pattern: "*.csv",
    //     saveAs: { fn -> "00-${user}-${file(fn).baseName}.csv" } 
    // )

    // sample_status.R finds *.assembly_stats.tsv files in the work dir via glob - no need to pass them explicitly
    // for wf-bacterial-genomes, sample_status is results.json
    input:
    tuple val(user), path(sample_status), path(samplesheet_validated), path(assembly_stats)

    output:
    path("*.csv"), emit: merged_sample_status_ch
    tuple val(user), path("*.csv"), emit: user_sample_status_ch

    script:
    def status_script = params.pipeline == 'wf-bacterial-genomes' ? 'bacterial_sample_status.R' : 'sample_status.R'
    """
    ${status_script} ${user} ${sample_status} ${samplesheet_validated}
    """
}

// https://www.nextflow.io/docs/latest/process.html#multiple-input-files
process SAMPLE_SUMMARY {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    errorStrategy 'ignore'
    publishDir "$params.outdir", mode: "copy", pattern: "*.html"

    input:
    path "sample-status*.csv"

    output:
    path("*.html")

    script:
    """
    gitcommit=\$(cat $workflow.projectDir/.git/logs/HEAD  | cut -f2 -d " " | tail -1 | cut -c-7)
    sample_summary.R "*.csv" $workflow.runName \$gitcommit ${params.pipeline} ${wfVersionMap[params.pipeline] ?: 'N/A'}
    """
}

process MAPPING {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    tag "$user - $sample"
    publishDir "$params.outdir/$user/03-mapping", mode: 'copy', pattern: "*.{bam,bai}"

    //[user, sample, [final.fasta, annotations.bed], fastq.gz]
    input:
    tuple val(user), val(sample), path(finalfasta), path(fastq)
    path mg1655

    output:
    path "*.{bam,bai,tsv}"
    tuple val(user), val(sample), path("*.{bam,bai,problems.tsv}"), emit: mapping_bam_ch
    tuple val(user), path("*mapping-counts.txt"), emit: mapping_counts_ch//, optional: true

    script:
    //def mg1655 = file("${projectDir}/assets/mg1655.fasta")
    """
    minimap2 -ax lr:hq ${finalfasta[0]} $fastq > mapping.sam

    samtools view -S -b -T ${finalfasta[0]} mapping.sam | \
    samtools sort -o ${sample}.bam -
    samtools index ${sample}.bam
    rm mapping.sam

    perbase base-depth --threads 4 ${sample}.bam -F 260 > ${sample}.perbase.tsv
    perpos_freq.sh ${sample}.perbase.tsv > ${sample}.problems.tsv

    # get mapping statistics for mapping to assembly and to e. coli genome
    # only for plasmids
    
    if [ ${params.pipeline} == 'wf-clone-validation' ]; then
        # replace fasta header to "sample"
        seqkit replace -p ".*" -r "assembly" ${finalfasta[0]} > target.fasta
        cat ${mg1655} >> target.fasta
        minimap2 -x lr:hq --secondary=no target.fasta $fastq > mapping_counts.paf
        allreads=\$(faster2 -l $fastq | wc -l | tr -d " ")
        
        # filter mapq > 60 and alen/qlen > 0.4 
        awk '{if (\$12 >= 60 && \$11/\$2 >= 0.4 ) print \$1,\$6}' mapping_counts.paf | sort | uniq | cut -d" " -f2 | sort | uniq -c > temp.txt
        echo "\$allreads allreads" >> temp.txt
        awk '{print \$2}' temp.txt | paste -sd ' ' - > $sample-mapping-counts.txt
        awk '{print \$1}' temp.txt | paste -sd ' ' - >> $sample-mapping-counts.txt  
    fi
    
    """

}

process MAPPING_COUNTS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    tag "$user"
    // publishDir(
    //     "$params.outdir/$user", 
    //     mode: "copy", 
    //     pattern: "*.csv",
    //     saveAs: { fn -> "00-${user}-${file(fn).baseName}.csv" } 
    // )

    input:
    tuple val(user), path(mapping_counts)

    output:
    path('*.csv'), emit: merged_mapping_counts_ch
    tuple val(user), path('*.csv'), emit: user_mapping_counts_ch

    script:
    """
    mapping_counts.R "*mapping-counts.txt" $user
    """
}

process MAPPING_SUMMARY {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    publishDir "$params.outdir", mode: 'copy'

    input:
    path('*mapping-summary.csv')

    output:
    path("*.html")

    script:
    """
    mapping_summary.R "*.csv" $workflow.runName
    """ 

}

process USER_REPORT {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    tag "$user"
    publishDir "$params.outdir/$user", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(user), path(sample_status), path(mapping_summary)

    output:
    path("*.html")

    script:
    def pipeline_label = wfVersionMap[params.pipeline] ?: 'N/A'
    def ms_arg = mapping_summary ? mapping_summary : 'dummy_mapping.csv'
    """
    user_report.R "$user" "$sample_status" "$ms_arg" "${params.pipeline} ${pipeline_label}" "${workflow.manifest.version ?: 'unknown'}"
    """
}

process IGV_REPORTS {
    //container 'docker.io/aangeloo/nxf-tgs:latest'
    container 'docker.io/aangeloo/igv-reports-image:latest'
    errorStrategy 'ignore'
    tag "$user - $sample"
    publishDir "$params.outdir/$user/04-igv-reports", mode: 'copy'
    //[user, sample, [bam, bam.bai, problems.tsv], [final.fasta, annotations2.bed]]
    input:
    tuple val(user), val(sample), path(mapping), path(fasta)

    output:
    path "*.igvreport.html"

    script:
    """
    # len=\$(faster2 -l ${fasta[0]})
    samtools faidx ${fasta[0]}
    len=\$(awk '{print \$2}' ${fasta[0]}.fai | head -n 1)

    header=\$(grep ">" ${fasta[0]} | cut -c 2-)
    # dynamic calculation for subsampling, subsample for > 500 alignments
    count=\$(samtools view -c ${mapping[0]})
    percent_aln=\$(samtools flagstat ${mapping[0]} | grep 'primary mapped' | cut -d"(" -f 2 | cut -d" " -f1)
    all_reads=\$(samtools flagstat ${mapping[0]} | grep 'primary\$' | cut -d" " -f1)

    subsample=\$(echo \$count | awk '{if (\$1 <500) {print 1} else {print 500/\$1}}')
    
    # construct bed file
    echo -e "\$header\t0\t\$len\tPrimary alignments: \$percent_aln of \$all_reads reads" > bedfile.bed
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
        VALIDATE_SAMPLESHEET(samplesheet_ch, reads_dir_ch) 
        .tap {validated_samplesheet_ch }
        .splitCsv(header: true)
        .filter{ it -> it.barcode =~ /^barcode*/ }
        .tap { all_rows_ch }  // capture before OK filter
        .filter{it -> it.validate =~ /OK/ }
        .tap { validated_rows_ch }
        .map { row -> tuple(row.sample, row.barcode, row.user) } 
        .combine(reads_dir_ch)
        .set { samples_ch }

        all_rows_ch
        .filter { it.validate != 'OK' }
        .view { row -> "WARN: Skipping sample '${row.sample}' (${row.barcode}, user: ${row.user}, reason: ${row.validate})" }
    } else if (params.samplesheet.endsWith(".xlsx")) {
        VALIDATE_SAMPLESHEET(READEXCEL(samplesheet_ch), reads_dir_ch)
        .tap {validated_samplesheet_ch }
        .splitCsv(header: true)
        .filter{it -> it.barcode =~ /^barcode*/}
        .tap { all_rows_ch }  // capture before OK filter
        .filter{it -> it.validate =~ /OK/ }
        .tap { validated_rows_ch }
        .map { row -> tuple(row.sample, row.barcode, row.user) }
        .combine(reads_dir_ch)
        .set { samples_ch }

        all_rows_ch
        .filter { it.validate != 'OK' }
        .view { row -> "WARN: Skipping sample '${row.sample}' (${row.barcode}, user: ${row.user}, reason: ${row.validate})" }
    } else {
        exit 'Please provide either a .csv or a .xlsx samplesheet'
    }
    // generate sample sheets per user and save as files
    user_samplesheet_ch = validated_rows_ch
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
    samples_ch
    user_samplesheet_ch
    validated_samplesheet_ch
}

workflow merge_reads {
    //prep_samplesheet()
    samples_ch = prep_samplesheet().samples_ch
    samples_ch | MERGE_READS
}

// check this for potential mixup of users and samples
workflow report {
    samples_ch = prep_samplesheet().samples_ch
    // ASSEMBLY (epi2me wf) ingests fastq or bam directly, so pass the original data (except wf-amplicon, which has no --bam)
    assembly_reads_ch = isBamInput && params.pipeline == 'wf-amplicon' ? BAM_TO_FASTQ(reads_dir_ch) : reads_dir_ch
    MERGE_READS(samples_ch)
    // Use top-level isBamInput to pick the right channel for HTMLREPORT
    htmlreport_ch = isBamInput ? MERGE_READS.out.merged_bam_ch : MERGE_READS.out.merged_fastq_ch
    htmlreport_ch \
    | groupTuple(by: 0) \
    | HTMLREPORT
    //| view()
    //| (REPORT & HTMLREPORT)

    emit:
    merged_fastq_ch = MERGE_READS.out.merged_fastq_ch
    user_samplesheet_ch = prep_samplesheet.out.user_samplesheet_ch
    validated_samplesheet_ch = prep_samplesheet.out.validated_samplesheet_ch
    assembly_reads_ch = assembly_reads_ch
}

//barcode,alias,approx_size are needed by epi2me/wf
workflow {
    report() 
    if (params.pipeline != 'report-only') {
        report.out.user_samplesheet_ch
        .combine(report.out.assembly_reads_ch)
        .combine(wf_ver.flatten().last())
        //.join(assembly_versions, by: [0,3]) \
        | ASSEMBLY

        report.out.merged_fastq_ch
        .map{ it -> [ it[0], it.toString().split("/").last().split("\\.")[0], it[1] ] }
        .set { sample_fastq_ch }
    
        ASSEMBLY.out.assembly_fasta_ch
        .transpose()
        .map{ it -> [ it[0], it.toString().split("/").last().split("\\.")[0], it[1..2] ] }
        .join(sample_fastq_ch, by:[0,1])
        //.view()
        .set { mapping_ch }
    }

    if (params.pipeline == 'wf-clone-validation') {
        ASSEMBLY.out.sample_status_ch
        .combine(report.out.validated_samplesheet_ch)
        .join(ASSEMBLY.out.assembly_stats_ch, remainder: true)
        .map { user, sample_status, samplesheet, assembly_stats ->
            return [user, sample_status, samplesheet, assembly_stats ?: []]
        }
        | SAMPLE_STATUS
        
        SAMPLE_STATUS.out.merged_sample_status_ch
        .collect()
        | SAMPLE_SUMMARY
        
        MAPPING(mapping_ch, reference_ch)
        
        MAPPING.out.mapping_bam_ch
        .join( mapping_ch, by: [0,1] )
        .map{ it -> it[0..3] } 
        //[user, sample, [bam, bam.bai, problems.tsv], [final.fasta, annotations2.bed]]
        //.view()
        | IGV_REPORTS
        
        MAPPING.out.mapping_counts_ch
        .groupTuple() // group by user
        //.view()
        | MAPPING_COUNTS
        
        MAPPING_COUNTS.out.merged_mapping_counts_ch
        .collect()
        //.view()
        | MAPPING_SUMMARY

        SAMPLE_STATUS.out.user_sample_status_ch
        .join(MAPPING_COUNTS.out.user_mapping_counts_ch, remainder: true)
        .filter { user, sample_status, mapping_summary ->
            // Skip entirely if sample_status is missing (SAMPLE_STATUS was ignored/failed)
            sample_status != null
        }
        .map { user, sample_status, mapping_summary -> 
            [user, sample_status, mapping_summary ?: []] 
        }
        | USER_REPORT
 
    } else if (params.pipeline == 'wf-amplicon') {
        ASSEMBLY.out.amplicon_status_ch
        .combine(report.out.validated_samplesheet_ch)
        .join(ASSEMBLY.out.assembly_stats_ch, remainder: true)
        .join(AMPLICON_ALN_STATS(ASSEMBLY.out.amplicon_bam_ch), remainder: true)
        // sample_status.R picks up *.assembly_stats.tsv and *.alignment_stats.tsv from the work dir
        .map { user, amplicon_status, samplesheet, assembly_stats, aln_stats ->
            return [user, amplicon_status, samplesheet, [assembly_stats, aln_stats].flatten().findAll { it }]
        }
        | SAMPLE_STATUS
 
        SAMPLE_STATUS.out.merged_sample_status_ch
        .collect()
        | SAMPLE_SUMMARY

        // no mapping for wf-amplicon and wf-bacterial-genomes
        SAMPLE_STATUS.out.user_sample_status_ch
        .map { user, sample_status -> [user, sample_status, []] }
        | USER_REPORT
    } else if (params.pipeline == 'wf-bacterial-genomes') {
        ASSEMBLY.out.bacterial_results_ch
        .combine(report.out.validated_samplesheet_ch)
        .map { user, results_json, samplesheet ->
            return [user, results_json, samplesheet, []]
        }
        | SAMPLE_STATUS

        SAMPLE_STATUS.out.merged_sample_status_ch
        .collect()
        | SAMPLE_SUMMARY

        // no mapping for wf-amplicon and wf-bacterial-genomes
        SAMPLE_STATUS.out.user_sample_status_ch
        .map { user, sample_status -> [user, sample_status, []] }
        | USER_REPORT
    }
}
