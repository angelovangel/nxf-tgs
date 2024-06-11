#!/usr/bin/env Rscript
#============
#
# this just renders the faster-report.Rmd file,
# note that this uses system pandoc for rendering, not the Rstudio one
#  # nolint: trailing_whitespace_linter.
#
#============
require(funr)
require(optparse)
require(rmarkdown)
#require(renv)

calldir <- getwd()
scriptdir  <-  dirname(funr::sys.script())
setwd(scriptdir)
#renv::load()

option_list <- list(
  make_option(c('--path', '-p'), help = 'path to folder with fastq files [%default]', type = 'character', default = NULL),
  make_option(c('--regex', '-r'), help = 'regex pattern to match fastq files [%default]', type = 'character', default = 'fast(q|q.gz)$'),
  make_option(c('--type', '-t'), help = "seq platform used, can be one of 'illumina', 'ont' or 'pacbio' [%default]", default = 'ont'),
  make_option(c('--rundate', '-d'), help = 'Run date', type = 'character', default = NULL),
  make_option(c('--flowcell', '-f'), help = 'Flow cell ID', type = 'character', default = NULL),
  make_option(c('--basecall', '-b'), help = 'Basecaller model', type = 'character', default = NULL),
  make_option(c('--user', '-u'), help = 'User', type = 'character', default = NULL),
  make_option(c('--save_raw', '-s'), help = 'save raw csv data used for plotting [%default]', type = 'logical', default = FALSE),
  make_option(c('--subsample', '-x'), help = 'subsample reads for kmers calculation [%default]', type = 'double', default = 1.0),
  make_option(c('--outfile','-o'), help = 'name of output report file [%default]', type = 'character', default = 'faster-report.html')
  )

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

if (is.null(opts$path)){
  print_help(opt_parser)
  stop("At least a path to a folder with fastq files is required (use option '-p path/to/folder')", call.=FALSE)
}

# complicated case to parse correct fastq path when calling and script directories are not the same
# check if abs or relative path was provided
if (R.utils::isAbsolutePath(opts$path)) {
  fastqpath <- opts$path
} else {
  fastqpath <- normalizePath(file.path(calldir, opts$path))
}

print(paste0("call dir: ", calldir))
print(paste0("fastq path: ", fastqpath))

# change to match parameter used in Rmd
if (opts$type == 'illumina') {
  opts$type <- 'Illumina'
  opts$flowcell <- 'NA'
  opts$basecall <- 'NA'
} else if (opts$type == 'ont') {
  opts$type <- 'Nanopore'
} else if (opts$type == 'pacbio') {
  opts$type <- 'PacBio'
}

# render the rmarkdown, using fastq-report.Rmd as template
rmarkdown::render(input = "faster-report.Rmd",
                  output_file = opts$outfile,
                  output_dir = calldir, # important when knitting in docker
                  knit_root_dir = scriptdir, # important when knitting in docker
                  #envir = new.env(),
                  params = list(
                    fastq_dir = fastqpath,
                    fastq_pattern = opts$regex,
                    sequencer = opts$type,
                    rundate = opts$rundate,
                    flowcell = opts$flowcell,
                    basecall = opts$basecall,
                    user = opts$user,
                    rawdata = opts$save_raw,
                    subsample = opts$subsample
                  )
)
