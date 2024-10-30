#!/usr/bin/env Rscript

# reads csv/tsv and checks if valid
# checks also existence of fastq_pass/barcodeXX and presence of reads in fastq_pass/barcodeXX
# writes back *-checked.csv if ok

# arg[1] is csv, arg[2] is path to fastq_pass

arg <- commandArgs(trailingOnly = T)

library(readr)
library(stringr)
library(dplyr)
library(fs)

# valid barcode names
bc_pattern <- '^barcode[0-9]+$' 

# to print dataframe to stdout pretty
print_and_capture <- function(x)
{
  paste(capture.output(print(x)), collapse = "\n")
}

df <- readr::read_delim(arg[1], col_names = T, trim_ws = T) %>% 
  mutate(validate = 'OK')


# CHECKS #####################################################
# colnames contain 'user', 'sample', 'barcode'
if (!all( c('user', 'sample', 'barcode') %in% colnames(df) )) {
  stop(
    paste0(
      '\n--------------------------------------\n',
      '\nSamplesheet must contain columns user,sample,barcode\n',
      '\nThe provided samplesheet has columns:\n',
      str_flatten(colnames(df), collapse = ", "),
      '\n--------------------------------------\n'
    )
  )
}

# remove rows where sample, user or barcode is NA
df <- df[complete.cases(df[ ,c('sample', 'user', 'barcode')]), ]

# remove white space first
df$sample <- str_replace_all(df$sample, " ", "")
df$user <- str_replace_all(df$user, " ", "")

# barcode unique
# get indices of duplicates:) stupid R
dups_vector <- duplicated(df$barcode) | duplicated(df$barcode, fromLast = T)

if (any(dups_vector)) {
  stop(
    paste0(
      '\n--------------------------------------\n',
      '\nBarcodes must be unique!\n',
      print_and_capture(df[which(dups_vector), ]),
      '\n--------------------------------------\n'
    )
  )
}

### valid barcode names

# max len of sample name is 24 - perbase weird crash 
sl_vector <- str_length(df$sample)
if (!all(sl_vector <= 24)) {
  stop(paste0(
   '\nSamples with too long names:',
   '\n--------------------------------------\n',
   print_and_capture( df[sl_vector > 24, ] ),
   '\n--------------------------------------\n',
   'Sample names have to be < 25 characters long',
   '\n--------------------------------------\n'
  ))
}

# special characters in sample names
sn_vector <- str_detect(df$sample, '^[a-zA-Z0-9\\_\\-]+$')
if (!all(sn_vector)) {
  stop(
    paste0(
      '\nSamples with special characters:',
      '\n--------------------------------------\n',
      print_and_capture( df[!sn_vector, ] ),
      '\n--------------------------------------\n'
    )
  )
}

# duplicate sample names per user
grouped_df <- df %>% group_by(user, sample) %>% summarise(n_samples = n())
snames_vector <- grouped_df$n_samples == 1

if (any(!snames_vector)) {
  stop(
    paste0(
      '\nDuplicate sample names: ',
      '\n--------------------------------------\n',
      print_and_capture( grouped_df[!snames_vector, ] ),
      '\n--------------------------------------\n'
    )
  )
}

# only numerics in sample names
num_vector <- str_detect(df$sample, '^[0-9]+$')
if (any(num_vector)) {
  stop(
    paste0(
      '\nNumeric sample names: ',
      '\n--------------------------------------\n',
       print_and_capture( df[num_vector, ] ),
      '\n--------------------------------------\n'
    )
  )
}


# check for fastq files in barcode directory
paths <- fs::path(arg[2], df$barcode)
m <- mapply(list.files, paths)
l <- mapply(length, m) > 0
if (!all(l)) {
  warning( names(which(!l)), ' has no files')
  #df <- df[l, ]
  df$validate[!l] <- 'no fastq files'
}

# check existance of barcode dir
paths <- fs::path(arg[2], df$barcode)

bc_exists <- fs::dir_exists(paths)
if (!all(bc_exists)) {
  warning( names(which(!bc_exists)), ' does not exist', call. = F)
  #df <- df[bc_exists, ]
  df$validate[!bc_exists] <- 'bc does not exist'
}

write_csv(df, file = 'samplesheet-validated.csv')

