#!/usr/bin/env Rscript

# reads csv/tsv and checks if valid
# writes back *-checked.csv if ok

arg <- commandArgs(trailingOnly = T)

library(readr)
library(stringr)
library(dplyr)

# valid barcode names
bc_pattern <- '^barcode[0-9]+$' 

# to print dataframe to stdout pretty
print_and_capture <- function(x)
{
  paste(capture.output(print(x)), collapse = "\n")
}

df <- readr::read_delim(arg[1], col_names = T, trim_ws = T)


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

# valid barcode names


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


write_csv(df, file = 'validated-samplesheet.csv')

