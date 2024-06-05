#!/usr/bin/env Rscript

# reads csv/tsv and checks if valid
# writes back *-checked.csv if ok

arg <- commandArgs(trailingOnly = T)

library(tools)
library(readr)
library(stringr)
library(dplyr)

# to print dataframe to stdout pretty
print_and_capture <- function(x)
{
  paste(capture.output(print(x)), collapse = "\n")
}

bname <- tools::file_path_sans_ext(basename(arg[1]))
df <- readr::read_delim(arg[1], col_names = T, trim_ws = T)

# check colnames contain 'user', 'sample', 'barcode'
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

# remove white space first
df$sample <- str_replace_all(df$sample, " ", "")
df$sample <- str_replace_na(df$sample)

# checks
if (length(df$barcode) != length(str_unique(df$barcode))) {
  stop('\nBarcodes must be unique!\n')
}

sn_vector <- str_detect(df$sample, '^[a-zA-Z0-9\\_\\-]+$')
if (!all(sn_vector)) {
  stop(
    paste0(
      '\nSample names with special characters:',
      '\n--------------------------------------\n',
      str_flatten(df$sample[!sn_vector], collapse = ", "),
      '\n--------------------------------------\n'
    )
  )
}

# check for duplicate sample names per user
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

write_csv(df, file = paste0(bname, '-checked.csv'))

