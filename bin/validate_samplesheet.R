#!/usr/bin/env Rscript

# reads csv/tsv and checks if valid
# writes back *-checked.csv if ok

arg <- commandArgs(trailingOnly = T)

library(tools)
library(readr)
library(stringr)

bname <- tools::file_path_sans_ext(basename(arg[1]))

df <- readr::read_delim(arg[1])
# remove white space first
df$sample <- str_replace_all(df$sample, " ", "")

# checks
if (length(df$barcode) != length(str_unique(df$barcode))) {
  stop('Barcodes must be unique!')
}

sn_vector <- str_replace_na(df$sample) %>% str_detect('^[a-zA-Z0-9\\_\\-]+$')
if (!all(sn_vector)) {
  stop(
    paste0(
      'Sample names with special characters: ', 
      str_flatten(df$sample[!sn_vector], collapse = ", ")
    )
  )
}

write_csv(df, file = paste0(bname, '-checked.csv'))

