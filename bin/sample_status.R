#!/usr/bin/env Rscript

# combine sample_status.txt and samplesheet_validated.csv
# produse sample_status.csv summary for every users

# arg[1] is user arg[2] is sample_status, arg[3] is samplesheet_validated
library(dplyr)
library(readr)

arg <- commandArgs(trailingOnly = T)

df1 <- readr::read_delim(arg[2], col_names = T, trim_ws = T)
colnames(df1) <- c('sample', 'pass_fail', 'length')

df2 <- readr::read_delim(arg[3], col_names = T, trim_ws = T)

df1 %>% 
  mutate(user = arg[1]) %>%
  left_join(df2, by = join_by(user, sample)) %>%
  dplyr::select(c('user', 'sample', 'barcode', user_size = 'dna_size', 'obs_size', assembly_size = 'length')) %>%
  write.csv(file = 'sample-status.csv', row.names = F, col.names = T, sep = ",")
  