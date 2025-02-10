#!/usr/bin/env Rscript

# combine sample_status.txt and samplesheet_validated.csv
# produce sample_status.csv summary for every users

# arg[1] is user 
# arg[2] is sample_status.txt
# arg[3] is samplesheet_validated
# 
library(dplyr)
library(vroom)

arg <- commandArgs(trailingOnly = T)
tsvfiles <- list.files(path = ".", pattern = "*.assembly_stats.tsv", full.names = T)

df1 <- vroom(
  arg[2], 
  col_names = c('sample', 'pass_fail', 'length'), 
  trim_ws = T, na = 'N/A'
  ) #! sample_status.txt

#samplesheet_validated
df2 <- vroom(
  arg[3], 
  col_names = T, 
  trim_ws = T, delim = ","
  ) %>% 
  dplyr::filter(user == arg[1])

if (length(tsvfiles > 0)) {
  df3 <- vroom(
    tsvfiles, 
    col_select = (c('sample_name', assembly_quality = 'mean_quality')), 
    show_col_types = F
  )
} else {
  df3 <- data.frame(
    sample_name = NA,
    assembly_quality = NA
  )
}

df2 %>% 
  #mutate(user = arg[1]) %>%
  left_join(df1, by = join_by(sample)) %>%
  dplyr::select(c('user', 'sample', 'validate','barcode', user_size = 'dna_size', 'obs_size', assembly_size = 'length')) %>%
  left_join(df3, by = c('sample' = 'sample_name')) %>%
  write.csv(file = 'sample-status.csv', row.names = F, col.names = T, sep = ",")
  