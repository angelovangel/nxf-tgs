#!/usr/bin/env Rscript

# combine mapping_counts.txt for all samples of one user
# generate one summary file per user

# arg[1] is a pattern to list files
# arg[2] is user
library(vroom)
library(dplyr)
library(stringr)

arg <- commandArgs(trailingOnly = T)

countsfiles <- list.files(pattern = arg[1], full.names = T)

# do not use vroom directly for it fails if the files have different number of columns, for example when no reads map to coli

df <- 
  lapply(countsfiles, vroom, id = 'sample', show_col_types = F) %>% 
  bind_rows() %>%
  #vroom(countsfiles, id = 'sample') %>% 
  mutate(
    user = arg[2],
    sample = str_remove(basename(sample), "-mapping-counts.txt")) %>%
  mutate(
    percent_assembly = assembly/allreads*100, 
    percent_coli = NC_000913.3/allreads*100, 
    percent_nonmapping = (allreads - (assembly+NC_000913.3))/allreads*100
    ) %>%
  dplyr::select(-c('assembly', 'NC_000913.3'))
  
df %>%
  dplyr::relocate('user') %>%
  write.csv(file = 'mapping-summary.csv', sep = ",", row.names = F, col.names = T)
