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
df <- vroom(countsfiles, id = 'sample') %>% 
  mutate(
    user = arg[2],
    sample = str_remove(basename(sample), "-mapping-counts.txt")) %>%
  mutate(
    percent_assembly = round(assembly/allreads*100, 2), 
    percent_coli = round(NC_000913.3/allreads*100, 2), 
    percent_nonmapping = round((allreads - (assembly+NC_000913.3))/allreads*100, 2)
    )
  
df %>%
  dplyr::relocate('user') %>%
  write.csv(file = 'mapping-summary.csv', sep = ",", row.names = F, col.names = T)
