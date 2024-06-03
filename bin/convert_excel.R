#!/usr/bin/env Rscript

# save Excel as csv reliably
# pass excel file as positional param
# write csv in same path as Excel
arg <- commandArgs(trailingOnly = T)

library(readxl)
library(tools)

if (length(arg) != 1) {
  stop('One argument needed (path to excel file)')
}

#bdir <- dirname(tools::file_path_as_absolute(arg[1]))
bname <- tools::file_path_sans_ext(basename(arg[1]))


df <- readxl::read_excel(arg[1])
write.csv(df, file = paste0(bname, '.csv'), row.names = F, quote = FALSE)