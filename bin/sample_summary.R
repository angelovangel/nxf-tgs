#!/usr/bin/env Rscript

# combine all sample_status.csv from all users and 
# generate an html table (with user search field) 

# arg is a list of csv files to process
library(vroom)
library(dplyr)
library(DT)
library(htmlwidgets)

arg <- commandArgs(trailingOnly = T)
csvfiles <- list.files(pattern = arg, full.names = T)

df <- vroom(file = csvfiles) 

finaltable <- 
  DT::datatable(
    dplyr::arrange(df, user, sample),
    class = 'compact',
    # style = 'bootstrap',
    escape = F, filter = 'top',
    extensions = 'Buttons', rownames = FALSE,
    options = list(
      autoWidth = TRUE, pageLength = 25,
      dom = 'Btp',
      paging = FALSE,
      buttons = c('copy', 'csv', 'excel')
    )
  )

htmlwidgets::saveWidget(finaltable, '00-sample-status-summary.html')
