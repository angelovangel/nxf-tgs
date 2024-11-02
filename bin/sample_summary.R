#!/usr/bin/env Rscript

# combine all sample_status.csv from all users and 
# generate an html table (with user search field) 

# arg[1] is a pattern to list files
# arg[2] is workflow id
library(vroom)
library(dplyr)
library(DT)
#library(reactable)
#library(htmlwidgets)

arg <- commandArgs(trailingOnly = T)
csvfiles <- list.files(pattern = arg[1], full.names = T)

df <- vroom(file = csvfiles) 
#df
# finaltable <- 
#   reactable(
#     dplyr::arrange(df, user, sample), 
#     filterable = TRUE, minRows = 25
#   )
finaltable <- 
DT::datatable(
  dplyr::arrange(df, user, sample),
  class = 'compact',
  caption = paste0("Run name:", arg[2], " | ", format.POSIXct(Sys.time())),
  # style = 'bootstrap',
  escape = F, filter = 'top',
  extensions = 'Buttons', rownames = FALSE,
  options = list(
    autoWidth = TRUE, pageLength = 25,
    dom = 'Btp',
    paging = FALSE,
    buttons = c('copy', 'csv', 'excel')
  )
) %>% 
DT::formatStyle('assembly_quality', color = styleInterval(c(30, 40), c('red', 'orange', 'black')))

DT::saveWidget(finaltable, '00-sample-status-summary.html')
