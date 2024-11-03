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

# format NA as NA in DT
# https://stackoverflow.com/questions/58526047/customizing-how-datatables-displays-missing-values-in-shiny
rowCallback <- c(
    "function(row, data){",
    "  for(var i=0; i<data.length; i++){",
    "    if(data[i] === null){",
    "      $('td:eq('+i+')', row).html('-')",
    "        .css({'color': 'rgb(151,151,151)', 'font-style': 'italic'});",
    "    }",
    "  }",
    "}"  
  )

finaltable <- 
DT::datatable(
  dplyr::arrange(df, user, sample),
  class = 'compact',
  caption = paste0("Run name: ", arg[2], " | Time: ", format.POSIXct(Sys.time())),
  # style = 'bootstrap',
  escape = F, filter = 'top',
  extensions = 'Buttons', rownames = FALSE,
  options = list(
    searchHighlight = TRUE,
    rowCallback = JS(rowCallback),
    autoWidth = TRUE, pageLength = 25,
    dom = 'Btp',
    paging = FALSE,
    buttons = c('copy', 'csv', 'excel')
  )
) %>% 
DT::formatRound('assembly_quality', 2) %>%
DT::formatStyle('assembly_quality', color = styleInterval(c(30, 40), c('red', 'orange', 'black')))

DT::saveWidget(finaltable, '00-sample-status-summary.html')
