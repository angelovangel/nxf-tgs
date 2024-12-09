#!/usr/bin/env Rscript

# combine all sample_status.csv from all users and 
# generate an html table (with user search field) 

# arg[1] is a pattern to list files
# arg[2] is workflow id
library(vroom)
library(dplyr)
library(DT)

arg <- commandArgs(trailingOnly = T)
csvfiles <- list.files(pattern = arg[1], full.names = T)

df <- vroom(file = csvfiles)

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
# lexocographical arrange of file
locale <- list(locale = "en_US", numeric = TRUE)

finaltable <- 
  DT::datatable(
    dplyr::arrange(df, user, stringi::stri_rank(sample, opts_collator = locale)),
    class = 'compact',
    # caption = paste0("Run name: ", arg[2], " | Time: ", format.POSIXct(Sys.time())),
    caption = htmltools::tags$caption(
      style = 'caption-side: bottom; text-align: left; color: grey;',
      paste0(arg[2], " | ", format.POSIXct(Sys.time()))
    ),
    # style = 'bootstrap',
    escape = F, filter = 'top',
    extensions = 'Buttons', rownames = FALSE,
    options = list(
      searchHighlight = TRUE,
      rowCallback = JS(rowCallback),
      autoWidth = TRUE, pageLength = 125,
      dom = 'Btp',
      paging = FALSE,
      buttons = c('copy', 'csv', 'excel')
    )
  ) %>% 
  # style user_size based on diff
  #DT::formatStyle('user_size', 'diff', color = styleInterval(c(0.5, 1), c('#1e8449', '#f5b041', '#e74c3c'))) %>%
  DT::formatRound(c('percent_assembly', 'percent_coli', 'percent_nonmapping'), digits = 2) %>%
  DT::formatStyle('percent_assembly', color = styleInterval(c(40, 60), c('#e74c3c', '#f5b041', '#1e8449'))) %>%
  DT::formatStyle(c('percent_coli', 'percent_nonmapping'), color = styleInterval(c(20, 40), c('#1e8449', '#f5b041', '#e74c3c')))

DT::saveWidget(finaltable, '00-sample-mapping-summary.html')