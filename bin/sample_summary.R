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

df <- vroom(file = csvfiles) %>% 
  mutate(diff = abs(log2(obs_size / user_size)))
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
      columnDefs = list(list(targets = 'diff', visible = FALSE)), # hide diff
      searchHighlight = TRUE,
      rowCallback = JS(rowCallback),
      autoWidth = TRUE, pageLength = 125,
      dom = 'Btp',
      paging = FALSE,
      buttons = c('copy', 'csv', 'excel')
    )
  ) %>% 
# style user_size based on diff
DT::formatStyle('user_size', 'diff', color = styleInterval(c(0.5, 1), c('#1e8449', '#f5b041', '#e74c3c'))) %>%
DT::formatRound('assembly_quality', 2) %>%
DT::formatStyle('assembly_quality', color = styleInterval(c(25, 35), c('#e74c3c', '#f5b041', '#37474F')))

DT::saveWidget(finaltable, file = '00-sample-status-summary.html', title = "sample-status-summary")
