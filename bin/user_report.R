#!/usr/bin/env Rscript

# Combine sample_status.csv and mapping_summary.csv for a specific user
# and generate an html report

# arg[1] is user name
# arg[2] is sample_status.csv path
# arg[3] is mapping_summary.csv path
# arg[4] is pipeline version or other info (optional)

library(vroom)
library(dplyr)
library(DT)
library(htmltools)
library(htmlwidgets)

args <- commandArgs(trailingOnly = TRUE)
user_name <- args[1]
sample_status_file <- args[2]
mapping_summary_file <- args[3]

if (length(args) >= 4) {
  pipeline_label <- args[4]
} else {
  pipeline_label <- ""
}

if (length(args) >= 5) {
  nxf_tgs_version <- args[5]
} else {
  nxf_tgs_version <- ""
}

# Read sample status
df_status <- vroom(sample_status_file, delim = ",", show_col_types = FALSE)

# Read mapping summary if it exists and is not empty
if (file.exists(mapping_summary_file)) {
  df_mapping <- vroom(mapping_summary_file, delim = ",", show_col_types = FALSE)
  
  # join by sample and user
  if (nrow(df_mapping) > 0 && "sample" %in% names(df_mapping) && "sample" %in% names(df_status)) {
    df <- df_status %>%
      left_join(df_mapping, by = c("user", "sample"))
  } else {
    df <- df_status
  }
} else {
  df <- df_status
}

# calculate difference or just format it as in sample_summary.R
if ("obs_size" %in% names(df) && "user_size" %in% names(df)) {
  df <- df %>% 
    mutate(diff = abs(log2(obs_size / user_size)))
}

n_samples <- nrow(df)

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

col_titles <- c(
  'user' = 'User',
  'sample' = 'Sample',
  'barcode' = 'Barcode',
  'validate' = 'Validate',
  'status' = 'Status',
  'nreads' = 'Reads',
  'user_size' = 'Expected Size',
  'obs_size' = 'Observed Size',
  'assembly_size' = 'Assembly Size',
  'assembly_quality' = 'Assembly Quality',
  'allreads' = 'Total Reads',
  'percent_assembly' = 'Assembly Mapped (%)',
  'percent_coli' = 'E. coli Mapped (%)',
  'percent_nonmapping' = 'Unmapped (%)',
  'diff' = 'Diff'
)

title_defs <- lapply(intersect(names(df), names(col_titles)), function(col) {
  list(targets = which(names(df) == col) - 1, title = col_titles[[col]])
})

finaltable <- 
  DT::datatable(
    df,
    class = c('compact', 'hover'),
    escape = FALSE,
    #extensions = c('Select', 'SearchPanes', 'Buttons'), 
    rownames = FALSE,
    options = list(
      searchHighlight = TRUE,
      rowCallback = JS(rowCallback),
      autoWidth = TRUE, pageLength = 125,
      dom = 'Brt',
      paging = FALSE,
      buttons = c('copy', 'colvis'),
      searchPanes = list(show = FALSE),
      columnDefs = title_defs
    )
  )

# apply conditional formatting if columns exist
if ("diff" %in% names(df)) {
  finaltable <- finaltable %>% 
    DT::formatStyle('user_size', 'diff', color = styleInterval(c(0.5, 1), c('#228B22', '#f5b041', '#e74c3c')))
}

# hide specific columns
cols_to_hide <- intersect(names(df), c("diff", "validate", "allreads", "percent_nonmapping", "user"))
if (length(cols_to_hide) > 0) {
  targets_to_hide <- which(names(df) %in% cols_to_hide) - 1
  finaltable$x$options$columnDefs <- c(
    finaltable$x$options$columnDefs,
    list(list(targets = targets_to_hide, visible = FALSE))
  )
}
if ("nreads" %in% names(df)) {
  finaltable <- finaltable %>%
    DT::formatStyle('nreads', color = styleInterval(c(200, 500), c('#e74c3c', '#f5b041', 'inherit')))
}
if ("status" %in% names(df)) {
  finaltable <- finaltable %>%
    DT::formatStyle('status', color = styleEqual(c('pass', 'fail'), c('#228B22', '#e74c3c')))
}
if ("assembly_quality" %in% names(df)) {
  finaltable <- finaltable %>%
    DT::formatRound('assembly_quality', 0) %>%
    DT::formatStyle('assembly_quality', color = styleInterval(c(25, 35), c('#e74c3c', '#f5b041', '#228B22')))
}
if ("percent_assembly" %in% names(df)) {
  finaltable <- finaltable %>%
    DT::formatRound('percent_assembly', 1)
}
if ("percent_coli" %in% names(df)) {
  finaltable <- finaltable %>%
    DT::formatRound('percent_coli', 1)
}
if ("percent_nonmapping" %in% names(df)) {
  finaltable <- finaltable %>%
    DT::formatRound('percent_nonmapping', 1)
}

# Create a modern header and combine with table
header <- tags$div(
  style = "font-family: system-ui, -apple-system, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; padding: 25px; margin-bottom: 25px; background-color: #f8f9fa; border-radius: 8px; border-left: 6px solid #2c3e50; box-shadow: 0 2px 4px rgba(0,0,0,0.05); display: flex; justify-content: space-between; align-items: center;",
  tags$div(
    tags$h1("ONT Plasmid Assembly Report", style = "margin-top: 0; color: #2c3e50; font-size: 24px; font-weight: 600; letter-spacing: -0.5px;"),
    tags$div(
      style = "display: flex; gap: 30px; margin-top: 15px; font-size: 15px;",
      tags$div(
        tags$strong("User: ", style = "color: #7f8c8d; font-weight: 500;"), 
        tags$span(user_name, style = "color: #2c3e50; font-weight: 600;")
      ),
      tags$div(
        tags$strong("Samples: ", style = "color: #7f8c8d; font-weight: 500;"), 
        tags$span(n_samples, style = "color: #2c3e50; font-weight: 600;")
      ),
      tags$div(
        tags$strong("Plasmid Assembly: ", style = "color: #7f8c8d; font-weight: 500;"),
        tags$span(pipeline_label, style = "color: #2c3e50; font-weight: 600;")
      ),
      tags$div(
        tags$strong("NXF-TGS pipeline: ", style = "color: #7f8c8d; font-weight: 500;"),
        tags$span(nxf_tgs_version, style = "color: #2c3e50; font-weight: 600;")
      )
    )
  ),
  tags$div(
    style = "text-align: right; color: #7f8c8d; font-size: 13px;",
    tags$div("KAUST Bioscience Core Labs"),
    tags$div(format.POSIXct(Sys.time()))
  )
)

report <- htmlwidgets::prependContent(finaltable, header)

out_file <- paste0("00-", user_name, "-assembly-report.html")
htmlwidgets::saveWidget(report, file = out_file, title = paste(user_name, "Report"))
