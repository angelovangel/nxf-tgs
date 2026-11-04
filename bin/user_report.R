#!/usr/bin/env Rscript

# Combine sample_status.csv and mapping_summary.csv for a specific user
# and generate an html report

# arg[1] is user name
# arg[2] is sample_status.csv path
# arg[3] is mapping_summary.csv path (only wf-clone-validation, may not exist)
# arg[4] is pipeline name and version, e.g. 'wf-amplicon v1.2.2' (optional)

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

# difference between expected and assembly size, used to color user_size
if ("assembly_size" %in% names(df) && "user_size" %in% names(df)) {
  df <- df %>% 
    mutate(diff = coalesce(abs(log2(assembly_size / user_size)), -1)) # -1 = no assembly, not colored
}

n_samples <- nrow(df)

# report title and assembly label depend on the epi2me workflow (pipeline_label starts with its name)
report_type <- dplyr::case_when(
  startsWith(pipeline_label, 'wf-amplicon') ~ 'Amplicon',
  startsWith(pipeline_label, 'wf-bacterial-genomes') ~ 'Bacterial Genome',
  TRUE ~ 'Plasmid'
)

# wf-bacterial-genomes: show predicted resistance classes as short badges, wrapping within the cell (full name as tooltip)
amr_abbrev <- c(
  'carbapenem' = 'CARB', 'beta-lactam' = 'BL', 'aminoglycoside' = 'AG', 'aminocyclitol' = 'AMC',
  'macrolide' = 'MAC', 'lincosamide' = 'LIN', 'streptogramin a' = 'SGA', 'streptogramin b' = 'SGB',
  'tetracycline' = 'TET', 'quinolone' = 'QUIN', 'folate pathway antagonist' = 'FPA', 'amphenicol' = 'PHE',
  'fosfomycin' = 'FOS', 'rifamycin' = 'RIF', 'polymyxin' = 'POL', 'glycopeptide' = 'GLY',
  'oxazolidinone' = 'OXZ', 'nitroimidazole' = 'NIT', 'fusidic acid' = 'FUS', 'pseudomonic acid' = 'MUP',
  'quaternary ammonium compound' = 'QAC'
)
abbrev <- function(cls) ifelse(cls %in% names(amr_abbrev), amr_abbrev[cls], toupper(substr(cls, 1, 4)))
# pale background per class, carbapenem is the only solid (red) badge
amr_colors <- c(
  'carbapenem' = 'background:#e74c3c;color:#fff;',
  'beta-lactam' = 'background:#fadbd8;color:#78281f;',
  'aminoglycoside' = 'background:#d6eaf8;color:#1b4f72;', 'aminocyclitol' = 'background:#d6eaf8;color:#1b4f72;',
  'macrolide' = 'background:#e8daef;color:#4a235a;', 'lincosamide' = 'background:#e8daef;color:#4a235a;',
  'streptogramin a' = 'background:#e8daef;color:#4a235a;', 'streptogramin b' = 'background:#e8daef;color:#4a235a;',
  'tetracycline' = 'background:#d5f5e3;color:#145a32;', 'quinolone' = 'background:#fdebd0;color:#7e5109;',
  'folate pathway antagonist' = 'background:#fcf3cf;color:#7d6608;', 'amphenicol' = 'background:#d1f2eb;color:#0e6251;',
  'fosfomycin' = 'background:#f6ddcc;color:#6e2c00;', 'rifamycin' = 'background:#fae5d3;color:#784212;',
  'polymyxin' = 'background:#f5b7b1;color:#641e16;', 'glycopeptide' = 'background:#d4e6f1;color:#154360;',
  'quaternary ammonium compound' = 'background:#e5e7e9;color:#424949;'
)
amr_color <- function(cls) ifelse(cls %in% names(amr_colors), amr_colors[cls], 'background:#e9ecef;color:#2c3e50;')
amr_badges <- function(r) {
  if (is.na(r)) return(NA_character_)
  cls <- strsplit(r, ', ', fixed = TRUE)[[1]]
  cls <- c(intersect('carbapenem', cls), setdiff(cls, 'carbapenem'))
  colors <- amr_color(cls)
  paste0('<span class="amr-tip" data-tip="', cls, '" style="display:inline-block;padding:0 4px;margin:1px;border-radius:3px;font-size:11px;', colors, '">', abbrev(cls), '</span>', collapse = '')
}
amr_legend <- NULL
if ("resistance" %in% names(df) && any(!is.na(df$resistance))) {
  present <- sort(unique(unlist(strsplit(na.omit(df$resistance), ', ', fixed = TRUE))))
  amr_legend <- tags$div(
    style = "font-family: system-ui, -apple-system, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; font-size: 12px; color: #7f8c8d; margin-top: 15px;",
    tags$strong("Predicted resistance: "),
    paste(paste0(abbrev(present), " = ", present), collapse = ", ")
  )
  df$resistance <- vapply(df$resistance, amr_badges, character(1), USE.NAMES = FALSE)
}

# css tooltip for resistance badges, a native title attribute only shows after a browser-controlled delay
amr_style <- tags$style(HTML(paste0(
  ".amr-tip{position:relative;cursor:default}.amr-tip:hover::after{content:attr(data-tip);position:absolute;bottom:100%;right:0;z-index:10;width:max-content;max-width:320px;padding:2px 6px;margin-bottom:2px;border-radius:3px;font-size:12px;text-align:left;white-space:normal;background:#2c3e50;color:#fff;pointer-events:none}",
  ".amr-toggle{border:0;padding:0;background:none;font:inherit;color:inherit;cursor:pointer;white-space:nowrap}.amr-toggle::before{content:'\\25B8';display:inline-block;margin-right:4px;color:#7f8c8d}.amr-toggle[aria-expanded=true]::before{transform:rotate(90deg)}",
  ".amr-details{margin-left:auto;border-collapse:collapse;font-size:12px}.amr-details th,.amr-details td{padding:2px 10px;text-align:left;border-bottom:1px solid #e5e7e9}.amr-details .num{text-align:right}"
)))

rowCallback <- c(
    "function(row, data){",
    "  for(var i=0; i<data.length; i++){",
    "    if(data[i] === null){",
    "      $(this.api().cell(row, i).node()).html('-')", # by data column, td:eq() would count only visible cells
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
  'diff' = 'Diff',
  'n_contigs' = 'Contigs',
  'n_circular' = 'Circular Contigs',
  'largest_contig' = 'Largest Contig',
  'mean_depth' = 'Mean Depth',
  'primary_ratio' = 'Primary Ratio',
  'n_plasmids' = 'Plasmid Contigs',
  'species' = 'Species (ANI)',
  'other_species' = 'Other Species',
  'mlst_st' = 'MLST',
  'resistance' = 'Predicted Resistance',
  'n_amr_genes' = 'AMR Genes'
)

title_defs <- lapply(intersect(names(df), names(col_titles)), function(col) {
  list(targets = which(names(df) == col) - 1, title = col_titles[[col]])
})

# wf-bacterial-genomes: arrow on the AMR gene count opens a second row with the per-gene details
has_gene_details <- all(c("n_amr_genes", "amr_genes") %in% names(df))
amr_details_callback <- if (has_gene_details) c(
  "table.on('click', 'button.amr-toggle', function(){",
  "  var btn = $(this), row = table.row(btn.closest('tr'));",
  "  if(row.child.isShown()){",
  "    row.child.hide();",
  "    btn.attr('aria-expanded', 'false');",
  "    return;",
  "  }",
  "  var cols = [['gene', 'Gene'], ['drugs', 'Associated drugs'], ['identity', 'Identity (%)'], ['coverage', 'Coverage (%)'], ['contigs', 'Contigs']];",
  "  var numeric = function(key){ return key === 'identity' || key === 'coverage'; };",
  "  var details = $('<table class=\"amr-details\"><thead><tr></tr></thead><tbody></tbody></table>');",
  "  cols.forEach(function(c){",
  "    $('<th>').text(c[1]).toggleClass('num', numeric(c[0])).appendTo(details.find('thead tr'));",
  "  });",
  paste0("  JSON.parse(row.data()[", which(names(df) == "amr_genes") - 1, "]).forEach(function(g){"),
  "    var tr = $('<tr>').appendTo(details.find('tbody'));",
  "    cols.forEach(function(c){",
  "      $('<td>').text(g[c[0]] === '' ? '-' : g[c[0]]).toggleClass('num', numeric(c[0])).appendTo(tr);",
  "    });",
  "  });",
  "  row.child(details).show();",
  "  btn.attr('aria-expanded', 'true');",
  "});"
) else "return table;"

finaltable <- 
  DT::datatable(
    df,
    class = c('compact', 'hover', 'row-border'),
    escape = FALSE,
    callback = JS(amr_details_callback),
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
    DT::formatStyle('user_size', 'diff', color = styleInterval(c(-0.5, 0.5, 1), c('inherit', '#228B22', '#f5b041', '#e74c3c')))
}

# hide specific columns
cols_to_hide <- intersect(names(df), c("diff", "validate", "allreads", "percent_nonmapping", "user", "amr_genes"))
# amplicon/bacterial: also hide columns without any values (e.g. obs_size, assembly_quality)
if (report_type != 'Plasmid') {
  empty_cols <- names(df)[sapply(df, function(x) all(is.na(x)))]
  cols_to_hide <- union(cols_to_hide, setdiff(empty_cols, c("sample", "status")))
}
# bacterial genomes: expected size is a rough placeholder in the samplesheet
if (report_type == 'Bacterial Genome') {
  cols_to_hide <- union(cols_to_hide, intersect(names(df), "user_size"))
}
if (length(cols_to_hide) > 0) {
  targets_to_hide <- which(names(df) %in% cols_to_hide) - 1
  finaltable$x$options$columnDefs <- c(
    finaltable$x$options$columnDefs,
    list(list(targets = targets_to_hide, visible = FALSE))
  )
}
# render the arrow only for display so the count still sorts as a number
if (has_gene_details) {
  finaltable$x$options$columnDefs <- c(
    finaltable$x$options$columnDefs,
    list(list(
      targets = which(names(df) == "n_amr_genes") - 1,
      render = JS(paste0(
        "function(data, type, row){",
        "  if(type !== 'display' || data === null || row[", which(names(df) == "amr_genes") - 1, "] === null) return data;",
        "  return '<button type=\"button\" class=\"amr-toggle\" aria-expanded=\"false\" title=\"Show AMR genes\">' + data + '</button>';",
        "}"
      ))
    ))
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
if ("primary_ratio" %in% names(df)) {
  finaltable <- finaltable %>%
    DT::formatRound('primary_ratio', 3)
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
    tags$h1(paste("ONT", report_type, "Assembly Report"), style = "margin-top: 0; margin-bottom: 4px; color: #2c3e50; font-size: 24px; font-weight: 600; letter-spacing: -0.5px;"),
    tags$div("For research use only. Not for use in diagnostic procedures.", style = "color: #e74c3c; font-size: 13px; font-style: italic;"),
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
        tags$strong(paste0(report_type, " Assembly: "), style = "color: #7f8c8d; font-weight: 500;"),
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
if (!is.null(amr_legend)) report <- htmlwidgets::appendContent(report, amr_legend)
if (report_type == 'Bacterial Genome') report <- htmlwidgets::appendContent(report, amr_style)

out_file <- paste0("02-", user_name, "-assembly-report.html")
htmlwidgets::saveWidget(report, file = out_file, title = paste(user_name, "Report"))
