#!/usr/bin/env Rscript

# combine wf-bacterial-genomes results.json and samplesheet_validated.csv
# produce sample_status.csv summary for every user (same base columns as sample_status.R)

# arg[1] is user
# arg[2] is results.json (from wf-bacterial-genomes out_dir)
# arg[3] is samplesheet_validated
#
# samples missing from results.json, with sample_pass false or without contigs are 'fail'
library(dplyr)
library(vroom)
library(jsonlite)

arg <- commandArgs(trailingOnly = T)

res <- fromJSON(arg[2], simplifyVector = F)

# NULL/empty safe accessor
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

parse_sample <- function(s) {
  r <- s$results
  contigs <- r$assembly$contig %||% list()
  flye <- r$assembly$flye_stats %||% list()
  lengths <- vapply(contigs, function(x) as.numeric(x$length %||% NA_real_), numeric(1))
  depths <- vapply(contigs, function(x) as.numeric(x$coverage$mean %||% NA_real_), numeric(1))
  n_contigs <- length(contigs)
  matches <- r$species_identification$detected_matches %||% list()
  data.frame(
    sample = as.character(s$alias),
    n_seqs = r$fastq$n_seqs %||% NA_real_,
    pass_fail = if (isTRUE(s$sample_pass) && n_contigs > 0) 'pass' else 'fail',
    assembly_size = if (n_contigs > 0) sum(lengths) else NA_real_,
    n_contigs = n_contigs,
    n_circular = sum(vapply(flye, function(x) identical(x$circular, 'Y'), logical(1))),
    largest_contig = if (n_contigs > 0) max(lengths) else NA_real_,
    # length-weighted mean depth over all contigs
    mean_depth = if (n_contigs > 0) round(sum(lengths * depths) / sum(lengths)) else NA_real_,
    n_plasmids = r$plasmid_identification$n_plasmid_contigs %||% NA_integer_,
    species = if (length(matches) > 0) matches[[1]]$species_call %||% NA_character_ else NA_character_,
    mlst_st = r$sequence_typing$sequence_type %||% NA_character_,
    n_amr_genes = length(r$antimicrobial_resistance$seq_regions %||% list()),
    stringsAsFactors = F
  )
}

extra_cols <- c('n_contigs', 'n_circular', 'largest_contig', 'mean_depth', 'n_plasmids', 'species', 'mlst_st', 'n_amr_genes')

df1 <- bind_rows(lapply(res$samples, parse_sample))
if (nrow(df1) == 0) {
  df1 <- data.frame(sample = character(), n_seqs = numeric(), pass_fail = character(), assembly_size = numeric())
  for (col in extra_cols) df1[[col]] <- logical()
}

#samplesheet_validated
df2 <- vroom(
  arg[3],
  col_names = T,
  trim_ws = T, delim = ","
  ) %>%
  dplyr::filter(user == arg[1]) %>%
  dplyr::mutate(sample = as.character(sample))

# get_maxbin.sh (obs_size, nreads) is not run for wf-bacterial-genomes, add the missing columns
for (col in c('dna_size', 'obs_size', 'nreads')) {
  if (!col %in% names(df2)) df2[[col]] <- NA_real_
}

df2 %>%
  left_join(df1, by = join_by(sample)) %>%
  dplyr::mutate(
    pass_fail = coalesce(pass_fail, 'fail'),
    nreads = coalesce(nreads, n_seqs), # reads after wf-bacterial-genomes length filtering
    assembly_quality = NA_real_ # no per-sample assembly quality in wf-bacterial-genomes
  ) %>%
  dplyr::select(
    c('user', 'sample', 'barcode', 'validate', status = 'pass_fail', 'nreads', user_size = 'dna_size', 'obs_size', 'assembly_size', 'assembly_quality'),
    all_of(extra_cols)
  ) %>%
  write.csv(file = 'sample-status.csv', row.names = F)
